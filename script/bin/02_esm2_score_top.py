#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import math
from collections import OrderedDict

import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM


# =========================
# 1. 路径设置
# =========================

# 输入候选序列 fasta
INPUT_FASTA = "3.candidates/01_proteinmpnn/sfGFP.fa"

# ESM-2 650M 模型路径
ESM2_DIR = "2.model/esm2_t33_650M_UR50D"

# 输出目录
OUTDIR = "3.candidates/02_esm2_t33_score/"

# 输出前 8 条
TOPN = 8

# 每次同时 mask 多少个位点
# 12G 显卡建议先用 8
# 如果显存不够，改成 4、2、1
# 如果显存很宽裕，可以试 16
MASK_BATCH_SIZE = 8


# =========================
# 2. 读取 fasta
# =========================

def read_fasta(path):
    """
    读取 fasta 文件。

    返回：
        OrderedDict:
        {
            "candidate_00001": "MXXX...",
            "candidate_00002": "MYYY..."
        }
    """

    records = OrderedDict()
    name = None
    seqs = []

    if not os.path.exists(path):
        raise FileNotFoundError(f"Input fasta not found: {path}")

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seqs).replace(" ", "").replace("\t", "").upper()

                header = line[1:].strip()

                # 如果 header 是 candidate_00001|mut=xxx 这种形式，只取前面的 candidate_00001
                name = header.split()[0].split("|")[0]

                seqs = []

            else:
                seqs.append(line)

        if name is not None:
            records[name] = "".join(seqs).replace(" ", "").replace("\t", "").upper()

    return records


# =========================
# 3. 写出 fasta
# =========================

def write_fasta(rows, out_fa):
    """
    写出 Top8 fasta。
    """

    with open(out_fa, "w") as f:
        for row in rows:
            # Top8_candidates.fa 的序列名只保留名称
            # 例如：>sample3|len=236|esm_mean_logp=...  改为  >sample3
            header = row["candidate_id"]

            seq = row["sequence"]

            f.write(f">{header}\n")

            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


# =========================
# 4. ESM-2 单条序列打分
# =========================

def esm2_score_sequence(seq, tokenizer, model, device):
    """
    对一条蛋白序列计算 ESM-2 分数。

    计算逻辑：
        对序列中每个位置依次 mask。
        例如序列 MSKGE：

        [MASK]SKGE  -> 计算真实氨基酸 M 的概率
        M[MASK]KGE  -> 计算真实氨基酸 S 的概率
        MS[MASK]GE  -> 计算真实氨基酸 K 的概率
        MSK[MASK]E  -> 计算真实氨基酸 G 的概率
        MSKG[MASK]  -> 计算真实氨基酸 E 的概率

    最后：
        esm_mean_logp = 所有位置 log 概率的平均值
        esm_pseudo_perplexity = exp(-esm_mean_logp)

    返回：
        esm_mean_logp
        esm_pseudo_perplexity
    """

    inputs = tokenizer(seq, return_tensors="pt")

    input_ids = inputs["input_ids"].to(device)
    attention_mask = inputs["attention_mask"].to(device)

    seq_len = len(seq)

    # 每个真实氨基酸对应的 token id
    aa_token_ids = [tokenizer.convert_tokens_to_ids(aa) for aa in seq]

    total_logp = 0.0

    # 分批 mask，减少显存压力
    for start in range(0, seq_len, MASK_BATCH_SIZE):
        end = min(start + MASK_BATCH_SIZE, seq_len)
        batch_size = end - start

        # 复制 batch_size 份输入
        batch_input_ids = input_ids.repeat(batch_size, 1)
        batch_attention_mask = attention_mask.repeat(batch_size, 1)

        # 每一行 mask 一个不同的位置
        for row_idx, pos0 in enumerate(range(start, end)):
            # ESM tokenizer 会在序列前面加一个特殊 token
            # 所以蛋白第 1 位，对应 token 位置 1
            token_pos = pos0 + 1
            batch_input_ids[row_idx, token_pos] = tokenizer.mask_token_id

        # ==============================
        # 这里是真正调用 ESM-2 的地方
        # ==============================
        with torch.no_grad():
            outputs = model(
                input_ids=batch_input_ids,
                attention_mask=batch_attention_mask
            )
        # ==============================

        logits = outputs.logits

        # logits 转成 log probability
        log_probs = torch.log_softmax(logits, dim=-1)

        # 取每个 mask 位置真实氨基酸的 log probability
        for row_idx, pos0 in enumerate(range(start, end)):
            token_pos = pos0 + 1
            true_aa_id = aa_token_ids[pos0]
            total_logp += log_probs[row_idx, token_pos, true_aa_id].item()

    esm_mean_logp = total_logp / seq_len
    esm_pseudo_perplexity = math.exp(-esm_mean_logp)

    return esm_mean_logp, esm_pseudo_perplexity


# =========================
# 5. 主程序
# =========================

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"

    print("Device:", device)
    print("Input fasta:", INPUT_FASTA)
    print("ESM-2 model:", ESM2_DIR)
    print("MASK_BATCH_SIZE:", MASK_BATCH_SIZE)

    if device == "cpu":
        print("Warning: 当前没有检测到 GPU，会用 CPU 跑，速度会很慢。")

    # =========================
    # 加载 ESM-2
    # =========================

    tokenizer = AutoTokenizer.from_pretrained(ESM2_DIR)

    model = AutoModelForMaskedLM.from_pretrained(ESM2_DIR)

    # 12G 显卡使用 fp16，降低显存占用
    if device == "cuda":
        model = model.half()

    model = model.to(device)
    model.eval()

    # 读取候选序列
    records = read_fasta(INPUT_FASTA)

    if len(records) == 0:
        raise ValueError("No sequences found in final_candidates.fa")

    print("Candidate count:", len(records))

    rows = []

    # 逐条候选序列打分
    for i, (candidate_id, seq) in enumerate(records.items(), start=1):
        print(f"Scoring {i}/{len(records)}: {candidate_id}, length={len(seq)}")

        esm_mean_logp, esm_pseudo_perplexity = esm2_score_sequence(
            seq=seq,
            tokenizer=tokenizer,
            model=model,
            device=device
        )

        rows.append({
            "candidate_id": candidate_id,
            "length": len(seq),
            "esm_mean_logp": esm_mean_logp,
            "esm_pseudo_perplexity": esm_pseudo_perplexity,
            "sequence": seq
        })

    # 按 esm_mean_logp 从大到小排序
    rows.sort(
        key=lambda x: x["esm_mean_logp"],
        reverse=True
    )

    # 输出所有候选序列打分
    all_score_tsv = os.path.join(OUTDIR, "All_candidates_ESM2_score.tsv")

    with open(all_score_tsv, "w") as f:
        f.write("rank\tcandidate_id\tlength\tesm_mean_logp\tesm_pseudo_perplexity\tsequence\n")

        for rank, row in enumerate(rows, start=1):
            f.write(
                f"{rank}\t"
                f"{row['candidate_id']}\t"
                f"{row['length']}\t"
                f"{row['esm_mean_logp']:.8f}\t"
                f"{row['esm_pseudo_perplexity']:.8f}\t"
                f"{row['sequence']}\n"
            )

    # 输出 Top 候选表格
    top_rows = rows[:TOPN]
    top_tsv = os.path.join(OUTDIR, "Top_candidates.tsv")

    with open(top_tsv, "w") as f:
        f.write("rank\tcandidate_id\tlength\tesm_mean_logp\tesm_pseudo_perplexity\tsequence\n")

        for rank, row in enumerate(top_rows, start=1):
            f.write(
                f"{rank}\t"
                f"{row['candidate_id']}\t"
                f"{row['length']}\t"
                f"{row['esm_mean_logp']:.8f}\t"
                f"{row['esm_pseudo_perplexity']:.8f}\t"
                f"{row['sequence']}\n"
            )

    # 输出 Top 候选 fasta
    top_fa = os.path.join(OUTDIR, "Top_candidates.fa")
    write_fasta(top_rows, top_fa)

    print("Done.")
    print("All score table:", all_score_tsv)
    print("Top table:", top_tsv)
    print("Top fasta:", top_fa)


if __name__ == "__main__":
    main()