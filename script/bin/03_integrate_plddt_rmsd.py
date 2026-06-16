#!/usr/bin/env python3
import os
import glob
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser, Superimposer

# ============================================================
# 输入文件
# ============================================================

REF_PDB = "1.pdb/sfGFP.pdb"
WT_FA = "1.pdb/sfGFP.fa"
CAND_FA = "3.candidates/02_esm2_t33_score/Top_candidates.fa"
COLABFOLD_DIR = "4.filter/01_colabfol"

# ============================================================
# 输出目录
# ============================================================

OUT_DIR = "4.filter/02_out"
OUT_ALL = f"{OUT_DIR}/integrated_plddt_rmsd_all.tsv"
OUT_PASS = f"{OUT_DIR}/integrated_plddt_rmsd_pass.tsv"
OUT_PASS_FA = f"{OUT_DIR}/pass_candidates.fa"
OUT_PASS_PDB_DIR = f"{OUT_DIR}/pass_candidates_pdb"

# ============================================================
# 筛选条件
# ============================================================

PLDDT_CUTOFF = 90.0
RMSD_CUTOFF = 2

# 如果只想保留 Top10，改成 10；不限制就用 None
TOP_N = None


def mkdir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def sample_sort_key(sample_id):
    if sample_id.startswith("sample"):
        n = sample_id.replace("sample", "")
        if n.isdigit():
            return int(n)
    return sample_id


def fmt(x):
    if x is None:
        return "NA"
    return f"{x:.4f}"


def get_ca_atoms(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_file)

    ca_atoms = []

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] == " " and "CA" in res:
                    ca_atoms.append(res["CA"])
            break
        break

    return ca_atoms


def mean_plddt_from_pdb(pdb_file):
    ca_atoms = get_ca_atoms(pdb_file)

    if len(ca_atoms) == 0:
        return None

    vals = [atom.get_bfactor() for atom in ca_atoms]
    return sum(vals) / len(vals)


def calc_ca_rmsd(ref_pdb, query_pdb):
    ref_atoms = get_ca_atoms(ref_pdb)
    query_atoms = get_ca_atoms(query_pdb)

    ref_len = len(ref_atoms)
    query_len = len(query_atoms)

    n = min(ref_len, query_len)

    if n == 0:
        return None, ref_len, query_len

    ref_use = ref_atoms[:n]
    query_use = query_atoms[:n]

    sup = Superimposer()
    sup.set_atoms(ref_use, query_use)

    return sup.rms, ref_len, query_len


def find_rank1_pdb(sample_id):
    pattern = os.path.join(COLABFOLD_DIR, f"{sample_id}_unrelaxed_rank_001_*.pdb")
    pdbs = glob.glob(pattern)

    if len(pdbs) == 0:
        return None

    return pdbs[0]


def count_mutations(wt_seq, seq):
    if len(wt_seq) != len(seq):
        return "NA"

    n = 0
    for a, b in zip(wt_seq, seq):
        if a != b:
            n += 1

    return n


# ============================================================
# 主程序
# ============================================================

mkdir(OUT_DIR)

# 每次重建通过筛选 PDB 目录，避免旧文件残留
if os.path.exists(OUT_PASS_PDB_DIR):
    shutil.rmtree(OUT_PASS_PDB_DIR)
mkdir(OUT_PASS_PDB_DIR)

for f in [REF_PDB, WT_FA, CAND_FA]:
    if not os.path.exists(f):
        raise FileNotFoundError(f"找不到文件: {f}")

if not os.path.exists(COLABFOLD_DIR):
    raise FileNotFoundError(f"找不到 ColabFold 目录: {COLABFOLD_DIR}")

wt_seq = str(next(SeqIO.parse(WT_FA, "fasta")).seq).strip().upper()

seq_dict = {}
for record in SeqIO.parse(CAND_FA, "fasta"):
    seq_dict[record.id] = str(record.seq).strip().upper()

sample_ids = sorted(seq_dict.keys(), key=sample_sort_key)

results = []

for sample_id in sample_ids:
    seq = seq_dict[sample_id]
    selected_pdb = find_rank1_pdb(sample_id)

    if selected_pdb is None:
        print(f"[警告] {sample_id} 没有找到 rank_001 PDB")
        mean_plddt = None
        rmsd = None
        ref_ca_len = "NA"
        query_ca_len = "NA"
        selected_pdb = "NA"
    else:
        mean_plddt = mean_plddt_from_pdb(selected_pdb)
        rmsd, ref_ca_len, query_ca_len = calc_ca_rmsd(REF_PDB, selected_pdb)

    mutation_count = count_mutations(wt_seq, seq)

    pass_plddt = mean_plddt is not None and mean_plddt >= PLDDT_CUTOFF
    pass_rmsd = rmsd is not None and rmsd <= RMSD_CUTOFF

    final_status = "PASS" if pass_plddt and pass_rmsd else "FAIL"

    fail_reason = []

    if not pass_plddt:
        fail_reason.append("low_pLDDT")

    if not pass_rmsd:
        fail_reason.append("high_RMSD")

    if len(fail_reason) == 0:
        fail_reason.append("PASS")

    results.append({
        "sample_id": sample_id,
        "mean_pLDDT": mean_plddt,
        "RMSD_CA": rmsd,
        "mutation_count": mutation_count,
        "ref_CA_len": ref_ca_len,
        "query_CA_len": query_ca_len,
        "pass_pLDDT": "YES" if pass_plddt else "NO",
        "pass_RMSD": "YES" if pass_rmsd else "NO",
        "final_status": final_status,
        "fail_reason": ";".join(fail_reason),
        "selected_pdb": selected_pdb,
        "seq": seq,
    })


def sort_key(r):
    status_rank = 0 if r["final_status"] == "PASS" else 1
    plddt_sort = -(r["mean_pLDDT"] if r["mean_pLDDT"] is not None else -999)
    rmsd_sort = r["RMSD_CA"] if r["RMSD_CA"] is not None else 999

    return (
        status_rank,
        plddt_sort,
        rmsd_sort,
        sample_sort_key(r["sample_id"])
    )


results_sorted = sorted(results, key=sort_key)

pass_results = [r for r in results_sorted if r["final_status"] == "PASS"]

if TOP_N is not None:
    pass_results = pass_results[:TOP_N]


header = [
    "sample_id",
    "mean_pLDDT",
    "RMSD_CA",
    "mutation_count",
    "ref_CA_len",
    "query_CA_len",
    "pass_pLDDT",
    "pass_RMSD",
    "final_status",
    "fail_reason",
    "selected_pdb",
]


def write_table(outfile, rows):
    with open(outfile, "w") as out:
        out.write("\t".join(header) + "\n")

        for r in rows:
            row = [
                r["sample_id"],
                fmt(r["mean_pLDDT"]),
                fmt(r["RMSD_CA"]),
                str(r["mutation_count"]),
                str(r["ref_CA_len"]),
                str(r["query_CA_len"]),
                r["pass_pLDDT"],
                r["pass_RMSD"],
                r["final_status"],
                r["fail_reason"],
                r["selected_pdb"],
            ]
            out.write("\t".join(row) + "\n")


write_table(OUT_ALL, results_sorted)
write_table(OUT_PASS, pass_results)


with open(OUT_PASS_FA, "w") as out:
    for r in pass_results:
        out.write(f">{r['sample_id']}\n")
        seq = r["seq"]

        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")


for r in pass_results:
    selected_pdb = r["selected_pdb"]

    if selected_pdb != "NA" and os.path.exists(selected_pdb):
        out_pdb = os.path.join(OUT_PASS_PDB_DIR, f"{r['sample_id']}.pdb")
        shutil.copy(selected_pdb, out_pdb)


print("完成")
print("输出目录:", OUT_DIR)
print("整合总表:", OUT_ALL)
print("过滤表:", OUT_PASS)
print("过滤后 fasta:", OUT_PASS_FA)
print("过滤后 PDB 目录:", OUT_PASS_PDB_DIR)
print("候选总数:", len(results))
print("通过筛选数量:", len(pass_results))
print(f"筛选条件: mean_pLDDT >= {PLDDT_CUTOFF}, RMSD_CA <= {RMSD_CUTOFF}")