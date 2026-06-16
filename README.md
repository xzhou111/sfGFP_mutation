```bash
mkdir 1.pdb 2.model 3.candidates 4.filter
```
```bash
  1.pdb
     ├── Exclusion_List.txt      #排除序列
     ├── fixed_positions.txt     #固定区
     ├── mutable_positions.txt   #可变区
     ├── sfGFP.fa                #蛋白序列
     └── sfGFP.pdb               #结构文件

    固定区原理:
      发色团SYG及附近60–72，保证能形成荧光发色团
      发色团口袋关键残基,保持亮度和局部环境
      β-barrel 关键核心位点,保持桶状结构稳定
```
  2.model:按照下面网址，提前下载蛋白模型ESM，安装ProteinMPNN
  ```bash
    https://huggingface.co/facebook/esm2_t30_150M_UR50D
    https://github.com/dauparas/ProteinMPNN
  ```
  3.candidates：
    利用ProteinMPNN根据一个已经给定的蛋白三维骨架，重新设计适合这个骨架的氨基酸序列<br>
```python
    python script/bin/01_make_fixed_jsonl.py    #得到可识别的fixed_positions.jsonl文件
    python 2.model/ProteinMPNN/protein_mpnn_run.py \
      --pdb_path 1.pdb/sfGFP.pdb \
      --pdb_path_chains A \
      --fixed_positions_jsonl 3.candidates/01_proteinmpnn/fixed_positions.jsonl \
      --out_folder 3.candidates/01_proteinmpnn \
      --num_seq_per_target 3 \
      --sampling_temp "0.1 0.2 0.3 0.5 0.8" \
      --model_name v_48_020 \
      --path_to_model_weights 2.model/ProteinMPNN/vanilla_model_weights \
      --batch_size 2 \
      --seed 2026 \
      --save_score 1
      num_seq_per_target：表示每个目标结构生成多少条候选序列
      sampling_temp：采样温度，控制序列变化程度，温度越低，生成的序列越保守；温度越高，生成的序列变化越大
      model_name：使用 ProteinMPNN 的模型版本
      path_to_model_weights：指定 ProteinMPNN 模型权重的位置
      batch_size ：每次同时生成几条序列，报显存不足或内存不足，就改小
      seed：设置随机种子，让结果可重复
      save_score：对每条序列的打分

    #序列名称重命名
      awk 'BEGIN{n=0; keep=0} /^>/{if($0 ~ /^>sfGFP/){keep=0; next}else{n++; print ">sample"n; keep=1; next}} keep{print}'       3.candidates/01_proteinmpnn/seqs/sfGFP.fa > 3.candidates/01_proteinmpnn/sfGFP.fa
    #对ProteinMPNN生成的 sfGFP 候选序列进行 ESM-2 打分，然后按分数从高到低排序，筛出 Top 候选序列
      python script/bin/02_esm2_score_top.py
```

  4.filter：再通过pLDDT和RMSD进一步筛选
    预测蛋白结构(也可以通过alphafold3预测)
```bash
    colabfold_batch 3.candidates/02_esm2_t33_score/Top_candidates.fa  4.filter/01_colabfol
    python  script/bin/03_integrate_plddt_rmsd.py        #colabfold_batch预测结果
    python  script/bin/03_integrate_plddt_rmsd_af3.py    #alphafold3预测结果
```
