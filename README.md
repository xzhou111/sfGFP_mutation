https://huggingface.co/facebook/esm2_t30_150M_UR50D
https://github.com/dauparas/ProteinMPNN


python script/bin/01_make_fixed_jsonl.py

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

awk 'BEGIN{n=0; keep=0} /^>/{if($0 ~ /^>sfGFP/){keep=0; next}else{n++; print ">sample"n; keep=1; next}} keep{print}' \
3.candidates/01_proteinmpnn/seqs/sfGFP.fa \
> 3.candidates/01_proteinmpnn/sfGFP.fa

python script/bin/02_esm2_score_top.py



/data3/Data_all/Software/miniconda3/envs/colabfold/bin/colabfold_batch 3.candidates/02_esm2_t33_score/Top_candidates.fa  4.filter/01_colabfol


python  script/bin/03_integrate_plddt_rmsd.py
