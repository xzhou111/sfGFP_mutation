#!/usr/bin/env python3
import os
import re
import json
from collections import OrderedDict

PDB_FILE = "1.pdb/sfGFP.pdb"
WT_FASTA = "1.pdb/sfGFP.fa"
FIXED_TXT = "1.pdb/fixed_positions.txt"
OUT_JSONL = "3.candidates/01_proteinmpnn/fixed_positions.jsonl"
CHAIN = "A"

def read_fasta(path):
    records = OrderedDict()
    name = None
    seqs = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seqs).upper()
            name = line[1:].strip()
            seqs = []
        else:
            seqs.append(line)
    if name is not None:
        records[name] = "".join(seqs).upper()
    return records

def read_positions(path):
    positions = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        line = line.replace(",", " ")
        for x in line.split():
            if re.fullmatch(r"\d+", x):
                positions.append(int(x))
            elif re.fullmatch(r"\d+-\d+", x):
                a, b = map(int, x.split("-"))
                positions.extend(range(min(a, b), max(a, b) + 1))
    return sorted(set(positions))

wt = list(read_fasta(WT_FASTA).values())[0]
fixed_pos = read_positions(FIXED_TXT)
pdb_name = os.path.splitext(os.path.basename(PDB_FILE))[0]

bad = [p for p in fixed_pos if p < 1 or p > len(wt)]
if bad:
    raise ValueError(f"These fixed positions are out of range 1-{len(wt)}: {bad[:20]}")

data = {
    pdb_name: {
        CHAIN: fixed_pos
    }
}

os.makedirs(os.path.dirname(OUT_JSONL), exist_ok=True)

with open(OUT_JSONL, "w") as f:
    f.write(json.dumps(data) + "\n")

print("PDB name:", pdb_name)
print("Chain:", CHAIN)
print("WT length:", len(wt))
print("Fixed positions:", len(fixed_pos))
print("Designed positions:", len(wt) - len(fixed_pos))
print("Output:", OUT_JSONL)
