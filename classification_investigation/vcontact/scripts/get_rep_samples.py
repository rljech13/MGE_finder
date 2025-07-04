#!/usr/bin/env python3
"""
Извлекает sample-ID из FASTA репрезентативных MGE.
Берёт именно rec.description, потому что rec.id содержит только locus-tag.

usage:
    get_rep_samples.py rep_mge_nt.fa out_samples.txt
"""
import sys, re
from Bio import SeqIO

if len(sys.argv) != 3:
    sys.exit("usage: get_rep_samples.py <rep_mge_nt.fa> <out_samples.txt>")

fa_in, txt_out = sys.argv[1:]
samples = set()

for rec in SeqIO.parse(fa_in, "fasta"):
    m = re.search(r"\|sample=([^\|\s]+)", rec.description)
    if m:
        samples.add(m.group(1))

with open(txt_out, "w") as f:
    for s in sorted(samples):
        f.write(f"{s}\n")

print(f"{len(samples)} unique samples written to {txt_out}")