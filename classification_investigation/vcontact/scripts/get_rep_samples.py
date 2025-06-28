#!/usr/bin/env python3
"""
Собирает список sample-ID из FASTA репрезентативных MGE.
Использование:
    get_rep_samples.py rep_mge_nt.fa out_samples.txt
"""
import sys, re, Bio.SeqIO as b

if len(sys.argv) != 3:
    sys.exit("usage: get_rep_samples.py <rep_mge_nt.fa> <out_samples.txt>")

fa_in, txt_out = sys.argv[1:]
samples = set()

for rec in b.parse(fa_in, "fasta"):
    m = re.search(r"\|sample=([^\|\s]+)", rec.id)
    if m:
        samples.add(m.group(1))

with open(txt_out, "w") as f:
    for s in sorted(samples):
        f.write(f"{s}\n")

print(f"{len(samples)} unique samples written to {txt_out}")