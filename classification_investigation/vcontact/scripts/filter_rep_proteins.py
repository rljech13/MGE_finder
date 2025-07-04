#!/usr/bin/env python3
"""
Фильтрует белки по списку sample-ID.

Сохраняет запись, если в её описании:
  • присутствует  |sample=<ID>
        ИЛИ
  • последний пробел-разделённый токен равен <ID>
"""
import sys, re
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

if len(sys.argv) != 4:
    sys.exit("usage: filter_rep_proteins.py <all_proteins.faa> <samples.txt> <out.faa>")

faa_all, txt_samples, faa_out = sys.argv[1:]
samples = {s.strip() for s in open(txt_samples)}
kept = 0

with open(faa_out, "w") as fo:
    writer = FastaWriter(fo, wrap=None)

    for rec in SeqIO.parse(faa_all, "fasta"):
        desc = rec.description

        # 1) формат |sample=<ID>
        m = re.search(r"\bsample=([^\s|;,]+)", desc)
        sample_id = m.group(1) if m else None

        # 2) иначе берём последний токен
        if sample_id is None:
            sample_id = desc.strip().split()[-1].rstrip("|;,")
        if sample_id in samples:
            writer.write_record(rec)
            kept += 1

print(f"kept {kept} proteins → {faa_out}")