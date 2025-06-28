#!/usr/bin/env python3
"""
Фильтрует белки по списку sample-ID.
Оставляет только те записи, у которых в заголовке есть |sample=<ID>.
Использование:
    filter_rep_proteins.py all_proteins.faa samples.txt out_proteins.faa
"""
import sys, re, Bio.SeqIO as b

if len(sys.argv) != 4:
    sys.exit("usage: filter_rep_proteins.py <all_proteins.faa> <samples.txt> <out.faa>")

faa_all, txt_samples, faa_out = sys.argv[1:]
samples = {s.strip() for s in open(txt_samples)}

kept = 0
with open(faa_out, "w") as fo:
    writer = b.FastaIO.FastaWriter(fo, wrap=None)
    for rec in b.parse(faa_all, "fasta"):
        m = re.search(r"\|sample=([^\|\s]+)", rec.id)
        if m and m.group(1) in samples:
            writer.write_record(rec)
            kept += 1

print(f"kept {kept} proteins → {faa_out}")