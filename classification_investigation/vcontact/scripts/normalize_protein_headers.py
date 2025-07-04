#!/usr/bin/env python3
"""
Из <rep_proteins_raw.faa> делает:
  • rep_proteins_normalized.faa   – header == protein_id
  • gene2genome.csv               – protein_id,contig_id,genome_id,keywords
                                    (заголовок ОБЯЗАТЕЛЕН)

* protein_id = первый «слово-токен» описания, очищенный.
* sample_id  =  |sample=<ID>  → <ID>
                 иначе         → последний токен описания
  sample_id пишется и в contig_id, и в genome_id.
"""
import argparse, csv, re, sys
from pathlib import Path
from Bio import SeqIO

ALLOWED = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789|_-.")

def safe(text: str) -> str:
    return "".join(c for c in text if c in ALLOWED)

SID_RE = re.compile(r"\bsample=([^\s|;,]+)")

def pull_ids(desc: str):
    pid = safe(desc.split()[0])               # protein_id
    m   = SID_RE.search(desc)
    sid = m.group(1) if m else desc.split()[-1].rstrip("|;,")
    return pid, safe(sid)

def main(inp: Path, out_faa: Path, out_csv: Path):
    dup, seen = 0, set()
    with out_faa.open("w") as fout, out_csv.open("w", newline="") as cout:
        w = csv.writer(cout, lineterminator="\n")
        w.writerow(["protein_id","contig_id","genome_id","keywords"])   # <--- ключевая строка

        for rec in SeqIO.parse(str(inp), "fasta"):
            pid, sid = pull_ids(rec.description)

            if pid in seen:
                dup += 1
                pid = f"{pid}__dup{dup}"
            seen.add(pid)

            # пишем FASTA
            rec.id = pid;  rec.name = rec.description = ""
            SeqIO.write(rec, fout, "fasta")

            # пишем CSV (keywords пустые)
            w.writerow([pid, sid, sid, ""])

    sys.stderr.write(f"{len(seen)} proteins written (duplicates {dup})\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("raw_faa",  type=Path)
    p.add_argument("out_faa",  type=Path)
    p.add_argument("gene2genome_csv", type=Path)
    args = p.parse_args()

    args.out_faa.parent.mkdir(parents=True, exist_ok=True)
    args.gene2genome_csv.parent.mkdir(parents=True, exist_ok=True)
    main(args.raw_faa, args.out_faa, args.gene2genome_csv)