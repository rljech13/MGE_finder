#!/usr/bin/env python3
"""
annotate_mge_region.py
Выбирает лучший BLAST-хит тРНК в окне рядом с интегразой и
сохраняет все «сырые» хиты в один общий файл.
"""
import os
import subprocess
import argparse
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from logger import Logger

logger = Logger(name="annotate_mge_region",
                level=Logger.Level.INFO).get_logger()

WINDOW_SIZE = 300_000          # длина окна вокруг интегразы
SHIFT = 3                      # сколько нуклеотидов от края тРНК считаем «совпадением»

# порядок полей из -outfmt 6
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def load_integrases(integrases_file: str) -> pd.DataFrame:
    return pd.read_csv(integrases_file, sep="\t")

def get_integrase_record(df: pd.DataFrame, integrase_id: str):
    rec = df[df["orf_id"] == integrase_id]
    if rec.empty:
        logger.error(f"Integrase ID {integrase_id} not found in {integrases_file}")
        return None
    return rec.iloc[0]

def extract_subject_region(ffn: str, contig: str,
                           region_start: int, region_end: int):
    recs = SeqIO.to_dict(SeqIO.parse(ffn, "fasta"))
    rec = recs.get(contig)
    if not rec:                        # иногда id в FASTA длиннее
        rec = next((r for k, r in recs.items() if contig in k), None)
        if rec:
            logger.info(f"Using contig {k} for match {contig}")
    if not rec:
        logger.error(f"Contig {contig} not found")
        return None
    L = len(rec.seq)
    s = max(1, region_start)
    e = min(L, region_end)
    return rec.seq[s - 1:e]

def run_blast_on_region(query_rec: SeqRecord, subject_seq,
                        tmp_dir=".") -> pd.DataFrame:
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        qfa = os.path.join(tmp, "query.fa")
        sfa = os.path.join(tmp, "subject.fa")
        bout = os.path.join(tmp, "blast.tsv")

        SeqIO.write([query_rec], qfa, "fasta")
        SeqIO.write([SeqRecord(subject_seq,
                               id="subject", description="")], sfa, "fasta")

        subprocess.run(["makeblastdb", "-in", sfa, "-dbtype", "nucl"],
                       check=True, capture_output=True, text=True)
        logger.info("makeblastdb OK")

        subprocess.run([
            "blastn", "-query", qfa, "-db", sfa,
            "-outfmt", f"6 {' '.join(BLAST_COLS)}",
            "-word_size", "4", "-dust", "no",
            "-out", bout
        ], check=True, capture_output=True, text=True)
        logger.info("blastn OK")

        try:
            return pd.read_csv(bout, sep="\t",
                               header=None, names=BLAST_COLS)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=BLAST_COLS)

def main(ffn, integrases_file, query_fa, out_tsv, tmp_dir="."):
    ibi = load_integrases(integrases_file)
    records = []          # отфильтрованные «лучшие» хиты
    all_raw = []          # все сырые хиты

    for rec in SeqIO.parse(query_fa, "fasta"):
        try:
            iid, contig, rng, strand = rec.id.strip().split(":")
            start, end = map(int, rng.split("-"))
        except ValueError:
            logger.error(f"Bad header {rec.id}")
            continue

        integrase = get_integrase_record(ibi, iid)
        if integrase is None:
            continue

        # Определяем окно вокруг интегразы
        if strand == "+":
            wstart, wend = end, end + WINDOW_SIZE
        else:
            wstart, wend = max(1, start - WINDOW_SIZE), start

        subj = extract_subject_region(ffn, contig, wstart, wend)
        if not subj:
            logger.error(f"No region window for {rec.id}")
            continue

        df_blast = run_blast_on_region(
            SeqRecord(rec.seq, id=rec.id, description=""),
            subj, tmp_dir
        )
        if df_blast.empty:
            logger.info(f"No BLAST for {rec.id}")
            continue

        # --- сохраняем все хиты ---
        df_blast["integrase_id"] = iid
        df_blast["contig"] = contig
        df_blast["wstart"] = wstart
        all_raw.append(df_blast.copy())

        # --- фильтрация по SHIFT ---
        trna_len = len(rec.seq)
        if strand == "+":
            df_hit = df_blast[(df_blast.qstart <= SHIFT) |
                              (df_blast.qend >= trna_len - SHIFT + 1)]
        else:
            df_hit = df_blast[(df_blast.qend >= trna_len - SHIFT + 1) |
                              (df_blast.qstart <= SHIFT)]

        if df_hit.empty:
            logger.info(f"No hits within shift for {rec.id}")
            continue

        # --- выбираем лучший хит: min evalue, max bitscore ---
        best = df_hit.sort_values(
            by=["evalue", "bitscore"],
            ascending=[True, False]
        ).iloc[0]

        # абсолютные координаты в исходном континге
        full_s = wstart + int(best.sstart) - 1
        full_e = wstart + int(best.send) - 1

        records.append({
            "integrase_id": iid,
            "contig": contig,
            "hit_start": full_s,
            "hit_end": full_e,
            "pident": best.pident,
            "length": int(best.length),
            "evalue": best.evalue,
            "bitscore": best.bitscore,
            "qstart": int(best.qstart),
            "qend": int(best.qend),
        })

    # ---------- запись результатов ----------
    # 1. отфильтрованные «лучшие» хиты
    if records:
        pd.DataFrame(records).to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"{len(records)} filtered hits saved to {out_tsv}")
    else:
        open(out_tsv, "w").close()
        logger.info("No valid BLAST hits; wrote empty file.")

    # 2. полный «сырый» вывод BLAST
    raw_path = out_tsv.replace(".tsv", "_raw.tsv")
    if all_raw:
        pd.concat(all_raw).to_csv(raw_path, sep="\t", index=False)
        logger.info(f"{sum(len(df) for df in all_raw)} raw hits saved to {raw_path}")
    else:
        logger.info("No raw BLAST hits collected.")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ffn", required=True)
    ap.add_argument("--integrases", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--tmp_dir", default=".")
    args = ap.parse_args()
    main(args.ffn, args.integrases, args.query, args.out_tsv, args.tmp_dir)