#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd

def load_target_ids(cluster_stats: Path) -> dict[str, str]:
    """
    Загружает список MGE ID из cluster_stats.tsv
    Возвращает dict: {mge_id: sample}
    """
    df = pd.read_csv(cluster_stats, sep="\t", header=None)
    mapping = {}
    for row in df[0]:
        mge_id, sample = row.split("|")
        sample = sample.replace("sample=", "")
        mapping[mge_id] = sample
    return mapping

def extract_mge_from_gbff(genome_gbff: Path, target_ids: dict[str, str], out_dir: Path):
    """
    Читает gbff по геному и сохраняет только нужные записи (по LOCUS/ACCESSION).
    """
    records = list(SeqIO.parse(genome_gbff, "genbank"))
    saved = 0
    for rec in records:
        if rec.name in target_ids or rec.id in target_ids:
            mge_id = rec.name if rec.name in target_ids else rec.id
            sample = target_ids[mge_id]
            out_path = out_dir / f"{sample}_{mge_id}.gbff"
            SeqIO.write(rec, out_path, "genbank")
            saved += 1
    return saved

def main():
    parser = argparse.ArgumentParser(description="Извлекает отдельные MGE из gbff по списку cluster_stats.tsv")
    parser.add_argument("--cluster_stats", required=True, help="Файл cluster_stats.tsv")
    parser.add_argument("--bakta_dir", required=True, help="Директория с results_bakta (по семплам)")
    parser.add_argument("--out_dir", required=True, help="Куда класть выбранные gbff")
    args = parser.parse_args()

    cluster_stats = Path(args.cluster_stats)
    bakta_dir = Path(args.bakta_dir)
    out_dir = Path(args.out_dir)

    target_ids = load_target_ids(cluster_stats)

    for sample_dir in bakta_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        for gbff_file in sample_dir.glob("*.gbff"):
            n = extract_mge_from_gbff(gbff_file, target_ids, out_dir)
            if n > 0:
                print(f"[OK] {sample_dir.name}: извлечено {n} MGE из {gbff_file.name}")

if __name__ == "__main__":
    main()