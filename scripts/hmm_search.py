import os
import subprocess
import argparse
from logger import Logger
from collections import defaultdict
import pandas as pd

log = Logger(name="hmm_search").get_logger()

def run_hmmscan(faa_path, hmm_path, output_tbl):
    log.info(f"Запускаем hmmscan для {faa_path}...")
    subprocess.run([
        "hmmscan",
        "--tblout", output_tbl,
        hmm_path,
        faa_path
    ], check=True)
    log.info(f"Результаты сохранены: {output_tbl}")

def parse_tblout(tbl_path):
    """Вернёт dict: query_id -> list of profile accession"""
    mapping = defaultdict(list)
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                profile_acc = parts[1]
                query_id = parts[2]
                if profile_acc not in mapping[query_id]:
                    mapping[query_id].append(profile_acc)
    return mapping

def summarize(mapping, output_stats, output_orfs):
    counts = {"PF00589": 0, "PF22022": 0}
    detailed = []

    for orf_id, profiles in mapping.items():
        accs = [p for p in profiles if p.startswith("PF00589") or p.startswith("PF22022")]
        accs_set = set(a[:7] for a in accs)
        for pf in accs_set:
            counts[pf] += 1
        detailed.append((orf_id, ",".join(accs)))

    counts["total"] = counts["PF00589"] + counts["PF22022"]

    df1 = pd.DataFrame([counts])
    df1.to_csv(output_stats, sep="\t", index=False)

    df2 = pd.DataFrame(detailed, columns=["orf_id", "pfam_hits"])
    df2.to_csv(output_orfs, sep="\t", index=False)

    log.info(f"[✓] Статистика по интегразам сохранена в {output_stats}")
    log.info(f"[✓] Найденные ORF записаны в {output_orfs}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--faa", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--orfs", required=True)
    parser.add_argument("--pfam", required=True, nargs="+")  # игнорируется, но требуется Snakemake
    parser.add_argument("--combined", required=True)
    args = parser.parse_args()

    run_hmmscan(args.faa, args.combined, args.out)
    mapping = parse_tblout(args.out)
    summarize(mapping, args.summary, args.orfs)

if __name__ == "__main__":
    main()