import os
import argparse
from logger import Logger

log = Logger(name="hmm_search").get_logger()

def run_hmmscan(faa_path, out_path, pfam_profiles, combined_hmm_path):
    log.info("Объединяем профили Pfam...")
    os.makedirs(os.path.dirname(combined_hmm_path), exist_ok=True)

    with open(combined_hmm_path, "w") as out_hmm:
        for pfam in pfam_profiles:
            with open(pfam) as f:
                out_hmm.write(f.read())

    log.info(f"Pfam объединён: {combined_hmm_path}")
    log.info("Индексируем hmm-базу через hmmpress...")
    os.system(f"hmmpress {combined_hmm_path}")

    log.info("Запускаем hmmscan...")
    os.system(f"hmmscan --tblout {out_path} {combined_hmm_path} {faa_path}")
    log.info(f"Результаты сохранены: {out_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--faa", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--pfam", required=True, nargs="+")
    parser.add_argument("--combined", required=True)
    args = parser.parse_args()

    run_hmmscan(args.faa, args.out, args.pfam, args.combined)

if __name__ == "__main__":
    main()
else:
    run_hmmscan(
        snakemake.input.faa,
        snakemake.output.hits,
        snakemake.config["pfam_profiles"],
        snakemake.output.hmm
    )