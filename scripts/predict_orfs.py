import argparse
import subprocess
import os
from logger import Logger

log = Logger(name="predict_orfs").get_logger()

def predict_with_prodigal(fna_path, gff_path, ffn_path):
    os.makedirs(os.path.dirname(gff_path), exist_ok=True)
    os.makedirs(os.path.dirname(ffn_path), exist_ok=True)

    cmd = [
        "prodigal",
        "-i", fna_path,
        "-o", gff_path,
        "-a", os.devnull,  # Не сохраняем .faa здесь
        "-d", ffn_path,
        "-f", "gff",
        "-p", "meta"  # режим для метагеномов
    ]

    log.info(f"Выполняем prodigal:\n{' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        log.info(f"ORF-аннотация завершена для: {fna_path}")
        log.info(f"GFF: {gff_path}")
        log.info(f"FFN: {ffn_path}")
    except subprocess.CalledProcessError as e:
        log.error(f"Ошибка при запуске prodigal: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Predict ORFs using Prodigal CLI")
    parser.add_argument("--fna", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--ffn", required=True)
    args = parser.parse_args()

    predict_with_prodigal(args.fna, args.gff, args.ffn)

if __name__ == "__main__":
    main()
else:
    predict_with_prodigal(
        snakemake.input.fna,
        snakemake.output.gff,
        snakemake.output.ffn
    )