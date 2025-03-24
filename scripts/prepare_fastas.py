import os
import shutil
import glob
from Bio import SeqIO
import argparse
import yaml
from logger import Logger

log = Logger(name="prepare_fasta", draw_progress=True).get_logger()


def process_input_sources(sources, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    processed = []
    task = None

    total_files = 0
    for src in sources:
        total_files += len(glob.glob(os.path.join(src, "GCA_*", "*_genomic.fna")))
        total_files += len(glob.glob(os.path.join(src, "barcode*.fastq_final.fasta")))

    if total_files > 0:
        task = Logger().progress_task("Copying genomes", total=total_files)

    for src in sources:
        # === NCBI-style .fna ===
        fna_files = glob.glob(os.path.join(src, "GCA_*", "*_genomic.fna"))
        for fna in fna_files:
            # Получаем родительскую папку GCA_*
            sample = os.path.basename(os.path.dirname(fna))  # GCA_000008125.1
            dst = os.path.join(out_dir, f"{sample}.fna")
            try:
                shutil.copy(fna, dst)
                processed.append(sample)
                log.info(f"[NCBI] {sample} → {dst}")
            except Exception as e:
                log.error(f"Ошибка при копировании {fna}: {e}")
            Logger().advance_progress(task)

        # === Hybracter-style .fastq_final.fasta
        fastq_files = glob.glob(os.path.join(src, "barcode*.fastq_final.fasta"))
        for fasta in fastq_files:
            sample = os.path.basename(fasta).split(".")[0]  # barcode01.fastq_final.fasta → barcode01
            dst = os.path.join(out_dir, f"{sample}.fna")
            try:
                with open(fasta) as fin, open(dst, "w") as fout:
                    count = SeqIO.write(SeqIO.parse(fin, "fasta"), fout, "fasta")
                processed.append(sample)
                log.info(f"[HYBRACTER] {sample} ({count} seqs) → {dst}")
            except Exception as e:
                log.error(f"Ошибка при конвертации {fasta}: {e}")
            Logger().advance_progress(task)

    Logger().finish_progress(task)
    log.info(f"[✓] Всего обработано: {len(processed)} образцов.")
    return processed


def main(config_path, output_done=None):
    log.info(f"Загружаем конфиг: {config_path}")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    sources = config["input_sources"]
    out_dir = config["genomes_dir"]

    processed = process_input_sources(sources, out_dir)

    if output_done:
        with open(output_done, "w") as f:
            f.write("done\n")
        log.info(f"Файл завершения записан: {output_done}")


# === CLI + Snakemake compatibility ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare genome FASTA files")
    parser.add_argument("--config", required=False, default="config.yaml",
                        help="Path to config.yaml")
    parser.add_argument("--done", required=False,
                        help="(Optional) Path to .complete output file (used in Snakemake)")
    args = parser.parse_args()
    main(args.config, args.done)
else:
    # Snakemake hook
    main("config.yaml", snakemake.output[0])