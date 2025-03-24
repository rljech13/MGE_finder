import os
import shutil
import glob
from Bio import SeqIO

sources = snakemake.config["input_sources"]
out_dir = snakemake.config["genomes_dir"]
os.makedirs(out_dir, exist_ok=True)

processed = []

for src in sources:
    # === NCBI-style .fna ===
    fna_files = glob.glob(os.path.join(src, "GCA_*", "*_genomic.fna"))
    for fna in fna_files:
        sample = os.path.basename(fna).split("_")[0]
        dst = os.path.join(out_dir, f"{sample}.fna")
        shutil.copy(fna, dst)
        processed.append(sample)
        print(f"[NCBI] {sample} → {dst}")

    # === Hybracter .fastq_final.fasta ===
    fastq_files = glob.glob(os.path.join(src, "barcode*.fastq_final.fasta"))
    for fasta in fastq_files:
        sample = os.path.basename(fasta).split(".")[0]
        dst = os.path.join(out_dir, f"{sample}.fna")
        with open(fasta) as fin, open(dst, "w") as fout:
            count = SeqIO.write(SeqIO.parse(fin, "fasta"), fout, "fasta")
            processed.append(sample)
            print(f"[HYBRACTER] {sample} ({count} seqs) → {dst}")

# Отмечаем завершение
with open(snakemake.output[0], "w") as f:
    f.write("done\n")

print(f"[✔] Всего обработано: {len(processed)} образцов.")