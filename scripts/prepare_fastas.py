import sys
import os
from Bio import SeqIO

sample = snakemake.wildcards.sample
ext = snakemake.config["samples"][sample]
input_file = snakemake.input[0]
output_file = snakemake.output[0]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

if ext == "fna":
    records = list(SeqIO.parse(input_file, "fasta"))
    SeqIO.write(records, output_file, "fasta")
    print(f"[INFO] Copied {len(records)} sequences from .fna")
elif ext == "fastq":
    records = list(SeqIO.parse(input_file, "fastq"))
    SeqIO.write(records, output_file, "fasta")
    print(f"[INFO] Converted {len(records)} sequences from .fastq → .fna")
else:
    print(f"[ERROR] Unsupported format: {ext}")
    sys.exit(1)