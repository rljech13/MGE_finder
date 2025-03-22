import pyrodigal
import os
from Bio import SeqIO

sample = snakemake.wildcards.sample
fna_path = snakemake.input[0]
gff_path = snakemake.output.gff
ffn_path = snakemake.output.ffn

gf = pyrodigal.GeneFinder(meta=True)

with open(gff_path, "w") as gff_out, open(ffn_path, "w") as ffn_out:
    for record in SeqIO.parse(fna_path, "fasta"):
        genes = gf.find_genes(str(record.seq))
        for i, gene in enumerate(genes):
            gene_id = f"{record.id}_orf{i+1}"
            gff_out.write(gene.as_gff(seq_id=record.id, source="pyrodigal", gene_id=gene_id))
            ffn_out.write(f">{gene_id}\n{gene.nucleotide_sequence()}\n")