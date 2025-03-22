from Bio import SeqIO
from BCBio import GFF
import os

sample = snakemake.wildcards.sample
fna_path = snakemake.input.fna
gff_path = snakemake.input.gff
faa_path = snakemake.output.faa
coords_path = snakemake.output.coords

seq_dict = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))
with open(faa_path, "w") as faa_out, open(coords_path, "w") as coords_out:
    for rec in GFF.parse(gff_path, base_dict=seq_dict):
        for feat in rec.features:
            if feat.type == "CDS":
                prot = feat.qualifiers["translation"][0]
                id = feat.qualifiers["ID"][0]
                start = feat.location.start + 1
                end = feat.location.end
                strand = "+" if feat.strand == 1 else "-"
                faa_out.write(f">{id}\n{prot}\n")
                coords_out.write(f"{id}\t{rec.id}\t{start}\t{end}\t{strand}\n")