from Bio import SeqIO
from BCBio import GFF
import argparse
import os
from logger import Logger

log = Logger(name="translate_orfs").get_logger()

def translate(fna_path, gff_path, faa_path, coords_path):
    log.info(f"Загружаем геном: {fna_path}")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))
    os.makedirs(os.path.dirname(faa_path), exist_ok=True)

    translated = 0
    with open(faa_path, "w") as faa_out, open(coords_path, "w") as coords_out:
        for rec in GFF.parse(gff_path, base_dict=seq_dict):
            for feat in rec.features:
                if feat.type == "CDS":
                    prot = feat.qualifiers.get("translation", [""])[0]
                    feat_id = feat.qualifiers.get("ID", ["unnamed"])[0]
                    start = feat.location.start + 1
                    end = feat.location.end
                    strand = "+" if feat.strand == 1 else "-"
                    faa_out.write(f">{feat_id}\n{prot}\n")
                    coords_out.write(f"{feat_id}\t{rec.id}\t{start}\t{end}\t{strand}\n")
                    translated += 1

    log.info(f"Переведено {translated} CDS → {faa_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fna", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--faa", required=True)
    parser.add_argument("--coords", required=True)
    args = parser.parse_args()
    translate(args.fna, args.gff, args.faa, args.coords)

if __name__ == "__main__":
    main()
else:
    translate(
        snakemake.input.fna,
        snakemake.input.gff,
        snakemake.output.faa,
        snakemake.output.coords
    )