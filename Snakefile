configfile: "/home/lam34/MGE_finder/config.yaml"

import yaml
import os
import glob

with open("/home/lam34/MGE_finder/config.yaml") as f:
    config = yaml.safe_load(f)

ENV = config["env"]
PFAM = config["pfam_profiles"]
GENOMES_DIR = config["genomes_dir"]
RESULTS_DIR = config["results_dir"]

# Получаем список образцов из data/genomes/*.fna
def get_samples():
    fasta_files = glob.glob(os.path.join(GENOMES_DIR, "*.fna"))
    return [os.path.splitext(os.path.basename(p))[0] for p in fasta_files]

SAMPLES = get_samples()

rule all:
    input:
        expand("{results}/{sample}/integrase_hits.txt", results=RESULTS_DIR, sample=SAMPLES)

rule prepare_fasta:
    output:
        touch("data/genomes/.complete")
    conda:
        ENV
    script:
        "scripts/prepare_fastas.py"

rule predict_orfs:
    input:
        f"{GENOMES_DIR}" + "/{sample}.fna"
    output:
        gff=f"{RESULTS_DIR}" + "/{sample}/orfs.gff",
        ffn=f"{RESULTS_DIR}" + "/{sample}/orfs.ffn"
    conda:
        ENV
    script:
        "scripts/predict_orfs.py"

rule translate_orfs:
    input:
        fna=f"{GENOMES_DIR}" + "/{sample}.fna",
        gff=f"{RESULTS_DIR}" + "/{sample}/orfs.gff"
    output:
        faa=f"{RESULTS_DIR}" + "/{sample}/orfs.faa",
        coords=f"{RESULTS_DIR}" + "/{sample}/orfs.tsv"
    conda:
        ENV
    script:
        "scripts/translate_orfs.py"

rule hmm_search:
    input:
        faa=f"{RESULTS_DIR}" + "/{sample}/orfs.faa"
    output:
        hits=f"{RESULTS_DIR}" + "/{sample}/integrase_hits.txt"
    conda:
        ENV
    script:
        "scripts/hmm_search.py"