import yaml
import os

configfile: "config.yaml"

with open(configfile) as f:
    config = yaml.safe_load(f)

SAMPLES = list(config["samples"].keys())
GENOMES_DIR = config["genomes_dir"]
RESULTS_DIR = config["results_dir"]
ENV = config["env"]
PFAM = config["pfam_profiles"]

rule all:
    input:
        expand("{results}/{sample}/integrase_hits.txt", results=RESULTS_DIR, sample=SAMPLES)

rule prepare_fasta:
    input:
        lambda wc: f"{GENOMES_DIR}/{wc.sample}.{config['samples'][wc.sample]}"
    output:
        f"{RESULTS_DIR}" + "/{sample}/genome.fna"
    conda:
        ENV
    script:
        "scripts/prepare_fasta.py"

rule predict_orfs:
    input:
        f"{RESULTS_DIR}" + "/{sample}/genome.fna"
    output:
        gff=f"{RESULTS_DIR}" + "/{sample}/orfs.gff",
        ffn=f"{RESULTS_DIR}" + "/{sample}/orfs.ffn"
    conda:
        ENV
    script:
        "scripts/predict_orfs.py"

rule translate_orfs:
    input:
        fna=f"{RESULTS_DIR}" + "/{sample}/genome.fna",
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