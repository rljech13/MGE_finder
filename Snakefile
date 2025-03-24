import os
import yaml
import glob

# Load config
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "MGE_finder/config.yaml")
configfile: CONFIG_PATH

GENOMES_DIR = config["paths"]["genomes_dir"]
RESULTS_DIR = config["paths"]["results_dir"]
PFAM = config["pfam_profiles"]
ENV = config["execution"]["conda_env"]

def get_samples():
    fasta_files = glob.glob(os.path.join(GENOMES_DIR, "*.fna"))
    return [os.path.splitext(os.path.basename(p))[0] for p in fasta_files]

SAMPLES = get_samples()

rule all:
    input:
        expand("{results}/{sample}/integrase_hits.txt", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/combined.hmm", results=RESULTS_DIR, sample=SAMPLES)

rule prepare_fasta:
    output:
        touch(os.path.join(GENOMES_DIR, ".complete"))
    log:
        "logs/prepare_fasta.log"
    conda:
        ENV
    shell:
        "python scripts/prepare_fastas.py --config config.yaml > {log} 2>&1"

rule predict_orfs:
    input:
        fna=lambda wildcards: os.path.join(GENOMES_DIR, f"{wildcards.sample}.fna")
    output:
        gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff"),
        ffn=os.path.join(RESULTS_DIR, "{sample}", "orfs.ffn")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "predict_orfs.log")
    conda:
        ENV
    shell:
        "python scripts/predict_orfs.py "
        "--fna {input.fna} --gff {output.gff} --ffn {output.ffn} "
        "> {log} 2>&1"

rule translate_orfs:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna"),
        gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff")
    output:
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa"),
        coords=os.path.join(RESULTS_DIR, "{sample}", "orfs.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "translate_orfs.log")
    conda:
        ENV
    shell:
        "python scripts/translate_orfs.py "
        "--fna {input.fna} --gff {input.gff} "
        "--faa {output.faa} --coords {output.coords} "
        "> {log} 2>&1"

rule hmm_search:
    input:
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa")
    output:
        hits=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits.txt"),
        hmm=os.path.join(RESULTS_DIR, "{sample}", "combined.hmm")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "hmm_search.log")
    conda:
        ENV
    shell:
        "python scripts/hmm_search.py "
        "--faa {input.faa} "
        "--out {output.hits} "
        "--pfam {config[pfam_profiles]} "
        "--combined {output.hmm} "
        "> {log} 2>&1"