import os
import yaml
import glob

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "/home/lam34/MGE_finder/config.yaml")
configfile: CONFIG_PATH

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

GENOMES_DIR = config["paths"]["genomes_dir"]
RESULTS_DIR = config["paths"]["results_dir"]
PFAM_LIST = config["pfam_profiles"]
PFAM_ARGS = " ".join(PFAM_LIST)
ENV = config["execution"]["conda_env"]

COMBINED_HMM = os.path.join(RESULTS_DIR, "combined", "pfam_combined.hmm")

def get_samples():
    fasta_files = glob.glob(os.path.join(GENOMES_DIR, "*.fna"))
    return [os.path.splitext(os.path.basename(p))[0] for p in fasta_files]

SAMPLES = get_samples()

rule all:
    input:
        COMBINED_HMM,
        expand("{results}/{sample}/integrase_hits.txt", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/integrase_hits_summary.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/integrase_orfs.tsv", results=RESULTS_DIR, sample=SAMPLES)

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
        ffn=os.path.join(RESULTS_DIR, "{sample}", "orfs.ffn"),
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "predict_orfs.log")
    conda:
        ENV
    shell:
        "python scripts/predict_orfs.py "
        "--fna {input.fna} --gff {output.gff} --ffn {output.ffn} --faa {output.faa} "
        "> {log} 2>&1"

rule build_combined_hmm:
    output:
        os.path.join(RESULTS_DIR, "combined", "pfam_combined.hmm")
    log:
        "logs/build_combined_hmm.log"
    conda:
        ENV
    run:
        import os

        hmm_output = output[0]
        hmm_dir = os.path.dirname(hmm_output)
        os.makedirs(hmm_dir, exist_ok=True)

        with open(hmm_output, "w") as out:
            for pfam in config["pfam_profiles"]:
                with open(pfam) as f:
                    out.write(f.read())

        shell("hmmpress {output} > {log} 2>&1")

rule hmm_search:
    input:
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa"),
        hmm=COMBINED_HMM
    output:
        hits=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits.txt"),
        stats=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        orfs=os.path.join(RESULTS_DIR, "{sample}", "integrase_orfs.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "hmm_search.log")
    conda:
        ENV
    shell:
        """
        python scripts/hmm_search.py \
            --faa {input.faa} \
            --out {output.hits} \
            --summary {output.stats} \
            --orfs {output.orfs} \
            --pfam {config[pfam_profiles]} \
            --combined {input.hmm} \
            > {log} 2>&1
        """