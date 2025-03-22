import os
import glob

GENOMES = glob.glob("data/genomes/*")
GENOME_IDS = [os.path.basename(p).rsplit('.', 1)[0] for p in GENOMES]

rule all:
    input:
        expand("results/{sample}/integrase_hits.txt", sample=GENOME_IDS)

rule convert_fastq_to_fasta:
    input:
        "data/genomes/{sample}.fastq"
    output:
        "results/{sample}/converted.fna"
    shell:
        """
        seqtk seq -A {input} > {output}
        """

rule copy_fna:
    input:
        "data/genomes/{sample}.fna"
    output:
        "results/{sample}/converted.fna"
    shell:
        """
        cp {input} {output}
        """

rule predict_orfs:
    input:
        fasta="results/{sample}/converted.fna"
    output:
        proteins="results/{sample}/orfs.faa",
        coords="results/{sample}/orfs.tsv"
    conda:
        "envs/pyrodigal.yaml"
    script:
        "scripts/predict_orfs.py"

rule hmmsearch_integrases:
    input:
        proteins="results/{sample}/orfs.faa",
        hmm1="pfam/PF00589.27.hmm",
        hmm2="pfam/PF22022.2.hmm"
    output:
        hits="results/{sample}/integrase_hits.txt"
    conda:
        "envs/pyrodigal.yaml"
    shell:
        """
        cat {input.hmm1} {input.hmm2} > results/{wildcards.sample}/combined.hmm
        hmmpress results/{wildcards.sample}/combined.hmm
        hmmscan --tblout {output.hits} results/{wildcards.sample}/combined.hmm {input.proteins}
        """