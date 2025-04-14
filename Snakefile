import os
import yaml
import glob

CONFIG_PATH = "/home/lam34/MGE_finder/config.yaml"
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
    if not fasta_files:
        raise FileNotFoundError(f"No .fna files found in {GENOMES_DIR}")
    return [os.path.splitext(os.path.basename(p))[0] for p in fasta_files]

SAMPLES = get_samples()

rule all:
    input:
        COMBINED_HMM,
        expand("{results}/{sample}/integrase_hits.txt", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/integrase_hits_summary.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/integrase_orfs.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/trna.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/integrase_trna.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/mge_query.fa", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/mge_blast.tsv", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/mge_region.fa", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/mge_annotated.gbk", results=RESULTS_DIR, sample=SAMPLES),
        expand("{results}/{sample}/attachment_sites.tsv", results=RESULTS_DIR, sample=SAMPLES)


rule prepare_fasta:
    output:
        os.path.join(GENOMES_DIR, ".complete")
    log:
        "logs/prepare_fasta.log"
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/prepare_fastas.py --config config.yaml > {log} 2>&1
        touch {output}
        """

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
        f"{ENV}"
    shell:
        """
        python scripts/predict_orfs.py --fna {input.fna} --gff {output.gff} --ffn {output.ffn} --faa {output.faa} > {log} 2>&1
        """

rule build_combined_hmm:
    output:
        os.path.join(RESULTS_DIR, "combined", "pfam_combined.hmm")
    log:
        "logs/build_combined_hmm.log"
    conda:
        f"{ENV}"
    run:
        hmm_output = output[0]
        os.makedirs(os.path.dirname(hmm_output), exist_ok=True)
        with open(hmm_output, "w") as out:
            for pfam in config["pfam_profiles"]:
                with open(pfam) as f:
                    out.write(f.read())
        shell(f"hmmpress {hmm_output} > {log}")

rule hmm_search:
    input:
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa"),
        gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff"),
        hmm=COMBINED_HMM
    output:
        hits=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits.txt"),
        stats=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        orfs=os.path.join(RESULTS_DIR, "{sample}", "integrase_orfs.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "hmm_search.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/hmm_search.py --faa {input.faa} --gff {input.gff} --out {output.hits} --summary {output.stats} --orfs {output.orfs} --combined {input.hmm} > {log} 2>&1
        """

rule predict_trna:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna")
    output:
        trna=os.path.join(RESULTS_DIR, "{sample}", "trna.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "predict_trna.log")
    conda:
        f"{ENV}"
    shell:
        """
        aragorn -w -t -o {output.trna} {input.fna} > {log} 2>&1
        """

rule trna_proximity:
    input:
        integrases=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        trna=os.path.join(RESULTS_DIR, "{sample}", "trna.tsv")
    output:
        proximity=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "trna_proximity.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/annotate_trna_proximity.py --integrases {input.integrases} --trna {input.trna} --output {output.proximity} --max_distance 500 > {log} 2>&1
        """

rule extract_trna_region:
    input:
        fasta=os.path.join(GENOMES_DIR, "{sample}.fna"),
        trnas=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv")
    output:
        out_fa=os.path.join(RESULTS_DIR, "{sample}", "mge_query.fa")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "extract_trna_region.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/extract_trna_region.py --ffn {input.fasta} --trnas {input.trnas} --out_fa {output.out_fa} > {log} 2>&1
        """

rule blast_mge:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna"),
        integrases=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        query=os.path.join(RESULTS_DIR, "{sample}", "mge_query.fa")
    output:
        blast_tsv=os.path.join(RESULTS_DIR, "{sample}", "mge_blast.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "blast_mge.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/annotate_mge_region.py --ffn {input.fna} --integrases {input.integrases} --query {input.query} --out_tsv {output.blast_tsv} --tmp_dir . > {log} 2>&1
        """

rule extract_mge_region:
    input:
        fna = os.path.join(GENOMES_DIR, "{sample}.fna"),
        blast = os.path.join(RESULTS_DIR, "{sample}", "mge_blast.tsv"),
        trna = os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv")
    output:
        mge_fa = os.path.join(RESULTS_DIR, "{sample}", "mge_region.fa")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "extract_mge_region.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/extract_mge_regions.py --fna {input.fna} --blast {input.blast} --trna {input.trna} --out_fa {output.mge_fa} > {log} 2>&1
        """

rule annotate_mge:
    input:
        fasta=os.path.join(RESULTS_DIR, "{sample}", "mge_region.fa"),
        trna=os.path.join(RESULTS_DIR, "{sample}", "mge_query.fa"),
        orf=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        blast=os.path.join(RESULTS_DIR, "{sample}", "mge_blast.tsv")
    output:
        gbk=os.path.join(RESULTS_DIR, "{sample}", "mge_annotated.gbk"),
        att=os.path.join(RESULTS_DIR, "{sample}", "attachment_sites.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "annotate_mge.log")
    conda:
        f"{ENV}"
    shell:
        """
        python scripts/annotate_and_orient_mge.py \
            --fasta {input.fasta} \
            --trna_fa {input.trna} \
            --integrase {input.orf} \
            --blast {input.blast} \
            --out_gbk {output.gbk} \
            --out_att {output.att} > {log} 2>&1
        """