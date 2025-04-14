## MGE Finder: A Snakemake Pipeline for Mobile Genetic Element Discovery

MGE Finder is a modular and extensible Snakemake pipeline for detecting Mobile Genetic Elements (MGEs), such as integrases near tRNAs, and annotating their genomic context. It integrates ORF prediction, HMM-based domain search, tRNA detection, BLAST-based boundary detection, and GenBank annotation.

⸻

### Project Structure
```
MGE_finder/
├── config.yaml                  # Main configuration file
├── Snakefile                    # Snakemake pipeline
├── scripts/                     # Auxiliary scripts used by pipeline steps
│   ├── prepare_fastas.py
│   ├── predict_orfs.py
│   ├── hmm_search.py
│   ├── annotate_trna_proximity.py
│   ├── extract_trna_region.py
│   ├── annotate_mge_region.py
│   ├── extract_mge_regions.py
│   └── annotate_and_orient_mge.py
├── results/                     # All output results are stored here
└── run.sh                       # Example script to launch the pipeline
```


⸻

### Installation

1. Make sure you have conda or mamba installed.
2. Create the environment:

```
mamba env create -f envs/MGE_finder.yaml
```
3. Activate the environment:

```bash
conda activate MGE_finder
```


⸻

### Running the Pipeline

bash run.sh

Optional flags:

	•	--dry-run or -n: preview steps without running them
	•	--unlock: unlock the working directory after interruption
	•	--force: force re-execution of rules
	•	--rerun-incomplete: rerun only incomplete jobs

⸻

Pipeline Overview

Step in Snakefile	Description
0. prepare_fasta	Converts and copies input .fna genome files to a consistent format
1. predict_orfs	Predicts ORFs using Prodigal (via a Python wrapper)
2. build_combined_hmm	Merges Pfam HMM profiles (e.g., PF00589, PF22022) into a single HMM database and runs hmmpress
3. hmm_search	Scans protein sequences for integrase domains using hmmscan and filters the best hits
4. predict_trna	Detects tRNAs using ARAGORN
5. trna_proximity	Identifies tRNAs located close to integrase hits (default max distance = 500 bp)
6. extract_trna_region	Extracts short regions around selected tRNAs to be used as MGE boundaries
7. blast_mge	Locates regions between integrase and nearby tRNAs by BLASTing extracted tRNA sites
8. extract_mge_region	Extracts full MGE candidate regions from genomes
9. annotate_mge	Annotates MGE regions with integrases, tRNAs, and attachment (att) sites in GenBank format



⸻

Example config.yaml

paths:
  genomes_dir: "data/genomes"       # Path to input genome FASTA files (*.fna)
  results_dir: "results"            # Where results will be saved

execution:
  conda_env: "MGE_finder"           # Conda environment name

pfam_profiles:
  - "pfam/PF00589.27.hmm"
  - "pfam/PF22022.2.hmm"



⸻

Logs

Each rule has an associated log file saved under the logs/ directory (or inside results/{sample}/ for sample-specific steps).

Example log message format:

[2025-04-14 10:32:10 - INFO] [hmm_search] 24 integrases detected on AE017221.1



⸻

Requirements
	•	Python 3.8+
	•	Snakemake ≥ 7.x
	•	Biopython
	•	Prodigal or Pyrodigal
	•	ARAGORN
	•	HMMER 3.x
	•	BLAST+
	•	BCBio.GFF

⸻

Outputs

For each genome sample, the pipeline produces:
	•	orfs.gff, orfs.faa, orfs.ffn: ORF predictions
	•	integrase_hits.txt, integrase_hits_summary.tsv, integrase_orfs.tsv: HMM results
	•	trna.tsv: ARAGORN-predicted tRNAs
	•	integrase_trna.tsv: Nearby integrase-tRNA pairs
	•	mge_query.fa: Query fragments for BLAST
	•	mge_blast.tsv: BLAST hits linking tRNA and integrase
	•	mge_region.fa: Extracted MGE regions
	•	mge_annotated.gbk: Annotated MGE regions (GenBank)
	•	attachment_sites.tsv: Detected attL/attR coordinates

⸻

TODO / Future Features
	-	Add plots/stats for integrase and MGE counts
	-	Enable large-scale comparative analysis across genomes
	-	Optional MCAPP-compatible output
	-	Docker/Singularity container for full reproducibility

⸻

Additionals

/scripts dir contains two more scripts:

logger.py - wrapper around standart python logger(more convieninet way; shotout to zkh-dot)

summary.py - right now under construction - script that would make general summary about all things that can be used, or seen during the work of pipeline

