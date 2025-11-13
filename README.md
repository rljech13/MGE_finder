# MGE Finder: Comprehensive Pipeline for Mobile Genetic Element Discovery and Analysis

MGE Finder is a modular and extensible Snakemake-based pipeline system for detecting, annotating, and analyzing Mobile Genetic Elements (MGEs) in prokaryotic genomes. The project integrates multiple analysis modules including ORF prediction, HMM-based domain search, tRNA detection, BLAST-based boundary detection, GenBank annotation, defense system detection, and comparative genomics.

---

## Table of Contents

- [Project Structure](#project-structure)
- [Installation](#installation)
- [Modules Overview](#modules-overview)
  - [1. Finder Pipeline](#1-finder-pipeline)
  - [2. Bakta Pipeline](#2-bakta-pipeline)
  - [3. PADLOC Pipeline](#3-padloc-pipeline)
  - [4. MGE Clustering](#4-mge-clustering)
  - [5. vContact Classification](#5-vcontact-classification)
  - [6. Protein Clustering](#6-protein-clustering)
  - [7. RefSeq Database](#7-refseq-database)
- [Additional Scripts](#additional-scripts)
- [Output Files](#output-files)
- [Configuration](#configuration)

---

##  Project Structure

```
MGE_finder/
├── finder_pipeline/          # Main MGE detection pipeline
│   ├── Snakefile
│   ├── finder_config.yaml
│   ├── run.sh
│   ├── scripts/              # Python scripts for pipeline steps
│   ├── data/genomes/         # Input genome FASTA files
│   └── results/              # Output results
│
├── bakta_pipeline/           # Bakta annotation pipeline
│   ├── Snakefile
│   ├── bakta_config.yaml
│   ├── run_bakta.sh
│   └── results_bakta/
│
├── padloc_pipeline/          # PADLOC defense system detection
│   ├── Snakefile
│   ├── padloc_config.yaml
│   ├── run_padloc.sh
│   └── results/
│
├── classification_investigation/
│   ├── clusterMGE/           # MGE clustering with MMseqs2
│   │   ├── Snakefile
│   │   ├── mge_cluster_config.yaml
│   │   └── run_mge_cluster.sh
│   └── vcontact/             # vContact2/vContact3 classification
│       ├── Snakefile
│       ├── vcontact_config.yaml
│       └── run_vcontact.sh
│
├── protein_clusterization/   # Protein clustering with MMseqs2
│   ├── Snakefile
│   ├── classification_config.yaml
│   └── run_clusterization.sh
│
├── refseq_prokaryote/        # RefSeq database download scripts
│   ├── download_genomes_only.sh
│   ├── bacteria/             # Downloaded bacterial genomes
│   └── archaea/             # Downloaded archaeal genomes
│
├── analysis/                 # Analysis scripts and results
│   ├── stats.py
│   ├── statistics/
│   └── plots/
│
├── additionals/              # Additional utility scripts
│   ├── run_padloc_parallel.sh
│   └── collect_reps.py
│
├── pfam/                     # Pfam HMM profiles
│   ├── PF00589.hmm
│   └── PF22022.hmm
│
├── envs/
│   └── MGE_finder.yaml       # Conda environment definition
│
└── README.md                 # This file
```

---

## 🔧 Installation

### Quick Start

For detailed setup instructions, see [SETUP.md](SETUP.md).

### Prerequisites

- Conda or Mamba
- Python 3.10+
- Snakemake ≥ 7.x
- At least 50GB free disk space (for results and databases)

### Setup

1. **Clone the repository:**

```bash
git clone <repository-url>
cd MGE_finder
```

2. **Create the main conda environment:**

```bash
mamba env create -f envs/MGE_finder.yaml
conda activate MGE_finder
```

3. **Install additional dependencies for specific modules:**

- **Bakta**: Requires separate `bakta` conda environment (see Bakta Pipeline section)
- **PADLOC**: Requires separate `padloc` conda environment (see PADLOC Pipeline section)
- **vContact**: Requires `vContact3` conda environment (see vContact section)

4. **Download Pfam HMM profiles:**

Place your Pfam HMM profiles (e.g., `PF00589.hmm`, `PF22022.hmm`) in the `pfam/` directory.

5. **Configure pipelines:**

Copy example configuration files and edit them with your paths:

```bash
# Finder pipeline
cd finder_pipeline
cp finder_config.yaml.example finder_config.yaml
# Edit finder_config.yaml with your paths

# Other pipelines similarly
cd ../bakta_pipeline
cp bakta_config.yaml.example bakta_config.yaml
# etc.
```

See [SETUP.md](SETUP.md) for detailed configuration instructions.

---

##  Modules Overview

### 1. Finder Pipeline

**Location:** `finder_pipeline/`

**Purpose:** Main pipeline for detecting MGEs integrating into tRNA genes. Screens genomes for integrases, identifies nearby tRNAs, extracts MGE regions, and annotates attachment sites.

**Workflow Steps:**

| Step | Description | Outputs |
|------|-------------|---------|
| `prepare_fasta` | Converts/copies input genomes to consistent `.fna` format | `genomes_dir/*.fna` |
| `predict_orfs` | Predicts ORFs using Prodigal | `orfs.gff`, `orfs.ffn`, `orfs.faa` |
| `build_combined_hmm` | Merges Pfam HMM profiles and runs `hmmpress` | `pfam_combined.hmm` |
| `hmm_search` | Scans proteins for integrase domains using `hmmscan` | `integrase_hits_summary.tsv`, `integrase_orfs.tsv` |
| `predict_trna` | Detects tRNAs using ARAGORN | `trna.tsv` |
| `trna_proximity` | Identifies tRNAs near integrases (≤500 bp) | `integrase_trna.tsv` |
| `extract_trna_region` | Extracts tRNA regions for BLAST queries | `mge_query.fa` |
| `blast_mge` | BLASTs tRNA sites to find MGE boundaries | `mge_blast.tsv` |
| `extract_mge_region` | Extracts full MGE candidate regions | `mge_region.fa` |
| `annotate_mge` | Annotates MGE regions in GenBank format | `mge_annotated.gbk`, `attachment_sites.tsv` |

**Configuration:** `finder_config.yaml`

Copy `finder_config.yaml.example` to `finder_config.yaml` and edit:

```yaml
paths:
  genomes_dir: "data/genomes"
  results_dir: "results"

execution:
  conda_env: "MGE_finder"

input_sources:
  - "/path/to/genomes/dir1"
  - "/path/to/genomes/dir2"

pfam_profiles:
  - "pfam/PF00589.hmm"
  - "pfam/PF22022.hmm"
```

**Run:**

```bash
cd finder_pipeline
bash run.sh
```

**Optional flags:**
- `--dry-run` or `-n`: Preview steps without running
- `--unlock`: Unlock working directory after interruption
- `--force` or `-f`: Force re-execution of rules
- `--rerun-incomplete`: Rerun only incomplete jobs
- `--cores=N` or `-j=N`: Set number of cores

**Outputs per sample:**
- `integrase_hits_summary.tsv` - Full table of candidate integrases
- `integrase_orfs.tsv` - Short list of integrase ORFs
- `trna.tsv` - ARAGORN-predicted tRNAs
- `integrase_trna.tsv` - Integrase-tRNA pairs
- `mge_query.fa` - Query fragments for BLAST
- `mge_blast.tsv` - BLAST hits linking tRNA and integrase
- `mge_region.fa` - Extracted MGE regions
- `mge_annotated.gbk` - Annotated MGE regions (GenBank)
- `attachment_sites.tsv` - Detected attL/attR coordinates

---

### 2. Bakta Pipeline

**Location:** `bakta_pipeline/`

**Purpose:** Annotates MGE regions using Bakta, then merges Bakta annotations with MGE-specific attachment site information.

**Prerequisites:**

```bash
conda create -n bakta -c bioconda -c conda-forge bakta
conda activate bakta
```

**Configuration:** `bakta_config.yaml`

Copy `bakta_config.yaml.example` to `bakta_config.yaml` and edit:

```yaml
paths:
  input_dir: ../finder_pipeline/results
  output_dir: ../bakta_pipeline/results_bakta
  db: /path/to/bakta/db
  genomes_dir: ../finder_pipeline/data/genomes

params:
  threads: 32
  prefix: annotation
  bakta_env: "bakta"
  merge_env: "MGE_finder"
```

**Run:**

```bash
cd bakta_pipeline
bash run_bakta.sh
```

**Outputs per sample:**
- `annotation.gbff` - Bakta-annotated GenBank file
- `annotation.faa` - Predicted protein sequences
- `merged.gbff` - Merged annotation with MGE attachment sites

---

### 3. PADLOC Pipeline

**Location:** `padloc_pipeline/`

**Purpose:** Detects defense systems (CRISPR, restriction-modification, etc.) in MGE regions using PADLOC.

**Prerequisites:**

```bash
conda create -n padloc -c bioconda padloc
conda activate padloc
padloc download_db  # Download PADLOC database
```

**Configuration:** `padloc_config.yaml`

Copy `padloc_config.yaml.example` to `padloc_config.yaml` and edit:

```yaml
paths:
  input_dir: ../finder_pipeline/results
  output_dir: padloc_pipeline/results

execution:
  conda_env: padloc
```

**Run:**

```bash
cd padloc_pipeline
bash run_padloc.sh
```

**Outputs per sample:**
- `results_padloc.csv` - Detected defense systems

**Alternative: Parallel execution**

For faster processing of many samples, use the parallel script:

```bash
bash additionals/run_padloc_parallel.sh
```

This script processes samples in parallel using GNU parallel and includes CRISPRDetect and Infernal preprocessing.

---

### 4. MGE Clustering

**Location:** `classification_investigation/clusterMGE/`

**Purpose:** Clusters MGE sequences using MMseqs2 based on nucleotide similarity to identify related elements.

**Configuration:** `mge_cluster_config.yaml`

Copy `mge_cluster_config.yaml.example` to `mge_cluster_config.yaml` and edit:

```yaml
paths:
  mge_gbk_pattern: "../finder_pipeline/results/*/mge_annotated.gbk"

execution:
  conda_env: "MGE_finder"

mmseqs:
  results_dir: "results_mge"
  min_seq_id: 0.9      # ≥90% identity
  coverage: 0.95       # ≥95% coverage
```

**Run:**

```bash
cd classification_investigation/clusterMGE
bash run_mge_cluster.sh
```

**Outputs:**
- `results_mge/mge_repdb.fasta` - Representative MGE sequences
- `results_mge/mge_cluster.tsv` - Cluster assignments

---

### 5. vContact Classification

**Location:** `classification_investigation/vcontact/`

**Purpose:** Classifies MGEs using vContact2/vContact3 based on protein content similarity.

**Prerequisites:**

```bash
conda create -n vContact3 -c bioconda vcontact2
# or
conda create -n vContact3 -c bioconda vcontact3
```

**Configuration:** `vcontact_config.yaml`

Copy `vcontact_config.yaml.example` to `vcontact_config.yaml` and edit:

```yaml
paths:
  rep_mge_fasta: "../clusterMGE/results_mge/rep_mge_nt.fa"
  output_dir: "results_vcontact"

execution:
  conda_env: "vContact3"
```

**Run:**

```bash
cd classification_investigation/vcontact
bash run_vcontact.sh
```

**Note:** If Snakemake complains about conda environments, ensure `--conda-prefix` is set correctly (see script).

---

### 6. Protein Clustering

**Location:** `protein_clusterization/`

**Purpose:** Clusters protein sequences from Bakta annotations using MMseqs2 to identify protein families.

**Configuration:** `classification_config.yaml`

Copy `classification_config.yaml.example` to `classification_config.yaml` and edit:

```yaml
conda_env: "MGE_finder"

proteins:
  pattern: "../bakta_pipeline/results_bakta/*/annotation.faa"

mmseqs:
  results_dir: "results"
  min_seq_id: 0.8
  coverage: 0.5
```

**Run:**

```bash
cd protein_clusterization
bash run_clusterization.sh
```

**Outputs:**
- `results/mmseqs_repdb.fasta` - Representative protein sequences
- `results/clusters_rep_seqs/` - Per-cluster FASTA files

---

### 7. RefSeq Database

**Location:** `refseq_prokaryote/`

**Purpose:** Download scripts for RefSeq prokaryote genomes.

**Download genomes:**

```bash
cd refseq_prokaryote
bash download_genomes_only.sh
```

This downloads only genomic FASTA files (`.fna.gz`) from RefSeq for bacteria and archaea. The download runs in the background and can be monitored via `download_genomes.log`.

**Monitor progress:**

```bash
tail -f download_genomes.log
find bacteria/ archaea/ -name "*.fna.gz" | wc -l
```

**Stop download:**

```bash
pkill -f download_genomes_only.sh
killall wget
```

---

##  Additional Scripts

### Statistics Collection

**Location:** `finder_pipeline/scripts/`

- `collect_mge_statistics.py` - Collects MGE statistics across all samples
- `defense_summary.py` - Summarizes PADLOC defense system detections

**Usage:**

```bash
# Collect MGE statistics
python finder_pipeline/scripts/collect_mge_statistics.py \
    --results_dir finder_pipeline/results \
    --output mge_summary.tsv

# Summarize defense systems
python finder_pipeline/scripts/defense_summary.py \
    --padloc_dir padloc_pipeline/results \
    --out_per_sample defense_per_sample.tsv \
    --out_overall defense_overall.tsv
```

### Analysis Scripts

**Location:** `analysis/`

- `stats.py` - General statistics collection
- `plot_breakpoint_distributions.py` - Plotting breakpoint distributions
- `prepare_integrases_fasta.py` - Prepare integrase sequences for analysis

---

## Output Files

### Finder Pipeline Outputs

Per sample in `finder_pipeline/results/{sample}/`:

- `orfs.gff`, `orfs.ffn`, `orfs.faa` - ORF predictions
- `integrase_hits.txt`, `integrase_hits_summary.tsv`, `integrase_orfs.tsv` - HMM results
- `trna.tsv` - ARAGORN-predicted tRNAs
- `integrase_trna.tsv` - Nearby integrase-tRNA pairs
- `mge_query.fa` - Query fragments for BLAST
- `mge_blast.tsv` - BLAST hits linking tRNA and integrase
- `mge_region.fa` - Extracted MGE regions
- `mge_annotated.gbk` - Annotated MGE regions (GenBank)
- `attachment_sites.tsv` - Detected attL/attR coordinates

### Summary Files (Root Directory)

- `mge_summary.tsv` - Overall MGE statistics
- `mge_overall.tsv` - Overall MGE counts
- `defense_per_sample.tsv` - Defense systems per sample
- `defense_overall.tsv` - Overall defense system frequencies

---

## Configuration

### Customizing Distance Thresholds

In `finder_pipeline/Snakefile`, modify the `trna_proximity` rule:

```python
rule trna_proximity:
    params:
        max_distance=500  # Change this value
```

### Customizing BLAST Window Size

In `finder_pipeline/scripts/annotate_mge_region.py`, modify:

```python
WINDOW_SIZE = 100000  # Change this value
```

### Re-generating Combined HMM

If you add/remove Pfam profiles, edit `finder_config.yaml` and rerun:

```bash
cd finder_pipeline
snakemake build_combined_hmm
```

---

##  Requirements

### Core Dependencies

- Python 3.10+
- Snakemake ≥ 7.x
- Biopython
- Prodigal or Pyrodigal
- ARAGORN
- HMMER 3.x
- BLAST+
- BCBio.GFF
- MMseqs2 (for clustering modules)

### Module-Specific Dependencies

- **Bakta**: Bakta annotation tool
- **PADLOC**: PADLOC defense system detector
- **vContact**: vContact2 or vContact3

---

##  Troubleshooting

### Snakemake Lock Issues

```bash
snakemake --unlock
```

### Conda Environment Issues

Ensure conda environments are created and activated before running pipelines. Some modules require specific conda environments (see module sections).

### Missing Input Files

Check that input paths in configuration files are correct and files exist.

### Log Files

Check log files in `logs/` directories for detailed error messages.

---

##  Additional Resources

- **Finder Pipeline Details**: See `finder_pipeline/README_finder.md`
- **vContact Notes**: See `classification_investigation/vcontact/README.md`

---

##  Workflow Order

Recommended execution order:

1. **Finder Pipeline** - Detect and annotate MGEs
2. **Bakta Pipeline** - Annotate MGE regions with Bakta
3. **PADLOC Pipeline** - Detect defense systems
4. **MGE Clustering** - Cluster MGE sequences
5. **vContact** - Classify MGEs (requires clustered representatives)
6. **Protein Clustering** - Cluster proteins from Bakta annotations

---

##  License


*Under condtruction*
---

##  Contributors

*Under constructruction*

---

##  References

- Prodigal: https://github.com/hyattpd/Prodigal
- HMMER: http://hmmer.org/
- ARAGORN: http://130.235.244.92/ARAGORN/
- Bakta: https://github.com/oschwengers/bakta
- PADLOC: https://github.com/padlocbio/padloc
- vContact: https://bitbucket.org/MAVERICLab/vcontact2
- MMseqs2: https://github.com/soedinglab/MMseqs2
