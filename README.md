# IE_finder

IE_finder is a Snakemake-based workflow for discovering **integrative elements (IEs)** in prokaryotic genomes. The repository also contains downstream analysis modules (Bakta annotation, defense-system screening, clustering, classification) built around the core IE discovery pipeline.

**Core pipeline documentation:** [IE_finder/README.md](IE_finder/README.md)

**Project map (directories, notebooks, archived notes):** see [docs/NAVIGATION.md](docs/NAVIGATION.md).

---

## Table of Contents

- [Project navigation (`docs/NAVIGATION.md`)](docs/NAVIGATION.md)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Modules Overview](#modules-overview)
  - [1. IE_finder Pipeline](#1-ie_finder-pipeline)
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
MGE_finder/                   # repository root (GitHub clone directory)
├── docs/
├── IE_finder/                # core IE discovery pipeline
│   ├── Snakefile, ie_finder_config.yaml, run.sh
│   ├── scripts/
│   ├── data/                 # input genomes (git-ignored)
│   └── results/              # pipeline outputs (git-ignored)
├── bakta_pipeline/
├── padloc_pipeline/
├── results_organized/
├── classification_investigation/
├── protein_clusterization/
├── data/
├── pfam/
├── envs/IE_finder.yaml
└── README.md, SETUP.md
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
mamba env create -f envs/IE_finder.yaml
conda activate IE_finder
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
# IE_finder pipeline
cd IE_finder
cp ie_finder_config.yaml.example ie_finder_config.yaml
# Edit ie_finder_config.yaml with your paths

# Other pipelines similarly
cd ../bakta_pipeline
cp bakta_config.yaml.example bakta_config.yaml
# etc.
```

See [SETUP.md](SETUP.md) for detailed configuration instructions.

---

##  Modules Overview

### 1. IE_IE_finder pipeline

**Location:** `IE_finder/`

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

**Configuration:** `ie_finder_config.yaml`

Copy `ie_finder_config.yaml.example` to `ie_finder_config.yaml` and edit:

```yaml
paths:
  genomes_dir: "data/genomes"
  results_dir: "results"

execution:
  conda_env: "IE_finder"

input_sources:
  - "/path/to/genomes/dir1"
  - "/path/to/genomes/dir2"

pfam_profiles:
  - "pfam/PF00589.hmm"
  - "pfam/PF22022.hmm"
```

**Run:**

```bash
cd IE_finder
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

Copy `bakta_config.yaml.example` to `bakta_config.yaml` and edit paths. Directories under `paths` are **relative to `bakta_pipeline/`** unless you use an absolute path. The Bakta database can be set as `paths.db` or via `BAKTA_DB` (the environment variable wins if set).

```yaml
paths:
  input_dir: ../IE_finder/results_deinococcales
  output_dir: results_bakta_deinococcales
  genomes_dir: ../IE_finder/data/genomes
  db: ""

params:
  threads: 32
  prefix: annotation
  bakta_env: "bakta"
  merge_env: "IE_finder"
```

Thermaceae: use `bakta_config_thermaceae.yaml` (e.g. `CONFIG=bakta_config_thermaceae.yaml bash run_bakta.sh`).

**Run:**

```bash
cd bakta_pipeline
export BAKTA_DB=/path/to/bakta_db   # if paths.db is empty
bash run_bakta.sh
```

**Outputs per sample:**
- `annotation.gbff` - Bakta-annotated GenBank file
- `annotation.faa` - Predicted protein sequences
- `annotation.gbff.merged` - Merged annotation with MGE attachment sites (`merge_att` rule)

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

Copy `padloc_config.yaml.example` to `padloc_config.yaml` and edit. Paths under `paths` are **relative to `padloc_pipeline/`** unless absolute.

```yaml
paths:
  input_dir: ../IE_finder/results_deinococcales
  output_dir: results_deinococcales
  logger_dir: ../IE_finder/scripts
  run_padloc_script: scripts/padloc_wrapper.py

execution:
  conda_env: padloc
```

Thermaceae / whole-genome variants: `padloc_config_thermaceae.yaml`, `padloc_config_thermaceae_wholegenome.yaml` (e.g. `CONFIG=padloc_config_thermaceae.yaml bash run_padloc.sh`).

**Run:**

```bash
cd padloc_pipeline
bash run_padloc.sh
```

**Outputs per sample:**
- `{basename_of_output_dir}_padloc.csv` (e.g. `results_deinococcales_padloc.csv`) — detected defense systems

---

### 4. MGE Clustering

**Location:** `classification_investigation/clusterMGE/`

**Purpose:** Clusters MGE sequences using MMseqs2 based on nucleotide similarity to identify related elements.

**Configuration:** `mge_cluster_config.yaml`

Copy `mge_cluster_config.yaml.example` to `mge_cluster_config.yaml` and edit. Glob patterns are **relative to `classification_investigation/clusterMGE/`** (or use absolute paths). Default inputs point at **`IE_finder/results_deinococcales`**.

MMseqs2: Snakemake uses the `mmseqs` executable on `PATH` inside the conda env. To force a binary, set **`MMSEQS_BIN`**.

```yaml
paths:
  mge_fasta_pattern: "../../IE_finder/results_deinococcales/*/mge_region.fa"
  mge_gbk_pattern: "../../IE_finder/results_deinococcales/*/mge_annotated.gbk"
  use_fasta: true

execution:
  conda_env: "IE_finder"

mmseqs:
  results_dir: "results_mge"
  min_seq_id: 0.9
  coverage: 0.95
```

For the merged **`results_organized`** layout, use `mge_cluster_config_organized.yaml` (e.g. `CONFIG=mge_cluster_config_organized.yaml bash run_mge_cluster.sh`).

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

Copy `vcontact_config.yaml.example` to `vcontact_config.yaml` and edit. Paths under `paths` are **relative to `classification_investigation/vcontact/`** unless absolute.

Set **`paths.vcontact_db_path`** or **`export VCONTACT_DB`** before running (environment variable is checked first in the Snakefile, then the YAML value).

```yaml
paths:
  rep_mge_fasta: "../clusterMGE/results_mge/unique_mge_for_vcontact.fa"
  output_dir: "results_vcontact"
  vcontact_db_path: "/path/to/vcontactdb"

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

Copy `classification_config.yaml.example` to `classification_config.yaml` and edit. The `proteins.pattern` glob is **relative to `protein_clusterization/`** unless absolute. MMseqs2: use the `mmseqs` from the conda env, or set **`MMSEQS_BIN`** to a specific binary (as in the MGE clustering Snakemake workflow).

```yaml
conda_env: "IE_finder"

proteins:
  pattern: "../bakta_pipeline/results_bakta_deinococcales/*/annotation.faa"

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

**Location:** `IE_finder/scripts/`

- `collect_mge_statistics.py` — агрегирует статистику MGE по всем сэмплам в каталоге результатов. Записывает три файла: `{prefix}_summary.tsv`, `{prefix}_details.tsv`, `{prefix}_overall.tsv`.

**Usage:**

```bash
python IE_finder/scripts/collect_mge_statistics.py \
    --results IE_finder/results \
    --out-prefix results_organized/summary/mge
```

Сводные таблицы PADLOC / defense (`defense_per_sample.tsv`, `defense_overall.tsv` и объединённые отчёты) хранятся рядом в [`results_organized/summary/`](results_organized/summary/). 
### Analysis Scripts

**Location:** `analysis/`

- `stats.py` - General statistics collection
- `plot_breakpoint_distributions.py` - Plotting breakpoint distributions
- `prepare_integrases_fasta.py` - Prepare integrase sequences for analysis

---

## Output Files

### IE_IE_finder pipeline Outputs

Per sample in `IE_finder/results/{sample}/`:

- `orfs.gff`, `orfs.ffn`, `orfs.faa` - ORF predictions
- `integrase_hits.txt`, `integrase_hits_summary.tsv`, `integrase_orfs.tsv` - HMM results
- `trna.tsv` - ARAGORN-predicted tRNAs
- `integrase_trna.tsv` - Nearby integrase-tRNA pairs
- `mge_query.fa` - Query fragments for BLAST
- `mge_blast.tsv` - BLAST hits linking tRNA and integrase
- `mge_region.fa` - Extracted MGE regions
- `mge_annotated.gbk` - Annotated MGE regions (GenBank)
- `attachment_sites.tsv` - Detected attL/attR coordinates

### Summary tables (repository default location)

Сводные TSV лежат в **`results_organized/summary/`** (их не держим в корне репозитория):

| File | Role |
|------|------|
| `mge_summary.tsv` | По-геномная сводка элементов (статистика пайплайна) |
| `mge_details.tsv` | Строка на каждый MGE |
| `mge_overall.tsv` | Агрегированные итоги |
| `mge_per_genome_stats.tsv` | Доп. агрегация по геномам (в т.ч. из аналитических скриптов) |
| `defense_per_sample.tsv`, `defense_overall.tsv` | Сводки по системам защиты (PADLOC) |
| `padloc_combined_results.tsv` | Объединённые результаты PADLOC |
| `mge_extended_stats.txt` | Текстовая расширенная сводка (скрипт `extend_mge_statistics.py`) |

---

## Configuration

### Customizing Distance Thresholds

In `IE_finder/Snakefile`, modify the `trna_proximity` rule:

```python
rule trna_proximity:
    params:
        max_distance=500  # Change this value
```

### Customizing BLAST Window Size

In `IE_finder/scripts/annotate_mge_region.py`, modify:

```python
WINDOW_SIZE = 100000  # Change this value
```

### Re-generating Combined HMM

If you add/remove Pfam profiles, edit `ie_finder_config.yaml` and rerun:

```bash
cd IE_finder
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

During runs, logs accumulate under each pipeline’s `logs/` folder and the repo root `logs/`. These directories are gitignored and safe to delete or rotate when debugging is done; see [docs/NAVIGATION.md](docs/NAVIGATION.md).

---

##  Additional Resources

- **Directory map (RU)**: [docs/NAVIGATION.md](docs/NAVIGATION.md)
- **IE_IE_finder pipeline Details**: See `IE_finder/README.md`
- **vContact Notes**: See `classification_investigation/vcontact/README.md`

---

##  Workflow Order

Recommended execution order:

1. **IE_IE_finder pipeline** - Detect and annotate MGEs
2. **Bakta Pipeline** - Annotate MGE regions with Bakta
3. **PADLOC Pipeline** - Detect defense systems
4. **MGE Clustering** - Cluster MGE sequences
5. **vContact** - Classify MGEs (requires clustered representatives)
6. **Protein Clustering** - Cluster proteins from Bakta annotations

---

##  License

*Under construction*

---

##  Contributors

*Under construction*

---

##  References

- Prodigal: https://github.com/hyattpd/Prodigal
- HMMER: http://hmmer.org/
- ARAGORN: http://130.235.244.92/ARAGORN/
- Bakta: https://github.com/oschwengers/bakta
- PADLOC: https://github.com/padlocbio/padloc
- vContact: https://bitbucket.org/MAVERICLab/vcontact2
- MMseqs2: https://github.com/soedinglab/MMseqs2
