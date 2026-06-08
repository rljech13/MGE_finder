# NCBI genome downloads and metadata

Local layer for downloading assemblies from NCBI, storing FASTA files, and building taxonomy tables. The primary consumer in this repository is **Deinococcales** preparation for `finder_pipeline` (`data/ncbi/Deinococcales/*.fna.gz`).

Thermus accession supplements live outside this directory:  
`Thermaceae_genomes/thermus_missing_accessions/` (see [Thermus accession supplement](#thermus-accession-supplement)).

> **Project name:** the tool is branded **IE_finder**. The repository directory on disk may still be named `MGE_finder` on some machines; paths below use `IE_finder/` as the repository root.

---

## Directory layout

```
data/
├── README.md              # this file
├── ncbi/                  # downloaded genomes (*.fna.gz or ncbi_dataset/)
│   ├── Deinococcales/     # main set for the pipeline
│   ├── Deinococcus/
│   ├── Thermococcus/
│   └── ...
├── metadata/              # Deinococcales taxonomy TSV/MD
├── scripts/               # utilities (see table below)
├── docs/                  # historical notes (STATUS, TAXONOMY_REPORT, …)
└── logs/                  # download_progress.log (not committed)
```

Downloaded `*.fna.gz` files and `logs/` are excluded by the root `.gitignore`; re-download them locally with the scripts below.

---

## Requirements

### Conda environment

```bash
source ~/miniforge3/etc/profile.d/conda.sh   # adjust to your conda install
conda activate IE_finder
```

Required packages: `biopython`, `wget`. Optional: `ncbi-datasets-cli` (`datasets`).

### NCBI CLI tools

Interactive installer:

```bash
cd IE_finder/data
./scripts/install_tools.sh
```

Options: `datasets` (recommended), `ncbi-genome-download`, or both.

### Entrez email (required for Python scripts)

```bash
export NCBI_EMAIL='your.email@example.com'
```

Without a valid email, NCBI may throttle or block requests.

---

## Which script to use

| Task | Script | Tool | Output |
|------|--------|------|--------|
| Full **Deinococcales** order (bulk, resumable) | `download_genomes.py` | Entrez + wget via FTP | `ncbi/Deinococcales/*.fna.gz` |
| Several **genera** (complete genomes, zip) | `download_genomes.sh` | `datasets` or `ncbi-genome-download` | `ncbi/<Taxon>/` |
| Full **Thermaceae** family (outside `data/`) | `download_thermaceae.py` | Entrez + wget | `$THERMACEAE_DOWNLOAD_ROOT` |
| **Accession list** (Thermus supplement) | `Thermaceae_genomes/.../run_download.sh` | `datasets download genome accession` | supplement zip → FASTA |
| Metadata for downloaded GCAs | `extract_taxonomy_metadata_fixed.py` | Entrez | `metadata/deinococcales_metadata_complete.tsv` |
| Metadata summary tables | `create_taxonomy_summary.py` | local | `metadata/deinococcales_taxonomy_*` |
| Progress snapshot | `monitor_progress.sh` | — | stdout |
| Pre-download taxonomy check | `check_taxonomy.py` | Entrez | stdout |
| Missing-metadata diagnostics | `check_missing_metadata.py` | Entrez | stdout |

**Practical recommendation:** for large downloads with resume after network drops, use **`download_genomes.py`** (skips existing files on the Entrez path). For quick accession batches, use **`datasets`** + `run_download.sh`.

---

## Scripts in `scripts/` (details)

### `download_genomes.py`

Downloads genomes for the **Deinococcales** order (all genera within the order).

**Tool priority inside `main()`:**

1. If `datasets` is available → `download_with_datasets` (taxon, complete, zip).
2. Else if `ncbi-genome-download` is available → genbank/fasta folder.
3. Else → **Entrez + wget** (tested fallback; ~467 `.fna.gz` in `ncbi/Deinococcales`).

**Run:**

```bash
cd IE_finder/data
export NCBI_EMAIL='your.email@example.com'
python3 scripts/download_genomes.py
```

**Output:**

- FASTA: `ncbi/Deinococcales/<assembly>_genomic.fna.gz`
- Log: `logs/download_progress.log`

**Re-run:** existing `.fna.gz` files are **skipped** on the Entrez path — safe to restart after a network failure.

**Note:** when `datasets` is installed, `main()` uses it instead of Entrez; the skip-existing-files resume behaviour applies only to the Entrez fallback.

**Long run in tmux:**

```bash
cd IE_finder/data
tmux new-session -d -s deinococcales_download \
  "bash -lc 'export NCBI_EMAIL=your.email@example.com && python3 scripts/download_genomes.py 2>&1 | tee logs/download_output.log'"
tmux attach -t deinococcales_download
```

---

### `download_genomes.sh`

Bash wrapper: downloads **by genus** via `datasets` (or `ncbi-genome-download`):

`Deinococcus`, `Thermococcus`, `Pyrococcus`, `Thermoplasma`, `Picrophilus` — assembly level `complete`, GFF3 included in zip.

**Run:**

```bash
cd IE_finder/data
./scripts/download_genomes.sh
```

Activates `IE_finder` if conda is not already in PATH. Result: `ncbi/<Taxon>/` with extracted `ncbi_dataset/`.

**tmux:**

```bash
tmux new-session -d -s genome_download \
  "bash -lc 'cd IE_finder/data && ./scripts/download_genomes.sh'"
```

---

### `download_thermaceae.py`

Thin wrapper around `download_genomes_via_entrez` for the **Thermaceae** taxon.

**Environment variables:**

```bash
export THERMACEAE_DOWNLOAD_ROOT=/path/to/output   # default: /mnt/data/procaryota_genomes/ncbi
export NCBI_EMAIL='your.email@example.com'
```

**Run:**

```bash
cd IE_finder/data
python3 scripts/download_thermaceae.py
```

Log: `$THERMACEAE_DOWNLOAD_ROOT/download_thermaceae.log`

---

### `extract_taxonomy_metadata_fixed.py`

For each `*.fna.gz` in `ncbi/Deinococcales/`, extracts the GCA, queries NCBI Assembly, and writes a full metadata table.

**Run** (after genome download):

```bash
cd IE_finder/data
export NCBI_EMAIL='your.email@example.com'
python3 scripts/extract_taxonomy_metadata_fixed.py
```

**Output:** `metadata/deinococcales_metadata_complete.tsv`

Takes ~0.34 s per genome (NCBI rate limit). For hundreds of files, use tmux.

---

### `create_taxonomy_summary.py`

Builds summary tables from `metadata/deinococcales_metadata_complete.tsv` (fallback: `deinococcales_metadata.tsv`).

```bash
cd IE_finder/data
python3 scripts/create_taxonomy_summary.py
```

**Output:**

- `metadata/deinococcales_taxonomy_summary.tsv`
- `metadata/deinococcales_taxonomy_detailed.md`

---

### `check_taxonomy.py`

Before bulk download: shows which genera/species will be included in **Deinococcales** (sample of up to 100 assemblies via Entrez).

```bash
cd IE_finder/data
export NCBI_EMAIL='your.email@example.com'
python3 scripts/check_taxonomy.py
```

Stdout only; does not modify files.

---

### `check_missing_metadata.py`

Compares files in `ncbi/Deinococcales/` with `metadata/*.tsv` and diagnoses the first 20 gaps via NCBI.

```bash
cd IE_finder/data
python3 scripts/check_missing_metadata.py
```

---

### `monitor_progress.sh`

Quick snapshot of Entrez download progress:

```bash
cd IE_finder/data
./scripts/monitor_progress.sh
# or in a loop:
watch -n 5 ./scripts/monitor_progress.sh
```

Shows the tail of `logs/download_progress.log` and `*.gz` counts per taxon in `ncbi/`.

---

### `install_tools.sh`

See [NCBI CLI tools](#ncbi-cli-tools).

---

## Thermus accession supplement

Directory: `Thermaceae_genomes/thermus_missing_accessions/` (sibling to the IE_finder repo on this machine).

| File | Purpose |
|------|---------|
| `missing_assembly_accessions.txt` | GCA/GCF list (one per line), 71 accessions |
| `missing_5_accessions.txt` | 5 recent records that `datasets` may not serve yet |
| `run_download.sh` | Main supplement download |
| `run_download_5.sh` | Retry for the remaining 5 |
| `download.log` / `download_5.log` | Run logs |

**Successful run (March 2026):** 66 assemblies in ~18 s →  
`Thermaceae_genomes/ncbi_thermus_supplement/ncbi_dataset/data/`

**Run:**

```bash
# background session so SSH drops do not kill the job
tmux new-session -d -s thermus_dl \
  "bash Thermaceae_genomes/thermus_missing_accessions/run_download.sh"
tmux attach -t thermus_dl
```

**Custom accession list:** put accessions in a text file (one per line) and edit `LIST=` in `run_download.sh`.

**List format:**

```
GCA_046097575.1
GCF_046293695.1
```

**Verify after download:**

```bash
ls Thermaceae_genomes/ncbi_thermus_supplement/ncbi_dataset/data/ | wc -l
tail Thermaceae_genomes/thermus_missing_accessions/download.log
```

---

## Bridge to `finder_pipeline`

NCBI downloads and pipeline input use different layouts:

| Stage | Location | Format |
|-------|----------|--------|
| Download (this directory) | `data/ncbi/Deinococcales/` | flat `GCA_*_genomic.fna.gz` |
| Pipeline input | `finder_pipeline/data/genomes/` | `{sample}.fna` (uncompressed) |
| Prepared Deinococcales set | `finder_pipeline/data/genomes_deinococcales/` | flat `.fna` (358 assemblies) |
| `datasets` layout | `ncbi_dataset/data/GCA_.../` | `GCA_*_genomic.fna` per directory |

`prepare_fastas.py` expects the NCBI directory layout (`*_genomic.fna`, sample name from parent folder). Point `input_sources` at a `datasets` tree or at `data/genomes_deinococcales/`.

**From scratch — Deinococcales for finder_pipeline:**

```bash
cd IE_finder/data
export NCBI_EMAIL='your.email@example.com'

# 1) Download (prefer tmux)
tmux new-session -d -s deinococcales_download \
  "bash -lc 'python3 scripts/download_genomes.py 2>&1 | tee logs/download_output.log'"

# 2) Monitor
tail -f logs/download_progress.log
./scripts/monitor_progress.sh

# 3) Metadata
python3 scripts/extract_taxonomy_metadata_fixed.py
python3 scripts/create_taxonomy_summary.py

# 4) Verify
find ncbi/Deinococcales -name '*.fna.gz' | wc -l
```

**Configure and run the pipeline:**

```bash
cd IE_finder/finder_pipeline
cp ie_finder_config.yaml.example ie_finder_config.yaml
# set input_sources to data/genomes_deinococcales or an ncbi_dataset path
./run.sh --cores=8
```

`DEINOCOCCALES_SOURCE_DIR` (if used by prepare scripts) defaults to `IE_finder/data/ncbi/Deinococcales`.

### B. Quick genus download (zip + NCBI metadata)

```bash
cd IE_finder/data
./scripts/download_genomes.sh
```

### C. Thermus supplement only

See [Thermus accession supplement](#thermus-accession-supplement). Merge with an existing set by **numeric assembly ID** (digits after `GCA_`/`GCF_`), not by GenBank/RefSeq prefix.

---

## Network failure protection

| Mechanism | Where | Effect |
|-----------|-------|--------|
| **tmux** | any long run | process survives SSH disconnect |
| **skip existing files** | `download_genomes.py` (Entrez path) | re-run downloads only missing genomes |
| **log to file** | `logs/download_progress.log`, `tee` in run_download | see where the run stopped |
| **rate limit 0.34 s** | Entrez in Python | fewer NCBI blocks |
| **keep zip archives** | `run_download.sh` | re-unzip with `-o` is idempotent |

**Limitation:** `download_genomes.py` uses wget **without** `-c`; on error, a partial `.fna.gz` is deleted. A re-run re-downloads that file in full (existing files are untouched).

---

## Environment variables

| Variable | Scripts | Default |
|----------|---------|---------|
| `NCBI_EMAIL` | all Entrez Python scripts | `your.email@example.com` (must be set!) |
| `THERMACEAE_DOWNLOAD_ROOT` | `download_thermaceae.py` | `/mnt/data/procaryota_genomes/ncbi` |
| `DEINOCOCCALES_SOURCE_DIR` | finder prepare scripts | `IE_finder/data/ncbi/Deinococcales` |

---

## Logs and monitoring

| File | Source |
|------|--------|
| `logs/download_progress.log` | `download_genomes.py` (Entrez) |
| `logs/download_output.log` | full stdout when run via `tee` |
| `Thermaceae_genomes/.../download.log` | `run_download.sh` |

```bash
# count downloaded files
find IE_finder/data/ncbi/Deinococcales -name '*.fna.gz' | wc -l
du -sh IE_finder/data/ncbi/Deinococcales
```

---

## Troubleshooting

**`datasets: command not found`**  
→ `./scripts/install_tools.sh` or `conda install -c conda-forge ncbi-datasets-cli`

**Empty Entrez response / few genomes**  
→ check `NCBI_EMAIL`; for Thermus, NCBI website “latest” counts may differ from `datasets summary taxon`.

**`datasets download` — `record_count: 0` for recent GCF**  
→ assembly exists in Assembly DB but the genome package is not in Datasets yet; wait or use FTP via Entrez.

**Resume after disconnect**  
→ run the same script again; check the log and `monitor_progress.sh`.

**Metadata missing for some GCAs**  
→ `python3 scripts/check_missing_metadata.py`, then re-run `extract_taxonomy_metadata_fixed.py`.

---

## Related files in `docs/`

Historical notes (optional reading):

- `docs/HOW_TO_MONITOR.md` — tmux and log tailing
- `docs/TAXONOMY_REPORT.md` — Deinococcales composition report
- `docs/STATUS.md` — snapshot of bulk download status

The canonical guide is **this README**.
