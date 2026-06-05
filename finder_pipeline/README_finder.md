# finder_pipeline

Snakemake workflow that screens bacterial and archaeal assemblies for integrative mobile genetic elements (MGEs / IEs) inserting into tRNA genes.

Starting from genome FASTA files, the pipeline predicts ORFs, scans for integrases with Pfam HMMs, finds the nearest opposite-strand tRNA within 500 bp, extracts the putative element region via BLAST, and produces an annotated GenBank file with attL/attR sites.

## Layout

```
finder_pipeline/
├── Snakefile
├── finder_config.yaml          # active config (copy from .example)
├── finder_config.yaml.example
├── run.sh                      # main entry point
├── e2e_test/                   # stress-test panel + config
├── scripts/
│   ├── prepare_fastas.py       # ingest NCBI / Hybracter / bin FASTA
│   ├── sequence_qc.py          # N-ambiguity QC helpers
│   ├── build_combined_hmm.py   # merge Pfam profiles + hmmpress
│   ├── predict_orfs.py
│   ├── hmm_search.py
│   ├── annotate_trna_proximity.py
│   ├── extract_trna_region.py
│   ├── annotate_mge_region.py
│   ├── extract_mge_regions.py
│   ├── annotate_and_orient_mge.py
│   ├── filter_confident_ie.py  # confident IE QC (V_3pS + length filters)
│   ├── dedup_ie_representatives.py  # cohort MMseqs2 dedup
│   ├── v3ps_filters.py         # shared filter logic (library)
│   ├── collect_mge_statistics.py  # optional post-run summary
│   ├── batch_runner.py         # optional manifest batch runner
│   └── logger.py
├── data/genomes/               # canonical input (git-ignored)
└── results/                    # output (git-ignored)
```

## Quick start

```bash
# From repo root: create conda env
conda env create -f envs/MGE_finder.yaml
conda activate MGE_finder

cd finder_pipeline
cp finder_config.yaml.example finder_config.yaml
# Edit input_sources and paths in finder_config.yaml

./run.sh                  # full run
./run.sh --dry-run        # validate workflow graph
CONFIG=finder_config.yaml ./run.sh --cores=4
```

## Configuration

All paths in `finder_config.yaml` are relative to `finder_pipeline/` unless absolute.

| Key | Description |
|-----|-------------|
| `paths.genomes_dir` | Directory for canonical `*.fna` inputs |
| `paths.results_dir` | Snakemake output directory |
| `input_sources` | Source dirs searched recursively |
| `ncbi_pattern` | NCBI-style assembly glob (default `*_genomic.fna`) |
| `hybracter_pattern` | Hybracter assembly glob |
| `bins_pattern` | Metagenomic bin glob (default `*.fa`) |
| `reject_ambiguous_n` | Skip genomes with N before ORF prediction |
| `reject_ambiguous_n_ie` | Flag for downstream IE N-gap filtering |
| `pfam_profiles` | Integrase HMM files (bundled in `../pfam/`) |

## Filters and confident IE output

The pipeline produces two layers of output per sample:

| Layer | Files | Meaning |
|-------|-------|---------|
| Candidates | `mge_region.fa`, `mge_annotated.gbk` | All IE candidates after detection |
| **Confident set** | `ie_confident.fa`, `ie_confident.gbk`, `ie_filter_audit.tsv` | High-confidence IE after full QC |

### Filter chain (confident set)

Mirrors publication stages 4–8 (Stage 7_h + Stage 8 FINAL):

| Step | Criterion | Module |
|------|-----------|--------|
| V_3pS strict | 3'-anchor (qend ≥ trna_len − 3), strand-aware BLAST hit | `annotate_mge_region.py`, `v3ps_filters.py` |
| attL length | best attL ≥ 15 bp (longest strict hit from raw BLAST) | `filter_confident_ie.py` |
| Intergenic | attL does not overlap Prodigal CDS on full contig | `filter_confident_ie.py` |
| Integrase | ORF span > 300 aa | `filter_confident_ie.py` |
| IE length | extracted region ≥ 2000 nt | `filter_confident_ie.py` |
| N gaps | no ambiguous N in `mge_region.fa` | `filter_confident_ie.py` |

Thresholds are configured in `finder_config.yaml` under `filters:`.

Each rejected IE is logged in `ie_filter_audit.tsv` with `reject_reason` (e.g. `attl_too_short`, `cds_partial_CDS`, `integrase_too_short`).

## Representatives (MMseqs dedup, Stage 9)

After all per-sample confident IE are collected, `dedup_ie_representatives` runs MMseqs2 easy-cluster on the combined confident set:

| Output | Description |
|--------|-------------|
| `dedup/ie_representatives.fa` | Non-redundant IE sequences |
| `dedup/ie_representatives.gbk` | Representative GenBank records |
| `dedup/rep_index.tsv` | Index of per-rep `.gbk` files under `dedup/representatives/` |
| `dedup/ie_dedup_mapping.tsv` | member → representative mapping |
| `dedup/ie_cluster_metadata.tsv` | Cluster sizes and membership |

Default parameters (configurable via `dedup:` in config): min-seq-id 0.9, coverage 0.8 — same as publication Stage 9.

```bash
snakemake dedup_ie_representatives --configfile finder_config.yaml
```

## Workflow

| Rule | Script / tool | Output |
|------|---------------|--------|
| `prepare_fasta` | `prepare_fastas.py` | `data/genomes/*.fna` |
| `predict_orfs` | Prodigal | `orfs.gff`, `orfs.ffn`, `orfs.faa` |
| `build_combined_hmm` | `build_combined_hmm.py` | `combined/pfam_combined.hmm` |
| `hmm_search` | `hmm_search.py` + hmmscan | `integrase_hits_summary.tsv` |
| `predict_trna` | Aragorn | `trna.tsv` |
| `trna_proximity` | `annotate_trna_proximity.py` | `integrase_trna.tsv` |
| `extract_trna_region` | `extract_trna_region.py` | `mge_query.fa` |
| `blast_mge` | `annotate_mge_region.py` | `mge_blast.tsv` |
| `extract_mge_region` | `extract_mge_regions.py` | `mge_region.fa` |
| `annotate_mge` | `annotate_and_orient_mge.py` | `mge_annotated.gbk`, `attachment_sites.tsv` |
| `filter_confident_ie` | `filter_confident_ie.py` | **`ie_confident.fa`**, **`ie_confident.gbk`**, `ie_filter_audit.tsv` |
| `dedup_ie_representatives` | `dedup_ie_representatives.py` | **`dedup/ie_representatives.fa`**, `rep_index.tsv`, mapping |

Per-sample outputs live in `results/<sample>/`. **Primary deliverable:** `results/dedup/ie_representatives.fa` + `rep_index.tsv`.

## Filters implemented in the pipeline

See **Filters and confident IE output** above. Genome-level N-filter runs at ingest; all other QC runs in `filter_confident_ie`. Cohort deduplication runs in `dedup_ie_representatives`.

## Batch processing

For large manifest-driven runs:

```bash
python scripts/batch_runner.py \
  --manifest /path/to/genomes.txt \
  --work-dir . \
  --batch-size 100
```

## Statistics

```bash
python scripts/collect_mge_statistics.py \
  --results results/ \
  --out-prefix mge_stats
```

## Customisation

- **tRNA distance**: change `--max_distance` in the `trna_proximity` Snakefile rule (default 500).
- **BLAST window**: edit `WINDOW_SIZE` in `annotate_mge_region.py` (default 300 000 bp).
- **BLAST edge tolerance**: edit `SHIFT` in `annotate_mge_region.py` (default 3 bp).
- **Pfam profiles**: edit `pfam_profiles` in config, then rerun `build_combined_hmm`.

## Tests

```bash
cd /path/to/MGE_finder
python -m pytest tests/finder_pipeline/ -q
```
