# How to monitor NCBI genome download progress

Working data directory: `IE_finder/data/`.

- **Genomes** are written to `ncbi/<taxon>/` (see `scripts/download_genomes.sh` / `download_genomes.py`).
- **Logs** go to `logs/download_progress.log` and optionally `logs/download_output.log` (created when the Python script runs).

## Quick commands

From `IE_finder/data`:

```bash
./scripts/monitor_progress.sh
```

Entrez download log:

```bash
tail -f logs/download_progress.log
```

Count `.fna.gz` for Deinococcales:

```bash
find ncbi/Deinococcales -name "*.fna.gz" 2>/dev/null | wc -l
```

## tmux (manual long runs)

```bash
cd IE_finder/data
tmux new-session -d -s genome_download "bash -lc 'export NCBI_EMAIL=... && ./scripts/download_genomes.sh'"
tmux attach -t genome_download
```
