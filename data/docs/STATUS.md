# Genome download status (historical note)

This file described the state during a bulk download run. Current directory layout:

- **`ncbi/`** — downloaded archives by taxon (`Deinococcales/`, `Deinococcus/`, …).
- **`metadata/`** — Deinococcales metadata tables and reports.
- **`scripts/`** — launch and check scripts.
- **`logs/`** — run logs (optional, not committed).

## Monitoring now

```bash
cd IE_finder/data
./scripts/monitor_progress.sh
tail -f logs/download_progress.log   # if a log is being written
```

## Example counts

```bash
cd IE_finder/data
find ncbi/Deinococcales -name "*.fna.gz" | wc -l
du -sh ncbi/Deinococcales
```
