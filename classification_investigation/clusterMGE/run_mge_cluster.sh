#!/usr/bin/env bash
set -euo pipefail

WD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$WD"        # → /classification_investigation/mge_cluster

echo -e "\033[0;34m[INFO] MGE clustering pipeline\033[0m"
echo -e "\033[0;34m[INFO] Working dir: $WD\033[0m"

CORES=$(nproc)
snakemake                                                 \
  --snakefile Snakefile                                   \
  --configfile mge_cluster_config.yaml                    \
  --cores "$CORES"                                        \
  --use-conda                                             \
  --printshellcmds                                        \
  --show-failed-logs                                       \
  --rerun-incomplete "$@"