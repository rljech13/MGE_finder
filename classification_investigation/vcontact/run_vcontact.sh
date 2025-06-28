#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_vcontact.sh ─ запуск Snakemake-конвейера vContact2 с «правильными» env-ами
###############################################################################


CONFIG="vcontact_config.yaml"
SNAKEFILE="Snakefile"
CORES=$(nproc)         

echo -e "\033[0;34m[INFO] Запуск с конфигом: ${CONFIG}\033[0m"

snakemake \
  --snakefile  "$SNAKEFILE" \
  --configfile "$CONFIG" \
  --cores      "$CORES" \
  --use-conda \
  --conda-prefix /home/lam34/miniforge3/envs \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
  "$@"        # передаём все флаги пользователя как есть