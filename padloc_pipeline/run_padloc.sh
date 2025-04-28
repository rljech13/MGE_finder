#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
CONFIG="padloc_config.yaml"
SNAKEFILE="Snakefile"
CORES=$(nproc)
WORKDIR=$(pwd)

# === COLORS ===
BLUE="\033[0;34m"
RESET="\033[0m"

echo -e "${BLUE}[INFO] Using config: ${CONFIG}${RESET}"
echo -e "${BLUE}[INFO] Running Padloc Snakemake pipeline in: ${WORKDIR}${RESET}"

# === PARSE FLAGS ===
EXTRA_ARGS=()
for arg in "$@"; do
  case $arg in
    --dry-run|-n) EXTRA_ARGS+=("--dry-run") ;;
    --unlock) EXTRA_ARGS+=("--unlock") ;;
    --force|-f) EXTRA_ARGS+=("--force") ;;
    --rerun-incomplete) EXTRA_ARGS+=("--rerun-incomplete") ;;
    --cores=*|-j=*) EXTRA_ARGS+=("$arg") ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

# === RUN SNAKEMAKE ===
snakemake \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIG" \
  --cores "$CORES" \
  --use-conda \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
  --keep-going \
  "${EXTRA_ARGS[@]}"