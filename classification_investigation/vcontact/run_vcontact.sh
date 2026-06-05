#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [ -z "$CONDA_BASE" ] || [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
fi
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx IE_finder; then
  conda activate IE_finder
elif conda env list | awk '{print $1}' | grep -qx snakemake; then
  conda activate snakemake
else
  echo "ERROR: Не найдено окружение с snakemake (IE_finder или snakemake)"
  exit 1
fi

CONFIG="${CONFIG:-vcontact_config.yaml}"
SNAKEFILE="Snakefile"
CORES=$(nproc)
CONDAP="${SNAKEMAKE_CONDA_PREFIX:-${MINIFORGE_ENVS:-$CONDA_BASE/envs}}"

echo -e "\033[0;34m[INFO] Запуск с конфигом: ${CONFIG}\033[0m"
echo -e "\033[0;34m[INFO] Используемое окружение: $(conda info --envs | grep '*' | awk '{print $1}')\033[0m"

snakemake \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIG" \
  --cores "$CORES" \
  --use-conda \
  --conda-prefix "$CONDAP" \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
  "$@"
