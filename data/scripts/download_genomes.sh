#!/bin/bash
# Download genomes from NCBI into data/ncbi/<taxon>/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
NCBI_DIR="$DATA_ROOT/ncbi"
LOGS_DIR="$DATA_ROOT/logs"
mkdir -p "$NCBI_DIR" "$LOGS_DIR"
cd "$NCBI_DIR"

echo "=========================================="
echo "Downloading genomes from NCBI -> $NCBI_DIR"
echo "=========================================="

if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "Activating conda environment IE_finder..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate IE_finder
fi

if command -v datasets &> /dev/null; then
    echo "✓ Found tool: datasets (NCBI)"
    echo ""
    for taxon in "Deinococcus" "Thermococcus" "Pyrococcus" "Thermoplasma" "Picrophilus"; do
        echo ""
        echo "=========================================="
        echo "Downloading genomes for: $taxon"
        echo "=========================================="

        taxon_dir="${taxon// /_}"
        mkdir -p "$taxon_dir"

        datasets download genome taxon "$taxon" \
            --assembly-level complete \
            --include-gff3 \
            --filename "${taxon_dir}_genomes.zip"

        if [ -f "${taxon_dir}_genomes.zip" ]; then
            echo "Extracting archive..."
            unzip -q -d "$taxon_dir" "${taxon_dir}_genomes.zip"
            echo "✓ Genomes for $taxon downloaded and extracted"
        fi
    done

elif command -v ncbi-genome-download &> /dev/null; then
    echo "✓ Found tool: ncbi-genome-download"
    echo ""
    for taxon in "Deinococcus" "Thermococcus" "Pyrococcus" "Thermoplasma" "Picrophilus"; do
        echo ""
        echo "=========================================="
        echo "Downloading genomes for: $taxon"
        echo "=========================================="

        taxon_dir="${taxon// /_}"
        mkdir -p "$taxon_dir"

        ncbi-genome-download \
            --taxon "$taxon" \
            --assembly-level complete \
            --format genbank,fasta \
            --output-folder "$taxon_dir" \
            --parallel 4 \
            "$taxon"
    done

else
    echo "⚠ Specialized tools not found"
    echo "Falling back to Python Entrez script..."
    echo ""
    python3 "$SCRIPT_DIR/download_genomes.py"
fi

echo ""
echo "=========================================="
echo "Download complete!"
echo "Genomes in: $NCBI_DIR"
echo "=========================================="
