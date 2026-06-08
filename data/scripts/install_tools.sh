#!/bin/bash
# Install NCBI genome download tools

set -e

echo "=========================================="
echo "Installing genome download tools"
echo "=========================================="

if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install conda/miniconda first."
    exit 1
fi

if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "Activating conda environment IE_finder..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate IE_finder
fi

echo ""
echo "Choose a tool to install:"
echo "1) datasets (NCBI, recommended — fast and modern)"
echo "2) ncbi-genome-download (alternative, via pip)"
echo "3) Both tools"
echo ""
read -p "Enter choice (1-3): " choice

case $choice in
    1)
        echo "Installing datasets..."
        conda install -c conda-forge ncbi-datasets-cli -y
        echo "✓ datasets installed"
        ;;
    2)
        echo "Installing ncbi-genome-download..."
        pip install ncbi-genome-download
        echo "✓ ncbi-genome-download installed"
        ;;
    3)
        echo "Installing both tools..."
        conda install -c conda-forge ncbi-datasets-cli -y
        pip install ncbi-genome-download
        echo "✓ Both tools installed"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Installation complete!"
D="$(cd "$(dirname "$0")/.." && pwd)"
echo "Next: cd \"$D\" && ./scripts/download_genomes.sh"
echo "=========================================="
