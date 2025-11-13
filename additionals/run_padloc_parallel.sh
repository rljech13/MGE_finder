#!/bin/bash

set -e

BASE_DIR="/home/lam34/MGE_finder"
SAMPLES_DIR="$BASE_DIR/results"
CRISPR_OUTDIR="$BASE_DIR/crisprdetect_results"
INFERNAL_OUTDIR="$BASE_DIR/infernal_results"
PADLOC_OUTDIR="$BASE_DIR/padloc_results"
PADLOC_DB="$HOME/.local/share/padloc"
THREADS=64 

process_sample() {
    sample_dir="$1"
    sample=$(basename "$sample_dir")
    input_fasta="$sample_dir/mge_region.fa"

    echo "📂 [$sample] File checkup $input_fasta..."

    if [[ ! -s "$input_fasta" ]]; then
        echo "[$sample] Skipping — file empty or missing."
        return
    fi

    mkdir -p "$CRISPR_OUTDIR" "$INFERNAL_OUTDIR" "$PADLOC_OUTDIR"

    echo "[$sample] Launch CRISPRDetect..."
    run-crisprdetect \
        --input "$input_fasta" \
        --output "$CRISPR_OUTDIR/${sample}_crispr"
    echo "[$sample] CRISPRDetect finished"

    crispr_gff="${CRISPR_OUTDIR}/${sample}_crispr.gff"

    echo "[$sample] Launch Infernal..."
    run-infernal \
        --input "$input_fasta" \
        --output "$INFERNAL_OUTDIR/${sample}_ncrna.tblout"
    echo "[$sample] Infernal finished"

    ncrna_tblout="${INFERNAL_OUTDIR}/${sample}_ncrna.tblout.formatted"

    echo "[$sample] Launch PADLOC..."
    padloc run \
        --db-dir "$PADLOC_DB" \
        --fna "$input_fasta" \
        --crispr "$crispr_gff" \
        --ncrna "$ncrna_tblout" \
        --output-dir "$PADLOC_OUTDIR/$sample"
    echo "[$sample] PADLOC finised"
    echo "[$sample] Sample completed"
    echo "--------------------------------------------------"
}

export -f process_sample

echo "Starting parallel with $THREADS threads..."
find "$SAMPLES_DIR"/*/ -type d | parallel -j "$THREADS" process_sample {}
echo "All samples completed."