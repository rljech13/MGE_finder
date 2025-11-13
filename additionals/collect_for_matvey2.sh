REP_FASTA="/home/lam34/MGE_finder/classification_investigation/clusterMGE/results_mge/rep_mge_nt.fa"
RES_DIR="/home/lam34/MGE_finder/bakta_pipeline/results_bakta"
OUT_DIR="/home/lam34/MGE_finder/analysis/selected_samples"

mkdir -p "$OUT_DIR"

grep '^>' "$REP_FASTA" \
| awk -F'sample=' '{split($2,a,/[\|[:space:]]/); print a[1]}' \
| sort -u \
| while read -r S; do
    d="$RES_DIR/$S"
    for f in annotation.gbff.merged; do
        if [[ -s "$d/$f" ]]; then
            cp "$d/$f" "$OUT_DIR/${S}_$f"
        else
            echo "WARN: нет файла $d/$f" >&2
        fi
    done
done