#!/bin/bash
# Monitor download progress (logs in data/logs/)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOGS_DIR="$DATA_ROOT/logs"
LOG_FILE="$LOGS_DIR/download_progress.log"
OUTPUT_LOG="$LOGS_DIR/download_output.log"
NCBI_DIR="$DATA_ROOT/ncbi"

mkdir -p "$LOGS_DIR"

echo "=========================================="
echo "Genome download monitor (NCBI -> $NCBI_DIR)"
echo "=========================================="
echo ""
echo "Progress (Python/Entrez): tail -f $LOG_FILE"
echo "Full output (if tee is used): tail -f $OUTPUT_LOG"
echo ""
echo "=========================================="
echo "Last 20 lines of progress log:"
echo "=========================================="
if [ -f "$LOG_FILE" ]; then
    tail -20 "$LOG_FILE"
else
    echo "Log file not created yet."
fi
echo ""
echo "=========================================="
echo "*.gz files per taxon in ncbi/:"
echo "=========================================="
for dir in Deinococcus Thermococcus Pyrococcus Thermoplasma Picrophilus Deinococcales; do
    d="$NCBI_DIR/$dir"
    if [ -d "$d" ]; then
        count=$(find "$d" -name "*.gz" 2>/dev/null | wc -l)
        size=$(du -sh "$d" 2>/dev/null | cut -f1)
        echo "  $dir: $count files, $size"
    fi
done
echo ""
echo "Refresh: watch -n 5 $0"
