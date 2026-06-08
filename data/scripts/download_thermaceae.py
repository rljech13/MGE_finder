#!/usr/bin/env python3
"""Download Thermaceae genomes from NCBI (output directory set via environment)."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from download_genomes import download_genomes_via_entrez

OUTPUT_DIR = Path(os.environ.get("THERMACEAE_DOWNLOAD_ROOT", "/mnt/data/procaryota_genomes/ncbi"))

def main():
    taxa = ["Thermaceae"]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log_file = OUTPUT_DIR / "download_thermaceae.log"
    log_file.unlink(missing_ok=True)

    download_genomes_via_entrez(taxa, OUTPUT_DIR, log_file=str(log_file))

if __name__ == "__main__":
    main()
