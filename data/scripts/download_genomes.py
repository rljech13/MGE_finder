#!/usr/bin/env python3
"""
Download genomes from NCBI.
Primary target: Deinococcales order assemblies for IE_finder / finder_pipeline.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

try:
    from Bio import Entrez
    from Bio import SeqIO
except ImportError:
    print("Install biopython: conda install -c bioconda biopython")
    sys.exit(1)

Entrez.email = os.environ.get("NCBI_EMAIL", "your.email@example.com")
if Entrez.email == "your.email@example.com":
    print("WARNING: Set NCBI_EMAIL or change the email in this script")
    print("Example: export NCBI_EMAIL='your.email@example.com'")
    print("Continuing with default email...")

DATA_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = DATA_ROOT / "ncbi"
LOGS_DIR = DATA_ROOT / "logs"
NCBI_DATASETS_CMD = "datasets"
NCBI_GENOME_DOWNLOAD_CMD = "ncbi-genome-download"

def check_tool_available(tool_name):
    """Return True if the CLI tool responds to --version."""
    try:
        subprocess.run([tool_name, "--version"],
                      capture_output=True,
                      check=True,
                      timeout=5)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False

def download_with_datasets(taxa, output_dir, assembly_level="complete"):
    """Download genomes using NCBI datasets."""
    print(f"Using datasets to download: {taxa}")

    for taxon in taxa:
        print(f"\nDownloading genomes for {taxon}...")
        taxon_dir = output_dir / taxon.replace(" ", "_")
        taxon_dir.mkdir(exist_ok=True)

        zip_file = taxon_dir / f"{taxon.replace(' ', '_')}_genomes.zip"

        cmd = [
            NCBI_DATASETS_CMD, "download", "genome",
            "taxon", taxon,
            "--assembly-level", assembly_level,
            "--include-gff3",
            "--filename", str(zip_file)
        ]

        try:
            subprocess.run(cmd, check=True, cwd=str(taxon_dir))
            print(f"Genomes for {taxon} saved to {zip_file}")
            print(f"Extract archive: unzip {zip_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {taxon}: {e}")

def download_with_ncbi_genome_download(taxa, output_dir, assembly_level="complete"):
    """Download genomes using ncbi-genome-download."""
    print(f"Using ncbi-genome-download to download: {taxa}")

    for taxon in taxa:
        print(f"\nDownloading genomes for {taxon}...")
        taxon_dir = output_dir / taxon.replace(" ", "_")
        taxon_dir.mkdir(exist_ok=True)

        cmd = [
            NCBI_GENOME_DOWNLOAD_CMD,
            "--taxon", taxon,
            "--assembly-level", assembly_level,
            "--format", "genbank,fasta",
            "--output-folder", str(taxon_dir),
            "--parallel", "4"
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"Genomes for {taxon} saved to {taxon_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {taxon}: {e}")

def get_genome_accessions_via_entrez(taxon_name):
    """Return assembly IDs for a taxon via Entrez."""
    print(f"Searching genomes for {taxon_name}...")

    search_term = f'"{taxon_name}"[Organism] AND "latest"[filter]'

    try:
        handle = Entrez.esearch(db="assembly", term=search_term, retmax=10000)
        record = Entrez.read(handle)
        handle.close()

        count = int(record["Count"])
        print(f"Found {count} genomes for {taxon_name}")

        if count > 0:
            ids = record["IdList"]
            if count > 10000:
                all_ids = ids[:]
                retstart = 10000
                while retstart < count:
                    handle = Entrez.esearch(db="assembly", term=search_term, retmax=10000, retstart=retstart)
                    record = Entrez.read(handle)
                    handle.close()
                    all_ids.extend(record["IdList"])
                    retstart += 10000
                return all_ids
            return ids
        return []

    except Exception as e:
        print(f"Error searching genomes for {taxon_name}: {e}")
        return []

def download_genomes_via_entrez(taxa, output_dir, log_file=None):
    """Download genomes via Entrez summaries and wget from NCBI FTP."""
    def log_print(message):
        print(message)
        if log_file:
            with open(log_file, 'a', encoding='utf-8') as f:
                f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

    log_print("Using Entrez + wget for genome download...")

    total_taxa = len(taxa)
    for taxon_idx, taxon in enumerate(taxa, 1):
        log_print(f"\n{'='*60}")
        log_print(f"[{taxon_idx}/{total_taxa}] Processing {taxon}...")
        log_print(f"{'='*60}")

        taxon_dir = output_dir / taxon.replace(" ", "_")
        taxon_dir.mkdir(exist_ok=True)

        log_print(f"Searching genomes for {taxon}...")
        assembly_ids = get_genome_accessions_via_entrez(taxon)

        if not assembly_ids:
            log_print(f"⚠ No genomes found for {taxon}")
            continue

        total_genomes = len(assembly_ids)
        log_print(f"✓ Found {total_genomes} genomes. Starting download...")

        max_genomes = total_genomes
        log_print(f"Will download: {max_genomes} genomes")

        downloaded_count = 0
        skipped_count = 0
        error_count = 0

        for i, assembly_id in enumerate(assembly_ids[:max_genomes], 1):
            try:
                log_print(f"\n[{i}/{max_genomes}] Processing assembly ID: {assembly_id}")

                handle = Entrez.esummary(db="assembly", id=assembly_id)
                record = Entrez.read(handle)
                handle.close()

                if record:
                    assembly_info = record["DocumentSummarySet"]["DocumentSummary"][0]
                    ftp_path = assembly_info.get("FtpPath_GenBank", "")
                    assembly_name = assembly_info.get("AssemblyAccession", "unknown")

                    if ftp_path:
                        file_name = ftp_path.split("/")[-1]
                        fna_file = f"{file_name}_genomic.fna.gz"
                        ftp_fna = f"{ftp_path}/{fna_file}"

                        output_file = taxon_dir / fna_file
                        if not output_file.exists():
                            log_print(f"  → Downloading: {fna_file}")
                            try:
                                cmd = ["wget", "--progress=bar:force", "-O", str(output_file), ftp_fna]
                                subprocess.run(cmd, check=True, capture_output=True, text=True)
                                file_size = output_file.stat().st_size / (1024*1024)
                                log_print(f"  ✓ {fna_file} downloaded ({file_size:.2f} MB)")
                                downloaded_count += 1
                            except subprocess.CalledProcessError as e:
                                log_print(f"  ✗ Error downloading {fna_file}: {e}")
                                error_count += 1
                                if output_file.exists():
                                    output_file.unlink()
                        else:
                            file_size = output_file.stat().st_size / (1024*1024)
                            log_print(f"  ⊗ {fna_file} already exists ({file_size:.2f} MB)")
                            skipped_count += 1

                        if i % 5 == 0:
                            log_print(f"\n📊 Progress: {i}/{max_genomes} | Downloaded: {downloaded_count} | Skipped: {skipped_count} | Errors: {error_count}")

                        time.sleep(0.34)
                    else:
                        log_print(f"  ⚠ FTP path not found for assembly {assembly_id}")
                        error_count += 1
                else:
                    log_print(f"  ⚠ Could not fetch assembly info for {assembly_id}")
                    error_count += 1

            except Exception as e:
                log_print(f"  ✗ Error processing assembly {assembly_id}: {e}")
                error_count += 1
                continue

        log_print(f"\n{'='*60}")
        log_print(f"✓ Finished {taxon}:")
        log_print(f"  - Processed: {max_genomes} genomes")
        log_print(f"  - Downloaded: {downloaded_count}")
        log_print(f"  - Skipped (existing): {skipped_count}")
        log_print(f"  - Errors: {error_count}")
        log_print(f"{'='*60}\n")

def main():
    taxa = [
        "Deinococcales",
    ]

    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOGS_DIR / "download_progress.log"
    log_file.unlink(missing_ok=True)

    def log_print(message):
        print(message)
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

    log_print("=" * 60)
    log_print("NCBI genome download")
    log_print("=" * 60)
    log_print(f"Output directory: {OUTPUT_DIR}")
    log_print(f"Log file: {log_file}")
    log_print(f"Taxa: {', '.join(taxa)}")
    log_print("=" * 60)

    start_time = time.time()

    if check_tool_available(NCBI_DATASETS_CMD):
        log_print("\n✓ Found tool: datasets")
        download_with_datasets(taxa, OUTPUT_DIR, assembly_level="complete")
    elif check_tool_available(NCBI_GENOME_DOWNLOAD_CMD):
        log_print("\n✓ Found tool: ncbi-genome-download")
        download_with_ncbi_genome_download(taxa, OUTPUT_DIR, assembly_level="complete")
    else:
        log_print("\n⚠ Specialized tools not found")
        log_print("Using Entrez + wget fallback...")

        if not check_tool_available("wget"):
            log_print("ERROR: wget not found. Install wget to download files.")
            sys.exit(1)

        download_genomes_via_entrez(taxa, OUTPUT_DIR, log_file=str(log_file))

    elapsed_time = time.time() - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)

    log_print("\n" + "=" * 60)
    log_print("✓ Download complete!")
    log_print(f"Elapsed time: {hours:02d}:{minutes:02d}:{seconds:02d}")
    log_print(f"Genomes saved in: {OUTPUT_DIR}")
    log_print(f"Log saved in: {log_file}")
    log_print("=" * 60)

if __name__ == "__main__":
    main()
