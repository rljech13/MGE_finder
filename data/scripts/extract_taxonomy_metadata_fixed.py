#!/usr/bin/env python3
"""
Extract taxonomy metadata for downloaded Deinococcales genomes.
Uses per-GCA Entrez lookups for reliable assembly matching.
"""

import os
import sys
import re
import time
from pathlib import Path
from collections import defaultdict
import csv

try:
    from Bio import Entrez
except ImportError:
    print("Install biopython: conda install -c bioconda biopython")
    sys.exit(1)

Entrez.email = os.environ.get("NCBI_EMAIL", "genome.download@example.com")

DATA_ROOT = Path(__file__).resolve().parent.parent
NCBI_DIR = DATA_ROOT / "ncbi"
METADATA_DIR = DATA_ROOT / "metadata"
DEINOCOCCALES_DIR = NCBI_DIR / "Deinococcales"
OUTPUT_TSV = METADATA_DIR / "deinococcales_metadata_complete.tsv"

def extract_gca_from_filename(filename):
    """Extract GCA accession from a filename."""
    match = re.search(r'GCA_\d+\.\d+', filename)
    if match:
        return match.group(0)
    return None

def get_metadata_for_gca(gca_accession):
    """Fetch metadata for one GCA accession."""
    try:
        search_term = f'"{gca_accession}"[Assembly Accession]'
        handle = Entrez.esearch(db="assembly", term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return None

        assembly_id = record["IdList"][0]

        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(handle)
        handle.close()

        if summary and "DocumentSummarySet" in summary:
            doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
            return {
                "GCA_Accession": gca_accession,
                "Assembly_Accession": gca_accession,
                "Assembly_Name": doc.get("AssemblyName", "Unknown"),
                "Organism": doc.get("Organism", "Unknown"),
                "Genus": doc.get("Genus", "Unknown"),
                "Species": doc.get("SpeciesName", "Unknown"),
                "TaxID": doc.get("Taxid", "Unknown"),
                "Assembly_Level": doc.get("AssemblyStatus", "Unknown"),
                "Total_Length": doc.get("TotalLength", 0),
                "Contig_N50": doc.get("ContigN50", 0),
                "Scaffold_N50": doc.get("ScaffoldN50", 0),
                "BioProject": doc.get("BioProjectAccn", "Unknown"),
                "BioSample": doc.get("BioSampleAccn", "Unknown"),
                "Submission_Date": doc.get("SubmissionDate", "Unknown"),
                "Release_Date": doc.get("ReleaseDate", "Unknown"),
            }
    except Exception:
        return None

    return None

def get_file_size_mb(filepath):
    """Return file size in MB."""
    try:
        return filepath.stat().st_size / (1024 * 1024)
    except OSError:
        return 0

def main():
    print("="*60)
    print("EXTRACTING TAXONOMY METADATA")
    print("="*60)

    if not DEINOCOCCALES_DIR.exists():
        print(f"Error: directory not found: {DEINOCOCCALES_DIR}")
        sys.exit(1)

    genome_files = list(DEINOCOCCALES_DIR.glob("*.fna.gz"))
    print(f"\nFound files: {len(genome_files)}")

    gca_files = {}
    for filepath in genome_files:
        gca = extract_gca_from_filename(filepath.name)
        if gca:
            gca_files[gca] = filepath

    print(f"Extracted GCA accessions: {len(gca_files)}")

    existing_gcas = set()
    old_metadata_file = METADATA_DIR / "deinococcales_metadata.tsv"
    if old_metadata_file.exists():
        with open(old_metadata_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                existing_gcas.add(row.get("GCA_Accession", ""))
        print(f"Already in metadata: {len(existing_gcas)}")

    all_metadata = []
    failed = []

    gca_list = list(gca_files.keys())
    total = len(gca_list)

    print(f"\nFetching metadata from NCBI...")
    print("Estimated time: {:.1f} minutes...".format(total * 0.34 / 60))

    for i, gca in enumerate(gca_list, 1):
        if gca in existing_gcas:
            continue

        print(f"[{i}/{total}] Processing {gca}...", end="\r")

        filepath = gca_files[gca]
        metadata = get_metadata_for_gca(gca)

        if metadata:
            metadata["File_Size_MB"] = round(get_file_size_mb(filepath), 2)
            metadata["Filename"] = filepath.name
            all_metadata.append(metadata)
        else:
            failed.append(gca)

        time.sleep(0.34)

    if old_metadata_file.exists():
        with open(old_metadata_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gca = row.get("GCA_Accession", "")
                if gca in gca_files:
                    filepath = gca_files[gca]
                    row["File_Size_MB"] = round(get_file_size_mb(filepath), 2)
                    row["Filename"] = filepath.name
                all_metadata.append(row)

    print(f"\n\nSuccessfully processed: {len(all_metadata)}")
    print(f"Failures: {len(failed)}")

    if failed and len(failed) <= 20:
        print("\nCould not fetch metadata for:")
        for gca in failed:
            print(f"  - {gca}")
    elif failed:
        print(f"\nCould not fetch metadata for {len(failed)} files")

    if all_metadata:
        fieldnames = [
            "GCA_Accession", "Assembly_Accession", "Assembly_Name",
            "Organism", "Genus", "Species", "TaxID",
            "Assembly_Level", "Total_Length", "Contig_N50", "Scaffold_N50",
            "BioProject", "BioSample", "Submission_Date", "Release_Date",
            "File_Size_MB", "Filename"
        ]

        METADATA_DIR.mkdir(parents=True, exist_ok=True)
        with open(OUTPUT_TSV, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(all_metadata)

        print(f"\n✓ Metadata saved to: {OUTPUT_TSV}")
        print(f"  Total records: {len(all_metadata)}")

        genera = defaultdict(int)
        species_dict = defaultdict(int)

        for meta in all_metadata:
            genus = meta.get("Genus", "Unknown")
            species = meta.get("Species", "Unknown")
            organism = meta.get("Organism", "Unknown")

            if genus == "Unknown" and organism != "Unknown":
                parts = organism.replace(" (bacteria)", "").split()
                if parts:
                    genus = parts[0]

            if genus != "Unknown":
                genera[genus] += 1

            if species != "Unknown":
                full_species = f"{genus} {species}" if genus != "Unknown" else species
                species_dict[full_species] += 1

        print(f"\nTotal genera: {len(genera)}")
        print(f"Total species: {len(species_dict)}")

        if genera:
            print("\nTop 15 genera:")
            sorted_genera = sorted(genera.items(), key=lambda x: x[1], reverse=True)
            for i, (genus, count) in enumerate(sorted_genera[:15], 1):
                print(f"  {i}. {genus}: {count} genomes")

    print("\n" + "="*60)
    print("DONE!")
    print("="*60)

if __name__ == "__main__":
    main()
