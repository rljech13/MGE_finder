#!/usr/bin/env python3
"""Diagnose genomes that lack metadata entries."""

import re
from pathlib import Path
from Bio import Entrez
import time

DATA_ROOT = Path(__file__).resolve().parent.parent
NCBI_DIR = DATA_ROOT / "ncbi"
METADATA_DIR = DATA_ROOT / "metadata"
DEINOCOCCALES_DIR = NCBI_DIR / "Deinococcales"
METADATA_FILE = METADATA_DIR / "deinococcales_metadata_complete.tsv"
if not METADATA_FILE.exists():
    METADATA_FILE = METADATA_DIR / "deinococcales_metadata.tsv"

Entrez.email = "genome.download@example.com"

def extract_gca_from_filename(filename):
    """Extract GCA accession from a filename."""
    match = re.search(r'GCA_\d+\.\d+', filename)
    if match:
        return match.group(0)
    return None

def get_assembly_info(gca_accession):
    """Look up assembly information for a GCA accession."""
    try:
        search_term = f'"{gca_accession}"[Assembly Accession]'
        handle = Entrez.esearch(db="assembly", term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            assembly_id = record["IdList"][0]

            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()

            if summary and "DocumentSummarySet" in summary:
                doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
                return {
                    "found": True,
                    "organism": doc.get("Organism", "Unknown"),
                    "taxid": doc.get("Taxid", "Unknown"),
                    "assembly_id": assembly_id,
                }

        return {"found": False, "reason": "Assembly ID not found"}
    except Exception as e:
        return {"found": False, "reason": f"Error: {str(e)}"}

def main():
    print("="*60)
    print("CHECKING GENOMES WITHOUT METADATA")
    print("="*60)

    genomes_with_metadata = set()
    if METADATA_FILE.exists():
        with open(METADATA_FILE, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("GCA_"):
                    parts = line.split('\t')
                    if parts:
                        genomes_with_metadata.add(parts[0])

    print(f"\nGenomes with metadata: {len(genomes_with_metadata)}")

    all_files = list(DEINOCOCCALES_DIR.glob("*.fna.gz"))
    print(f"Total files: {len(all_files)}")

    missing_metadata = []
    for filepath in all_files:
        gca = extract_gca_from_filename(filepath.name)
        if gca and gca not in genomes_with_metadata:
            missing_metadata.append((gca, filepath))

    print(f"Genomes without metadata: {len(missing_metadata)}")

    if not missing_metadata:
        print("\nAll genomes have metadata!")
        return

    print(f"\nChecking first 20 genomes without metadata...")
    print("This may take a while...\n")

    results = []
    for i, (gca, filepath) in enumerate(missing_metadata[:20], 1):
        print(f"[{i}/20] Checking {gca}...", end="\r")
        result = get_assembly_info(gca)
        result["gca"] = gca
        result["filename"] = filepath.name
        results.append(result)
        time.sleep(0.34)

    print("\n" + "="*60)
    print("CHECK RESULTS:")
    print("="*60)

    found_count = sum(1 for r in results if r["found"])
    not_found_count = len(results) - found_count

    print(f"\nFound in NCBI: {found_count}/{len(results)}")
    print(f"Not found: {not_found_count}/{len(results)}")

    if found_count > 0:
        print("\n✓ Genomes FOUND in NCBI (but missing from metadata):")
        for r in results:
            if r["found"]:
                print(f"  - {r['gca']}: {r['organism']} (TaxID: {r['taxid']})")

    if not_found_count > 0:
        print("\n✗ Genomes NOT FOUND in NCBI:")
        for r in results:
            if not r["found"]:
                print(f"  - {r['gca']}: {r.get('reason', 'Unknown reason')}")

    print("\n" + "="*60)
    print("FILENAME PATTERN ANALYSIS (missing metadata):")
    print("="*60)

    gca_patterns = {}
    for gca, filepath in missing_metadata:
        pattern = re.match(r'GCA_(\d+)', gca)
        if pattern:
            first_digits = pattern.group(1)[:3]
            gca_patterns.setdefault(first_digits, []).append(gca)

    print(f"\nDistribution by leading GCA digits:")
    for pattern, gcas in sorted(gca_patterns.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        print(f"  GCA_{pattern}...: {len(gcas)} genomes")

    print("\n" + "="*60)
    print("SUMMARY:")
    print("="*60)
    print(f"\n1. Total genomes without metadata: {len(missing_metadata)}")
    print(f"2. Checked: {len(results)}")
    print(f"3. Found in NCBI: {found_count}")
    print(f"4. Not found in NCBI: {not_found_count}")

    if found_count > 0:
        print("\n⚠ WARNING: Some genomes exist in NCBI but are missing from metadata!")
        print("   Possible causes:")
        print("   - batch request issues")
        print("   - incorrect GCA to assembly ID mapping")
        print("   - errors while processing results")

    if not_found_count > 0:
        print("\n⚠ Some genomes were not found in NCBI.")
        print("   Possible causes:")
        print("   - outdated or removed GCA accession")
        print("   - incorrect accession format")
        print("   - genome not yet published in Assembly database")

if __name__ == "__main__":
    main()
