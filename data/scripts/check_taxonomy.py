#!/usr/bin/env python3
"""
Check taxonomy of Deinococcales genomes before bulk download.
Shows which genera and species will be included.
"""

import os
import sys
from pathlib import Path
from Bio import Entrez
import time

Entrez.email = os.environ.get("NCBI_EMAIL", "genome.download@example.com")
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "ncbi"

def get_taxonomy_from_assembly(assembly_id):
    """Return taxonomy fields from an assembly summary."""
    try:
        handle = Entrez.esummary(db="assembly", id=assembly_id)
        record = Entrez.read(handle)
        handle.close()

        if record and "DocumentSummarySet" in record:
            doc = record["DocumentSummarySet"]["DocumentSummary"][0]
            return {
                "genus": doc.get("Genus", "Unknown"),
                "species": doc.get("SpeciesName", "Unknown"),
                "organism": doc.get("Organism", "Unknown"),
                "taxid": doc.get("Taxid", "Unknown"),
            }
    except Exception as e:
        print(f"Error fetching taxonomy for {assembly_id}: {e}")

    return None

def get_all_assembly_ids_for_taxon(taxon_name):
    """Return all assembly IDs for a taxon."""
    print(f"Searching all genomes for {taxon_name}...")

    search_term = f'"{taxon_name}"[Organism] AND "latest"[filter]'

    try:
        handle = Entrez.esearch(db="assembly", term=search_term, retmax=10000)
        record = Entrez.read(handle)
        handle.close()

        count = int(record["Count"])
        print(f"Found {count} genomes")

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
    except Exception as e:
        print(f"Error: {e}")

    return []

def analyze_taxonomy(assembly_ids, sample_size=50):
    """Analyze taxonomy for a sample of assemblies."""
    print(f"\nAnalyzing taxonomy for {min(sample_size, len(assembly_ids))} genomes...")

    genera = {}
    species_list = []
    sample_ids = assembly_ids[:sample_size]

    for i, assembly_id in enumerate(sample_ids, 1):
        print(f"[{i}/{len(sample_ids)}] Processing assembly {assembly_id}...", end="\r")

        tax_info = get_taxonomy_from_assembly(assembly_id)

        if tax_info:
            genus = tax_info["genus"]
            species = tax_info["species"]

            if genus not in genera:
                genera[genus] = []
            genera[genus].append(species)
            species_list.append(f"{genus} {species}")

        time.sleep(0.34)

    print("\n" + "="*60)
    print("TAXONOMY ANALYSIS RESULTS")
    print("="*60)

    print(f"\nAnalyzed: {len(sample_ids)} genomes")
    print(f"Unique genera: {len(genera)}")

    print("\n" + "="*60)
    print("GENERA AND SPECIES (by frequency):")
    print("="*60)

    sorted_genera = sorted(genera.items(), key=lambda x: len(x[1]), reverse=True)

    for genus, species in sorted_genera:
        unique_species = list(set(species))
        print(f"\n{genus} ({len(unique_species)} species, {len(species)} genomes):")
        for sp in sorted(unique_species):
            count = species.count(sp)
            print(f"  - {sp} ({count} genomes)")

    print("\n" + "="*60)
    print("ALL UNIQUE SPECIES:")
    print("="*60)
    for sp in sorted(set(species_list)):
        print(f"  - {sp}")

    return genera

def main():
    print("="*60)
    print("DEINOCOCCALES TAXONOMY ANALYSIS")
    print("="*60)

    assembly_ids = get_all_assembly_ids_for_taxon("Deinococcales")

    if not assembly_ids:
        print("No genomes found!")
        return

    print(f"\nTotal found: {len(assembly_ids)} genomes")

    sample_size = min(100, len(assembly_ids))
    analyze_taxonomy(assembly_ids, sample_size)

    print("\n" + "="*60)
    print(f"Note: analyzed a sample of {sample_size} genomes")
    print(f"Total available for download: {len(assembly_ids)} genomes")
    print("="*60)

if __name__ == "__main__":
    main()
