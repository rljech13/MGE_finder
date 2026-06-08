#!/usr/bin/env python3
"""Build taxonomy summary tables from Deinococcales metadata."""

import csv
from collections import defaultdict
from pathlib import Path

DATA_ROOT = Path(__file__).resolve().parent.parent
METADATA_DIR = DATA_ROOT / "metadata"
METADATA_FILE = METADATA_DIR / "deinococcales_metadata_complete.tsv"
if not METADATA_FILE.exists():
    METADATA_FILE = METADATA_DIR / "deinococcales_metadata.tsv"
SUMMARY_FILE = METADATA_DIR / "deinococcales_taxonomy_summary.tsv"
DETAILED_SUMMARY = METADATA_DIR / "deinococcales_taxonomy_detailed.md"

def parse_organism_name(organism):
    """Parse organism name to extract genus and species."""
    if not organism or organism == "Unknown":
        return "Unknown", "Unknown"

    organism = organism.replace(" (bacteria)", "").strip()

    if "uncultured" in organism.lower():
        parts = organism.split()
        if len(parts) >= 2:
            genus_part = parts[1] if parts[0].lower() == "uncultured" else parts[0]
            if genus_part.startswith("f__"):
                genus = genus_part[3:]
            elif genus_part.endswith("aceae"):
                genus = genus_part
            elif genus_part.endswith("ales"):
                genus = genus_part
            else:
                genus = genus_part
            species = " ".join(parts[2:]) if len(parts) > 2 else "bacterium"
            return genus, species

    parts = organism.split()
    if len(parts) >= 2:
        genus = parts[0]
        if genus.startswith("f__"):
            genus = genus[3:]
        species = " ".join(parts[1:3]) if len(parts) >= 3 else parts[1]
        return genus, species
    if len(parts) == 1:
        return parts[0], "sp."
    return "Unknown", organism

def main():
    print("Reading metadata...")

    all_data = []
    with open(METADATA_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            all_data.append(row)

    print(f"Loaded records: {len(all_data)}")

    genera = defaultdict(lambda: {"count": 0, "species": defaultdict(int), "total_size_mb": 0})
    species_full = defaultdict(lambda: {"count": 0, "total_size_mb": 0, "genus": ""})

    for row in all_data:
        organism = row.get("Organism", "Unknown")
        genus_field = row.get("Genus", "")
        species_field = row.get("Species", "")

        if genus_field and genus_field != "Unknown":
            genus = genus_field
        else:
            genus, _ = parse_organism_name(organism)

        if species_field and species_field != "Unknown":
            species = species_field
        else:
            _, species = parse_organism_name(organism)

        if not species or species == "Unknown":
            species = "sp."

        full_species = f"{genus} {species}".strip()

        try:
            size_mb = float(row.get("File_Size_MB", 0))
        except ValueError:
            size_mb = 0

        genera[genus]["count"] += 1
        genera[genus]["species"][full_species] += 1
        genera[genus]["total_size_mb"] += size_mb

        species_full[full_species]["count"] += 1
        species_full[full_species]["total_size_mb"] += size_mb
        species_full[full_species]["genus"] = genus

    summary_data = []
    for genus, data in sorted(genera.items(), key=lambda x: x[1]["count"], reverse=True):
        summary_data.append({
            "Genus": genus,
            "Number_of_Genomes": data["count"],
            "Number_of_Species": len(data["species"]),
            "Total_Size_MB": round(data["total_size_mb"], 2),
            "Avg_Size_MB": round(data["total_size_mb"] / data["count"], 2) if data["count"] > 0 else 0,
            "Species_List": "; ".join(sorted(data["species"].keys()))
        })

    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(SUMMARY_FILE, 'w', newline='', encoding='utf-8') as f:
        fieldnames = ["Genus", "Number_of_Genomes", "Number_of_Species", "Total_Size_MB", "Avg_Size_MB", "Species_List"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_data)

    print(f"✓ Summary table saved: {SUMMARY_FILE}")

    with open(DETAILED_SUMMARY, 'w', encoding='utf-8') as f:
        f.write("# Deinococcales genome taxonomy analysis\n\n")
        f.write(f"**Total genomes with metadata:** {len(all_data)}\n\n")

        f.write("## Statistics by genus\n\n")
        f.write("| Genus | Number of genomes | Number of species | Total size (MB) | Average size (MB) |\n")
        f.write("|-------|-------------------|-------------------|-----------------|-------------------|\n")

        for row in summary_data:
            f.write(f"| {row['Genus']} | {row['Number_of_Genomes']} | {row['Number_of_Species']} | "
                   f"{row['Total_Size_MB']} | {row['Avg_Size_MB']} |\n")

        f.write("\n## Species breakdown\n\n")
        f.write("| Species | Number of genomes | Total size (MB) |\n")
        f.write("|---------|-------------------|-----------------|\n")

        sorted_species = sorted(species_full.items(), key=lambda x: x[1]["count"], reverse=True)
        for species, data in sorted_species:
            f.write(f"| {species} | {data['count']} | {round(data['total_size_mb'], 2)} |\n")

        f.write("\n## All species by genus\n\n")
        for genus, data in sorted(genera.items(), key=lambda x: x[1]["count"], reverse=True):
            f.write(f"### {genus} ({data['count']} genomes)\n\n")
            sorted_sp = sorted(data["species"].items(), key=lambda x: x[1], reverse=True)
            for species, count in sorted_sp:
                f.write(f"- **{species}**: {count} genomes\n")
            f.write("\n")

    print(f"✓ Detailed report saved: {DETAILED_SUMMARY}")

    print("\n" + "="*60)
    print("SUMMARY:")
    print("="*60)
    print(f"Total genera: {len(genera)}")
    print(f"Total species: {len(species_full)}")
    print("\nTop 10 genera:")
    for i, row in enumerate(summary_data[:10], 1):
        print(f"  {i}. {row['Genus']}: {row['Number_of_Genomes']} genomes, {row['Number_of_Species']} species")

    print("\nTop 10 species:")
    sorted_species = sorted(species_full.items(), key=lambda x: x[1]["count"], reverse=True)
    for i, (species, data) in enumerate(sorted_species[:10], 1):
        print(f"  {i}. {species}: {data['count']} genomes")

if __name__ == "__main__":
    main()
