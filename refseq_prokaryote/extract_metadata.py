#!/usr/bin/env python3
"""
Extract metadata (ID and species) from RefSeq prokaryote genomes.
Creates a TSV table with file_id, species, and other metadata.
"""

import os
import gzip
import re
import glob
from pathlib import Path
from collections import defaultdict
import argparse


def parse_fasta_header(header):
    """
    Parse FASTA header to extract species information.
    Headers typically look like:
    >NZ_JBDYLZ010000005.1 Hafnia sp. KE9867 NODE_5_length_346550_cov_30.270611, whole genome shotgun sequence
    >NZ_CAYKQK010000032.1 MAG Limosilactobacillus vaginalis isolate ERR11768928_concoct_15 ERZ24910109.214, whole genome shotgun sequence
    """
    # Remove the sequence ID part (everything before first space)
    parts = header.split(' ', 1)
    if len(parts) < 2:
        return None, None, None
    
    seq_id = parts[0].lstrip('>')
    description = parts[1]
    
    # Try to extract species name
    # Pattern: Genus species or Genus sp.
    species_match = re.search(r'([A-Z][a-z]+(?:\s+[a-z]+)?(?:\s+sp\.)?)', description)
    if species_match:
        species = species_match.group(1)
    else:
        # Fallback: try to get first few words
        words = description.split()[:3]
        species = ' '.join(words) if words else 'Unknown'
    
    # Extract strain/isolate if present
    strain_match = re.search(r'(isolate|strain)\s+([A-Za-z0-9_\-]+)', description, re.IGNORECASE)
    strain = strain_match.group(2) if strain_match else None
    
    # Extract assembly type
    assembly_type = None
    if 'whole genome shotgun' in description.lower():
        assembly_type = 'WGS'
    elif 'complete genome' in description.lower():
        assembly_type = 'Complete'
    elif 'chromosome' in description.lower():
        assembly_type = 'Chromosome'
    
    return species, strain, assembly_type


def extract_metadata_from_file(fasta_path):
    """Extract metadata from a single FASTA file."""
    file_id = Path(fasta_path).stem.replace('.genomic.fna', '')
    
    species_list = []
    strains = []
    assembly_types = []
    seq_ids = []
    
    try:
        # Open gzipped or regular file
        if fasta_path.endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        with opener(fasta_path, mode) as f:
            for line in f:
                if line.startswith('>'):
                    species, strain, assembly_type = parse_fasta_header(line.strip())
                    if species:
                        species_list.append(species)
                    if strain:
                        strains.append(strain)
                    if assembly_type:
                        assembly_types.append(assembly_type)
                    # Extract sequence ID
                    seq_id = line.split()[0].lstrip('>')
                    seq_ids.append(seq_id)
                    
                    # Only read first few headers to get representative info
                    if len(species_list) >= 5:
                        break
    
    except Exception as e:
        print(f"Error reading {fasta_path}: {e}", file=os.sys.stderr)
        return None
    
    if not species_list:
        return None
    
    # Get most common species (or first if all unique)
    from collections import Counter
    species_counter = Counter(species_list)
    most_common_species = species_counter.most_common(1)[0][0]
    
    # Get unique values
    unique_strains = list(set(strains)) if strains else []
    unique_assembly_types = list(set(assembly_types)) if assembly_types else []
    
    return {
        'file_id': file_id,
        'filename': Path(fasta_path).name,
        'species': most_common_species,
        'strain': unique_strains[0] if unique_strains else None,
        'assembly_type': unique_assembly_types[0] if unique_assembly_types else None,
        'num_contigs': len(seq_ids),
        'path': fasta_path
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extract metadata from RefSeq prokaryote genomes'
    )
    parser.add_argument(
        '--input_dir',
        type=str,
        default='/home/lam34/MGE_finder/refseq_prokaryote',
        help='Base directory with bacteria/ and archaea/ subdirectories'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='/home/lam34/MGE_finder/refseq_prokaryote/metadata.tsv',
        help='Output TSV file path'
    )
    parser.add_argument(
        '--domain',
        type=str,
        choices=['bacteria', 'archaea', 'both'],
        default='both',
        help='Which domain to process'
    )
    
    args = parser.parse_args()
    
    # Find all FASTA files
    fasta_files = []
    
    if args.domain in ['bacteria', 'both']:
        bacteria_dir = os.path.join(args.input_dir, 'bacteria')
        if os.path.exists(bacteria_dir):
            fasta_files.extend(glob.glob(os.path.join(bacteria_dir, '*.fna.gz')))
            fasta_files.extend(glob.glob(os.path.join(bacteria_dir, '*.fna')))
    
    if args.domain in ['archaea', 'both']:
        archaea_dir = os.path.join(args.input_dir, 'archaea')
        if os.path.exists(archaea_dir):
            fasta_files.extend(glob.glob(os.path.join(archaea_dir, '*.fna.gz')))
            fasta_files.extend(glob.glob(os.path.join(archaea_dir, '*.fna')))
    
    print(f"Found {len(fasta_files)} FASTA files to process")
    
    # Extract metadata
    metadata_list = []
    for i, fasta_file in enumerate(sorted(fasta_files), 1):
        if i % 100 == 0:
            print(f"Processing {i}/{len(fasta_files)}...", file=os.sys.stderr)
        
        metadata = extract_metadata_from_file(fasta_file)
        if metadata:
            metadata_list.append(metadata)
    
    print(f"Extracted metadata from {len(metadata_list)} files")
    
    # Write TSV
    import csv
    with open(args.output, 'w', newline='') as f:
        if metadata_list:
            fieldnames = ['file_id', 'filename', 'species', 'strain', 'assembly_type', 'num_contigs', 'path']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for metadata in sorted(metadata_list, key=lambda x: x['file_id']):
                writer.writerow(metadata)
    
    print(f"Metadata written to {args.output}")
    print(f"Total entries: {len(metadata_list)}")


if __name__ == '__main__':
    main()

