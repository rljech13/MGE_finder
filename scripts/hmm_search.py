import os
import re
import subprocess
import argparse
from logger import Logger

logger = Logger(name="hmm_search").get_logger()


def run_hmmscan(faa_path, hmm_path, output_tbl):
    """Execute hmmscan on the given FAA file.

    This function runs the hmmscan command with the specified combined HMM file and
    writes the results to the output table file.

    Args:
        faa_path (str): Path to the FAA file (protein sequences).
        hmm_path (str): Path to the combined HMM file.
        output_tbl (str): Path to the output table file (tblout format).
    """
    logger.info(f"Executing hmmscan for {faa_path} using {hmm_path}...")
    subprocess.run([
        "hmmscan",
        "--tblout", output_tbl,
        hmm_path,
        faa_path
    ], check=True)
    logger.info(f"Results saved in {output_tbl}")


def parse_tblout(tbl_path):
    """Parse the tblout file and create a mapping from ORF ID to model accession.

    The tblout file (e.g., integrase_hits.txt) is expected to have lines in the format:
      [Some header text]
      Phage_integrase      PF00589.27 JBKBIM010000027.1_8  - 4.5e-35 108.1 ...
    where the query name (ORF ID) is the third column and the model accession is the second column.

    Args:
        tbl_path (str): Path to the tblout file.

    Returns:
        dict: A dictionary mapping ORF IDs (str) to model accessions (str).
    """
    mapping = {}
    with open(tbl_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            orf_id = parts[2]
            model_acc = parts[1]
            mapping[orf_id] = model_acc
    return mapping


def parse_faa(faa_path):
    """Parse the FAA file to extract ORF coordinates.

    The FAA file header is expected to have the format:
      >JBKBIM010000001.1_1 # 3 # 1430 # -1 # ID=1_1;partial=10;...
    The function returns a dictionary mapping the ORF ID to a tuple of (start, end, strand).

    Args:
        faa_path (str): Path to the FAA file.

    Returns:
        dict: A dictionary mapping ORF IDs (str) to tuples (start (int), end (int), strand (str)).
    """
    orf_coords = {}
    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split(' # ')
                if len(parts) < 4:
                    logger.error(f"Incorrect FAA header format: {line.strip()}")
                    continue
                orf_id = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    logger.error(f"Error converting coordinates in FAA: {line.strip()}")
                    continue
                strand = parts[3]
                orf_coords[orf_id] = (start, end, strand)
    return orf_coords


def parse_gff(gff_path):
    """Parse the GFF file to extract contig lengths.

    The GFF file is expected to contain a line with the following format:
      # Sequence Data: seqnum=1;seqlen=128375;seqhdr="JBKBIM010000027.1 MAG: Thermus sp. isolate ..."
    This function extracts the sequence length (seqlen) and the contig ID (the first token in seqhdr).

    Args:
        gff_path (str): Path to the GFF file.

    Returns:
        dict: A dictionary mapping contig IDs (str) to their lengths (int).
    """
    contig_lengths = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                m = re.search(r'seqlen=(\d+);seqhdr="([^"]+)"', line)
                if m:
                    seqlen = int(m.group(1))
                    seqhdr = m.group(2)
                    contig_id = seqhdr.split()[0]
                    contig_lengths[contig_id] = seqlen
                else:
                    logger.error(f"Incorrect GFF line: {line.strip()}")
    return contig_lengths


def write_outputs(mapping, faa_coords, contig_lengths, summary_file, orfs_file):
    """Write the integrase summary and ORF coordinate files.

    The summary file contains the following columns:
      ORF ID, model accession, start, end, strand, contig ID, contig length.
    The ORFs file contains the following columns:
      ORF ID, start, end.

    Args:
        mapping (dict): Mapping from ORF ID to model accession.
        faa_coords (dict): Mapping from ORF ID to (start, end, strand).
        contig_lengths (dict): Mapping from contig ID to contig length.
        summary_file (str): Path to the output summary file.
        orfs_file (str): Path to the output ORFs file.
    """
    with open(summary_file, 'w') as summ_f, open(orfs_file, 'w') as orfs_f:
        # Write header for summary file
        summ_f.write("orf_id\tmodel_accession\tstart\tend\tstrand\tcontig_id\tcontig_length\n")
        # Write header for ORFs file
        orfs_f.write("orf_id\tstart\tend\n")
        for orf_id, model in mapping.items():
            if orf_id not in faa_coords:
                logger.error(f"ORF {orf_id} not found in FAA")
                continue
            start, end, strand = faa_coords[orf_id]
            contig_id = orf_id.split('_')[0]
            if contig_id not in contig_lengths:
                logger.error(f"Contig {contig_id} not found in GFF")
                continue
            contig_len = contig_lengths[contig_id]
            summ_f.write(f"{orf_id}\t{model}\t{start}\t{end}\t{strand}\t{contig_id}\t{contig_len}\n")
            orfs_f.write(f"{orf_id}\t{start}\t{end}\n")
    logger.info("Data successfully written to output files.")


def main():
    parser = argparse.ArgumentParser(
        description="Match HMM scan results with ORF and contig data"
    )
    parser.add_argument("--faa", required=True, help="Path to FAA file (orfs.faa)")
    parser.add_argument("--gff", required=True, help="Path to GFF file (orfs.gff)")
    parser.add_argument("--out", required=True, help="Path to tblout file (integrase_hits.txt)")
    parser.add_argument("--summary", required=True, help="Output file for integrase summary table")
    parser.add_argument("--orfs", required=True, help="Output file for ORF coordinates")
    parser.add_argument("--combined", required=True, help="Path to combined HMM file")
    args = parser.parse_args()

    # Step 1: Run hmmscan
    run_hmmscan(args.faa, args.combined, args.out)
    # Step 2: Parse tblout file to get mapping (ORF ID -> model accession)
    mapping = parse_tblout(args.out)
    # Step 3: Parse FAA file to get ORF coordinates
    faa_coords = parse_faa(args.faa)
    # Step 4: Parse GFF file to get contig lengths
    contig_lengths = parse_gff(args.gff)
    # Step 5: Write output summary and ORF coordinate files
    write_outputs(mapping, faa_coords, contig_lengths, args.summary, args.orfs)


if __name__ == "__main__":
    main()