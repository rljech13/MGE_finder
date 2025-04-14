import os
import re
import subprocess
import argparse
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from logger import Logger

logger = Logger(name="annotate_mge_region", level=Logger.Level.INFO).get_logger()

WINDOW_SIZE = 100000


def load_integrases(integrases_file):
    """Load integrase data from a TSV file.

    Reads a TSV file containing integrase information and returns a DataFrame.
    The file is expected to contain columns such as 'orf_id', 'start', 'end',
    'contig_length', and 'contig_id'.

    Args:
        integrases_file (str): Path to the integrases TSV file.

    Returns:
        pd.DataFrame: DataFrame containing integrase records.
    """
    df = pd.read_csv(integrases_file, sep="\t")
    return df


def get_integrase_record(df, integrase_id):
    """Retrieve an integrase record from the DataFrame by integrase_id.

    Searches the DataFrame for a row where the 'orf_id' column matches the given
    integrase_id. If no record is found, logs an error and returns None.

    Args:
        df (pd.DataFrame): DataFrame with integrase data.
        integrase_id (str): The integrase ID to search for.

    Returns:
        pd.Series or None: The integrase record as a Series, or None if not found.
    """
    rec = df[df["orf_id"] == integrase_id]
    if rec.empty:
        logger.error(f"Record for integrase {integrase_id} not found")
        return None
    return rec.iloc[0]


def extract_subject_region(genome_fasta, contig, region_start, region_end):
    """Extract a region from a genome FASTA file for a given contig and coordinates.

    Loads genome sequences from the FASTA file, finds the sequence corresponding to the given contig,
    and extracts the region specified by region_start and region_end. If the contig is not found, attempts
    to find a sequence whose key contains the contig string. The coordinates are adjusted to be within bounds.

    Args:
        genome_fasta (str): Path to the genome FASTA file.
        contig (str): The contig identifier.
        region_start (int): The start coordinate (1-indexed).
        region_end (int): The end coordinate (1-indexed).

    Returns:
        Bio.Seq.Seq or None: The extracted sequence region, or None if the contig is not found.
    """
    records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    if contig in records:
        rec = records[contig]
    else:
        rec = None
        for key in records:
            if contig in key:
                rec = records[key]
                logger.info(f"Using {key} for contig {contig}")
                break
    if rec is None:
        logger.error(f"Sequence for contig {contig} not found")
        return None
    # Adjust coordinates to be within bounds.
    region_start = max(region_start, 1)
    region_end = min(region_end, len(rec.seq))
    return rec.seq[region_start - 1:region_end]


def run_blast_on_region(query_seq, subject_seq, tmp_dir):
    """Run BLAST on the provided query sequence against a subject region.

    Creates temporary files for the query and subject sequences, builds a BLAST database
    from the subject sequence, and then runs blastn to compare the query against the subject.
    Returns the BLAST results as a DataFrame. If the output is empty, returns an empty DataFrame.

    Args:
        query_seq (Bio.SeqRecord.SeqRecord): The query sequence.
        subject_seq (Bio.Seq.Seq): The subject sequence from the genome.
        tmp_dir (str): Temporary directory to use for intermediate files.

    Returns:
        pd.DataFrame: DataFrame containing the BLAST results.
    """
    with tempfile.TemporaryDirectory(dir=tmp_dir) as temp_dir:
        query_file = os.path.join(temp_dir, "tmp_query.fa")
        subject_file = os.path.join(temp_dir, "tmp_subject.fa")
        blast_out = os.path.join(temp_dir, "tmp_blast.tsv")
        SeqIO.write([query_seq], query_file, "fasta")
        with open(subject_file, "w") as sf:
            SeqIO.write([SeqRecord(subject_seq, id="subject", description="")], sf, "fasta")
        try:
            proc = subprocess.run(
                ["makeblastdb", "-in", subject_file, "-dbtype", "nucl"],
                check=True, capture_output=True, text=True
            )
            logger.info(f"makeblastdb output: {proc.stdout}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating BLAST database: {e.stderr}")
            raise
        try:
            proc = subprocess.run(
                [
                    "blastn",
                    "-query", query_file,
                    "-db", subject_file,
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                    "-out", blast_out,
                    "-word_size", "4",         # <-- минимальное значение
                    "-dust", "no"              # <-- отключить low complexity filter
                ],
                check=True, capture_output=True, text=True
            )
            logger.info(f"blastn output: {proc.stdout}")
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error: {e.stderr}")
            raise
        try:
            df = pd.read_csv(blast_out, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()
        return df


def main(genome_fasta, integrases_file, query_fa, out_tsv, tmp_dir):
    """Annotate BLAST hits within a search window around integrases.

    Loads integrase data from a TSV file and query sequences from a FASTA file.
    For each query, it extracts the integrase record, calculates the search window
    based on tRNA orientation, extracts the subject region from the genome FASTA file,
    runs BLAST on that region, and collects the BLAST hit data. Finally, writes all
    results to an output TSV file.

    Args:
        genome_fasta (str): Path to the full genome FASTA file.
        integrases_file (str): Path to the integrases TSV file (e.g., integrase_hits_summary.tsv).
        query_fa (str): Path to the FASTA file with query sequences (e.g., mge_query.fa).
        out_tsv (str): Path to the output TSV file for BLAST hits.
        tmp_dir (str): Temporary directory to use for intermediate files.

    Returns:
        None
    """
    integrases_df = load_integrases(integrases_file)
    queries = list(SeqIO.parse(query_fa, "fasta"))
    results = []
    for query in queries:
        try:
            # Expected header format: integrase_id:contig:trna_start-trna_end:trna_strand
            integrase_id, contig, trna_range, trna_strand = query.id.split(":")
            trna_start_str, trna_end_str = trna_range.split("-")
            trna_start = int(trna_start_str)
            trna_end = int(trna_end_str)
        except Exception as e:
            logger.error(f"Error parsing header {query.id}: {e}")
            continue
        rec = get_integrase_record(integrases_df, integrase_id)
        if rec is None:
            continue

        # Determine search window based on tRNA orientation.
        if trna_strand.strip() == "+":
            # For tRNA on '+' strand, window = [trna_end, trna_end + WINDOW_SIZE]
            region_start = trna_end
            region_end = trna_end + WINDOW_SIZE
        else:
            # For tRNA on '-' strand, window = [max(trna_start - WINDOW_SIZE, 1), trna_start]
            region_start = max(trna_start - WINDOW_SIZE, 1)
            region_end = trna_start

        logger.info(f"For integrase {integrase_id} on contig {contig} ({trna_strand}): window = [{region_start}, {region_end}]")
        subj_seq = extract_subject_region(genome_fasta, contig, region_start, region_end)
        if subj_seq is None or len(subj_seq) == 0:
            logger.error(f"Failed to extract subject for {integrase_id} on contig {contig}")
            continue

        # Query sequence is used as is (it is assumed that mge_query.fa is correctly prepared)
        q_seq = query.seq
        try:
            df_blast = run_blast_on_region(SeqIO.SeqRecord(q_seq, id=query.id, description=""), subj_seq, tmp_dir)
            if not df_blast.empty:
                logger.info(f"Found {len(df_blast)} BLAST hits for {query.id}")
                for i, row in df_blast.iterrows():
                    full_start = region_start + int(row[8]) - 1
                    full_end = region_start + int(row[9]) - 1
                    results.append({
                        "integrase_id": integrase_id,
                        "contig": contig,
                        "hit_start": full_start,
                        "hit_end": full_end,
                        "pident": row[2],
                        "length": row[3],
                        "evalue": row[10],
                        "bitscore": row[11]
                    })
            else:
                logger.info(f"No BLAST hits found for {query.id}")
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error for {query.id}: {e}")
            continue

    if results:
        df_res = pd.DataFrame(results)
        df_res.to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"Saved {len(results)} BLAST hits to {out_tsv}")
    else:
        open(out_tsv, "w").close()
        logger.info(f"No BLAST hits; created empty file {out_tsv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate BLAST hits in the window around an integrase")
    parser.add_argument("--ffn", required=True, help="Path to the full genome FASTA file")
    parser.add_argument("--integrases", required=True, help="TSV file with integrase data (e.g., integrase_hits_summary.tsv)")
    parser.add_argument("--query", required=True, help="FASTA file with query sequences (e.g., mge_query.fa)")
    parser.add_argument("--out_tsv", required=True, help="Output TSV file for BLAST hits")
    parser.add_argument("--tmp_dir", default=".", help="Temporary directory for intermediate files")
    args = parser.parse_args()

    main(args.ffn, args.integrases, args.query, args.out_tsv, args.tmp_dir)