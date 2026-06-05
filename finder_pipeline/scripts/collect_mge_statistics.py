"""Aggregate per-sample and per-element statistics from finder pipeline outputs."""

import argparse
import glob
import os
import re

import pandas as pd
from Bio import SeqIO

MIN_MGE_SIZE_BP = 200


def collect_file(results_dir, sample, filename):
    """Return the path to a sample output file if it exists.

    Args:
        results_dir: Root directory containing one subdirectory per sample.
        sample: Sample identifier (subdirectory name).
        filename: Output file name within the sample directory.

    Returns:
        Absolute path to the file, or None if the file is missing.
    """
    path = os.path.join(results_dir, sample, filename)
    return path if os.path.exists(path) else None


def get_sample_list(results_dir):
    """List sample identifiers present under a results directory.

    Args:
        results_dir: Root directory containing one subdirectory per sample.

    Returns:
        Sorted list of sample names (immediate subdirectory basenames).
    """
    samples = [
        os.path.basename(d)
        for d in glob.glob(os.path.join(results_dir, "*"))
        if os.path.isdir(d)
    ]
    return sorted(samples)


def parse_summary_file(file_path):
    """Load a tab-separated summary table.

    Args:
        file_path: Path to a TSV file.

    Returns:
        Parsed DataFrame, or an empty DataFrame when reading fails.
    """
    try:
        return pd.read_csv(file_path, sep="\t")
    except Exception:
        return pd.DataFrame()


def parse_fasta_lengths(fasta_path):
    """Extract sequence lengths from a FASTA file.

    Args:
        fasta_path: Path to a FASTA file (for example ``mge_region.fa``).

    Returns:
        Dictionary mapping sequence identifier to length in base pairs.
    """
    lengths = {}
    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            lengths[rec.id] = len(rec.seq)
    except Exception:
        pass
    return lengths


def parse_blast_lengths(blast_path):
    """Extract BLAST hit lengths grouped by query identifier.

    Supports both legacy columns (``qseqid``, ``qstart``, ``qend``) and the
    current format (``integrase_id``, ``qstart``, ``qend``).

    Args:
        blast_path: Path to a tab-separated BLAST results file.

    Returns:
        Dictionary mapping query identifier to a list of hit lengths in base pairs.
    """
    lengths = {}
    try:
        df = pd.read_csv(blast_path, sep="\t")

        if {"qstart", "qend"}.issubset(df.columns):
            start_col, end_col = "qstart", "qend"
        elif {"hit_start", "hit_end"}.issubset(df.columns):
            start_col, end_col = "hit_start", "hit_end"
        else:
            return lengths

        if "qseqid" in df.columns:
            id_col = "qseqid"
        elif "integrase_id" in df.columns:
            id_col = "integrase_id"
        else:
            return lengths

        for _, row in df.iterrows():
            query_id = str(row[id_col])
            hit_len = abs(int(row[end_col]) - int(row[start_col])) + 1
            lengths.setdefault(query_id, []).append(hit_len)
    except Exception:
        pass
    return lengths


def parse_gff(gff_path):
    """Parse Prodigal GFF contig metadata and compute total assembly size.

    Reads ``# Sequence Data:`` header lines emitted by Prodigal, for example::

        # Sequence Data: seqnum=1;seqlen=128375;seqhdr="CONTIG_ID ..."

    Args:
        gff_path: Path to a Prodigal GFF file.

    Returns:
        Tuple ``(contig_lengths, total_genome_size)`` where ``contig_lengths``
        maps contig identifier to length and ``total_genome_size`` is the sum
        of all contig lengths.
    """
    contig_lengths = {}
    total_size = 0
    try:
        with open(gff_path) as f:
            for line in f:
                if line.startswith("# Sequence Data:"):
                    match = re.search(r'seqlen=(\d+);seqhdr="([^"]+)"', line)
                    if match:
                        seqlen = int(match.group(1))
                        contig_id = match.group(2).split()[0]
                        contig_lengths[contig_id] = seqlen
                        total_size += seqlen
    except Exception:
        pass
    return contig_lengths, total_size


def parse_trna_lengths(trna_path):
    """Extract tRNA lengths grouped by integrase identifier.

    Args:
        trna_path: Path to ``integrase_trna.tsv`` (preferred) or a legacy
            table with ``mge_id``, ``start``, and ``end`` columns.

    Returns:
        Dictionary mapping integrase or MGE identifier to a list of tRNA
        lengths in base pairs.
    """
    lengths = {}
    try:
        df = pd.read_csv(trna_path, sep="\t")
        if {"integrase_id", "trna_start", "trna_end"}.issubset(df.columns):
            for _, row in df.iterrows():
                integrase_id = str(row["integrase_id"])
                trna_len = abs(int(row["trna_end"]) - int(row["trna_start"])) + 1
                lengths.setdefault(integrase_id, []).append(trna_len)
        elif {"mge_id", "start", "end"}.issubset(df.columns):
            for _, row in df.iterrows():
                mge_id = str(row["mge_id"])
                trna_len = abs(int(row["end"]) - int(row["start"])) + 1
                lengths.setdefault(mge_id, []).append(trna_len)
    except Exception:
        pass
    return lengths


def _integrase_ids_from_mge_headers(mge_lengths):
    """Extract integrase identifiers embedded in MGE FASTA headers.

    MGE headers use the format ``integrase_id:contig:start-end``; only the
    first colon-separated field is returned.

    Args:
        mge_lengths: Dictionary keyed by MGE FASTA sequence identifier.

    Returns:
        Set of integrase identifiers.
    """
    integrase_ids = set()
    for mge_id in mge_lengths:
        parts = mge_id.split(":")
        if parts:
            integrase_ids.add(parts[0])
    return integrase_ids


def main(results_dir, out_prefix):
    """Write per-genome, per-element, and cohort-level summary tables.

    Args:
        results_dir: Root directory containing per-sample pipeline outputs.
        out_prefix: Path prefix for three output TSV files:
            ``{out_prefix}_summary.tsv``, ``{out_prefix}_details.tsv``, and
            ``{out_prefix}_overall.tsv``.
    """
    samples = get_sample_list(results_dir)

    summary = []
    details = []

    all_trna_lengths = []
    all_mge_lengths = []
    all_blast_lengths = []

    for sample in samples:
        integrase_file = collect_file(results_dir, sample, "integrase_hits_summary.tsv")
        trna_file = collect_file(results_dir, sample, "integrase_trna.tsv")
        blast_file = collect_file(results_dir, sample, "mge_blast.tsv")
        region_file = collect_file(results_dir, sample, "mge_region.fa")

        integrase_total = len(parse_summary_file(integrase_file)) if integrase_file else 0

        mge_lengths = parse_fasta_lengths(region_file) if region_file else {}
        blast_lengths = parse_blast_lengths(blast_file) if blast_file else {}
        trna_lengths_all = parse_trna_lengths(trna_file) if trna_file else {}

        mge_integrase_ids = _integrase_ids_from_mge_headers(mge_lengths)
        trna_lengths_mge = {
            k: v for k, v in trna_lengths_all.items() if k in mge_integrase_ids
        }

        gff_file = collect_file(results_dir, sample, "orfs.gff")
        genome_size = 0
        if gff_file:
            _, genome_size = parse_gff(gff_file)

        filtered_mge_lengths = {
            k: v for k, v in mge_lengths.items() if v >= MIN_MGE_SIZE_BP
        }
        filtered_integrase_ids = _integrase_ids_from_mge_headers(filtered_mge_lengths)
        trna_lengths_mge = {
            k: v for k, v in trna_lengths_all.items() if k in filtered_integrase_ids
        }

        elements_total = len(filtered_mge_lengths)
        trna_total_all = sum(len(v) for v in trna_lengths_all.values())
        trna_total_mge = sum(len(v) for v in trna_lengths_mge.values())
        blast_total = sum(len(v) for v in blast_lengths.values())

        sample_trna_lengths_all = []
        for values in trna_lengths_all.values():
            sample_trna_lengths_all.extend(values)
        trna_mean_len_all = (
            sum(sample_trna_lengths_all) / len(sample_trna_lengths_all)
            if sample_trna_lengths_all
            else 0
        )

        sample_trna_lengths_mge = []
        for values in trna_lengths_mge.values():
            sample_trna_lengths_mge.extend(values)
        trna_mean_len_mge = (
            sum(sample_trna_lengths_mge) / len(sample_trna_lengths_mge)
            if sample_trna_lengths_mge
            else 0
        )

        sample_blast_lengths = []
        for values in blast_lengths.values():
            sample_blast_lengths.extend(values)
        blast_mean_len = (
            sum(sample_blast_lengths) / len(sample_blast_lengths)
            if sample_blast_lengths
            else 0
        )

        sample_mge_lengths = list(filtered_mge_lengths.values())
        mge_mean_len = (
            sum(sample_mge_lengths) / len(sample_mge_lengths)
            if sample_mge_lengths
            else 0
        )

        all_mge_lengths.extend(sample_mge_lengths)
        all_trna_lengths.extend(sample_trna_lengths_all)
        all_blast_lengths.extend(sample_blast_lengths)

        summary.append({
            "genome": sample,
            "genome_size": genome_size,
            "integrases_total": integrase_total,
            "trna_total_all": trna_total_all,
            "trna_total_mge": trna_total_mge,
            "trna_mean_len_all": trna_mean_len_all,
            "trna_mean_len_mge": trna_mean_len_mge,
            "blast_hits_total": blast_total,
            "blast_mean_len": blast_mean_len,
            "elements_total": elements_total,
            "mge_mean_len": mge_mean_len,
        })

        for mge_id, mge_len in filtered_mge_lengths.items():
            parts = mge_id.split(":")
            integrase_id = parts[0] if parts else ""
            trna_lens = trna_lengths_mge.get(integrase_id, []) if integrase_id else []
            blast_lens = blast_lengths.get(integrase_id, []) if integrase_id else []
            details.append({
                "genome": sample,
                "mge_id": mge_id,
                "mge_len": mge_len,
                "trna_len": ",".join(map(str, trna_lens)) if trna_lens else "",
                "blast_hit_len": ",".join(map(str, blast_lens)) if blast_lens else "",
            })

    summary_df = pd.DataFrame(summary)
    summary_out = f"{out_prefix}_summary.tsv"
    summary_df.to_csv(summary_out, sep="\t", index=False)
    print(f"Per-genome summary written to {summary_out}")

    details_df = pd.DataFrame(details)
    details_out = f"{out_prefix}_details.tsv"
    details_df.to_csv(details_out, sep="\t", index=False)
    print(f"Per-MGE details written to {details_out}")

    overall = {
        "total_genomes": len(samples),
        "total_trna": len(all_trna_lengths),
        "trna_mean_len": (
            sum(all_trna_lengths) / len(all_trna_lengths) if all_trna_lengths else 0
        ),
        "total_mge": len(all_mge_lengths),
        "mge_mean_len": (
            sum(all_mge_lengths) / len(all_mge_lengths) if all_mge_lengths else 0
        ),
        "total_blast_hits": len(all_blast_lengths),
        "blast_mean_len": (
            sum(all_blast_lengths) / len(all_blast_lengths) if all_blast_lengths else 0
        ),
    }
    overall_df = pd.DataFrame.from_dict(overall, orient="index", columns=["value"]).reset_index()
    overall_df.columns = ["metric", "value"]
    overall_out = f"{out_prefix}_overall.tsv"
    overall_df.to_csv(overall_out, sep="\t", index=False)
    print(f"Overall summary written to {overall_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect MGE statistics with per-element detail tables."
    )
    parser.add_argument("--results", required=True, help="Path to the results directory.")
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output path prefix for summary, details, and overall TSV files.",
    )
    args = parser.parse_args()
    main(args.results, args.out_prefix)
