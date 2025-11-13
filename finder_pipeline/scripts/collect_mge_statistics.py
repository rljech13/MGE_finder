import os
import glob
import argparse
import re
import pandas as pd
from Bio import SeqIO


def collect_file(results_dir, sample, filename):
    """Collect a file path for a given sample."""
    path = os.path.join(results_dir, sample, filename)
    return path if os.path.exists(path) else None


def get_sample_list(results_dir):
    """Retrieve a sorted list of sample folder names from the results directory."""
    samples = [os.path.basename(d) for d in glob.glob(os.path.join(results_dir, "*")) if os.path.isdir(d)]
    return sorted(samples)


def parse_summary_file(file_path):
    """Parse a TSV file and return DataFrame (empty if fail)."""
    try:
        return pd.read_csv(file_path, sep="\t")
    except Exception:
        return pd.DataFrame()


def parse_fasta_lengths(fasta_path):
    """Return dict: MGE_id -> length (from FASTA)."""
    lengths = {}
    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            lengths[rec.id] = len(rec.seq)
    except Exception:
        pass
    return lengths


def parse_blast_lengths(blast_path):
    """
    Return dict: MGE_id -> list of hit lengths.
    Expect columns: qseqid, qstart, qend
    """
    lengths = {}
    try:
        df = pd.read_csv(blast_path, sep="\t")
        if not {"qseqid", "qstart", "qend"}.issubset(df.columns):
            return lengths
        for _, row in df.iterrows():
            mge_id = str(row["qseqid"])
            l = abs(int(row["qend"]) - int(row["qstart"])) + 1
            lengths.setdefault(mge_id, []).append(l)
    except Exception:
        pass
    return lengths


def parse_trna_lengths(trna_path):
    """
    Return dict: MGE_id -> list of tRNA lengths.
    Expect columns: mge_id, start, end
    """
    lengths = {}
    try:
        df = pd.read_csv(trna_path, sep="\t")
        if not {"mge_id", "start", "end"}.issubset(df.columns):
            return lengths
        for _, row in df.iterrows():
            mge_id = str(row["mge_id"])
            l = abs(int(row["end"]) - int(row["start"])) + 1
            lengths.setdefault(mge_id, []).append(l)
    except Exception:
        pass
    return lengths


def main(results_dir, out_prefix):
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
        trna_lengths = parse_trna_lengths(trna_file) if trna_file else {}

        elements_total = len(mge_lengths)
        trna_total = sum(len(v) for v in trna_lengths.values())
        blast_total = sum(len(v) for v in blast_lengths.values())

        # копим все длины
        all_mge_lengths.extend(mge_lengths.values())
        for v in trna_lengths.values():
            all_trna_lengths.extend(v)
        for v in blast_lengths.values():
            all_blast_lengths.extend(v)

        summary.append({
            "genome": sample,
            "integrases_total": integrase_total,
            "trna_total": trna_total,
            "trna_mean_len": sum(all_trna_lengths)/len(all_trna_lengths) if all_trna_lengths else 0,
            "blast_hits_total": blast_total,
            "blast_mean_len": sum(all_blast_lengths)/len(all_blast_lengths) if all_blast_lengths else 0,
            "elements_total": elements_total,
            "mge_mean_len": sum(all_mge_lengths)/len(all_mge_lengths) if all_mge_lengths else 0,
        })

        # детали по каждому MGE
        for mge_id, mge_len in mge_lengths.items():
            details.append({
                "genome": sample,
                "mge_id": mge_id,
                "mge_len": mge_len,
                "trna_len": ",".join(map(str, trna_lengths.get(mge_id, []))) if mge_id in trna_lengths else "",
                "blast_hit_len": ",".join(map(str, blast_lengths.get(mge_id, []))) if mge_id in blast_lengths else ""
            })

    # сохраняем summary
    summary_df = pd.DataFrame(summary)
    summary_out = f"{out_prefix}_summary.tsv"
    summary_df.to_csv(summary_out, sep="\t", index=False)
    print(f"Per-genome summary written to {summary_out}")

    # сохраняем details
    details_df = pd.DataFrame(details)
    details_out = f"{out_prefix}_details.tsv"
    details_df.to_csv(details_out, sep="\t", index=False)
    print(f"Per-MGE details written to {details_out}")

    # общий summary
    overall = {
        "total_genomes": len(samples),
        "total_trna": len(all_trna_lengths),
        "trna_mean_len": sum(all_trna_lengths)/len(all_trna_lengths) if all_trna_lengths else 0,
        "total_mge": len(all_mge_lengths),
        "mge_mean_len": sum(all_mge_lengths)/len(all_mge_lengths) if all_mge_lengths else 0,
        "total_blast_hits": len(all_blast_lengths),
        "blast_mean_len": sum(all_blast_lengths)/len(all_blast_lengths) if all_blast_lengths else 0,
    }
    overall_df = pd.DataFrame.from_dict(overall, orient='index', columns=['value']).reset_index()
    overall_df.columns = ['metric', 'value']
    overall_out = f"{out_prefix}_overall.tsv"
    overall_df.to_csv(overall_out, sep="\t", index=False)
    print(f"Overall summary written to {overall_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect MGE statistics with per-MGE details")
    parser.add_argument("--results", required=True, help="Path to the results directory")
    parser.add_argument("--out-prefix", required=True, help="Prefix for output TSV files")
    args = parser.parse_args()
    main(args.results, args.out_prefix)