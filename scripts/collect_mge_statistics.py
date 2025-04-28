import os
import glob
import argparse
import re
import pandas as pd


def collect_file(results_dir, sample, filename):
    """Collect a file path for a given sample."""
    path = os.path.join(results_dir, sample, filename)
    return path if os.path.exists(path) else None


def get_sample_list(results_dir):
    """Retrieve a sorted list of sample folder names from the results directory."""
    samples = [os.path.basename(d) for d in glob.glob(os.path.join(results_dir, "*")) if os.path.isdir(d)]
    return sorted(samples)


def parse_summary_file(file_path):
    """Parse a TSV file and return the number of rows."""
    try:
        df = pd.read_csv(file_path, sep="\t")
        return len(df)
    except Exception:
        return 0


def parse_trna_all(trna_path):
    """
    Count all tRNA entries in the ARAGORN output trna.tsv.
    Lines describing tRNAs start with a number followed by 'tRNA-'.
    """
    count = 0
    try:
        with open(trna_path) as f:
            for line in f:
                if re.match(r"^\s*\d+\s+tRNA-", line):
                    count += 1
    except Exception:
        pass
    return count


def parse_fasta_count(fasta_path):
    """Count the number of records in a FASTA file by counting '>' lines."""
    try:
        return sum(1 for line in open(fasta_path) if line.startswith('>'))
    except Exception:
        return 0


def main(results_dir, output_file, global_out_file=None):
    samples = get_sample_list(results_dir)
    num_genomes = len(samples)
    summary = []

    # totals
    total_integrases = 0
    total_integrases_in = 0
    total_trna = 0
    total_trna_near = 0
    total_trna_in = 0
    total_blast = 0
    total_blast_in = 0
    total_elements = 0

    for sample in samples:
        # file paths
        integrase_file = collect_file(results_dir, sample, "integrase_hits_summary.tsv")
        trna_all_file = collect_file(results_dir, sample, "trna.tsv")
        trna_near_file = collect_file(results_dir, sample, "integrase_trna.tsv")
        blast_file = collect_file(results_dir, sample, "mge_blast.tsv")
        region_file = collect_file(results_dir, sample, "mge_region.fa")

        # counts
        integrase_total = parse_summary_file(integrase_file) if integrase_file else 0
        trna_total = parse_trna_all(trna_all_file) if trna_all_file else 0
        trna_near = parse_summary_file(trna_near_file) if trna_near_file else 0
        blast_total = parse_summary_file(blast_file) if blast_file else 0
        elements_total = parse_fasta_count(region_file) if region_file else 0

        # define in-elements counts as number of elements
        integ_in = elements_total
        trna_in = elements_total
        blast_in = elements_total

        # accumulate
        total_integrases += integrase_total
        total_integrases_in += integ_in
        total_trna += trna_total
        total_trna_near += trna_near
        total_trna_in += trna_in
        total_blast += blast_total
        total_blast_in += blast_in
        total_elements += elements_total

        summary.append({
            "genome": sample,
            "integrases_total": integrase_total,
            "integrases_in_elements": integ_in,
            "trna_total": trna_total,
            "trna_near_integrase": trna_near,
            "trna_in_elements": trna_in,
            "blast_hits_total": blast_total,
            "blast_hits_in_elements": blast_in,
            "elements_total": elements_total
        })

    # per-genome summary
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_file, sep="\t", index=False)
    print("Per-genome summary:")
    print(summary_df.to_markdown(index=False))

    # overall totals
    overall = {
        "total_genomes": num_genomes,
        "total_integrases": total_integrases,
        "integrases_in_elements": total_integrases_in,
        "total_trna": total_trna,
        "trna_near_integrase": total_trna_near,
        "trna_in_elements": total_trna_in,
        "total_blast_hits": total_blast,
        "blast_hits_in_elements": total_blast_in,
        "total_elements": total_elements
    }
    print("\nOverall totals:")
    for k, v in overall.items():
        print(f"{k}: {v}")

    # write overall summary file
    if not global_out_file:
        global_out_file = output_file.replace('.tsv', '_overall.tsv')
    overall_df = pd.DataFrame.from_dict(overall, orient='index', columns=['count']).reset_index()
    overall_df.columns = ['metric', 'count']
    overall_df.to_csv(global_out_file, sep="\t", index=False)
    print(f"\nOverall summary written to {global_out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect detailed MGE statistics" )
    parser.add_argument("--results", required=True, help="Path to the results directory")
    parser.add_argument("--out", required=True, help="Path to save per-genome summary TSV file")
    parser.add_argument("--global-out", help="(Optional) Path to save overall summary TSV file")
    args = parser.parse_args()
    main(args.results, args.out, args.global_out)
