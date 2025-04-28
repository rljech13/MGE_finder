#!/usr/bin/env python3
"""
Summarize Padloc defense systems across samples.
Produces per-sample and overall frequency tables of detected systems.
"""
import os
import glob
import argparse
import pandas as pd

def get_padloc_files(results_dir):
    """Recursively find all results_padloc.csv files under results_dir."""
    pattern = os.path.join(results_dir, "**", "results_padloc.csv")
    return glob.glob(pattern, recursive=True)


def main(padloc_dir, out_per_sample, out_overall):
    # Collect all padloc result files
    padloc_files = get_padloc_files(padloc_dir)
    if not padloc_files:
        print(f"No Padloc files found under {padloc_dir}")
        return

    # Per-sample and overall aggregation
    per_sample_records = []
    overall_counts = {}

    for file in padloc_files:
        sample = os.path.basename(os.path.dirname(file))
        # Read CSV if not empty
        try:
            df = pd.read_csv(file)
        except pd.errors.EmptyDataError:
            # Skip empty files
            continue
        # Count occurrences of each system
        counts = df['system'].value_counts().to_dict()
        for system, cnt in counts.items():
            per_sample_records.append({
                'sample': sample,
                'system': system,
                'count': cnt
            })
            overall_counts[system] = overall_counts.get(system, 0) + cnt

    # Write per-sample summary
    per_sample_df = pd.DataFrame(per_sample_records)
    per_sample_df.sort_values(['sample', 'count'], ascending=[True, False], inplace=True)
    per_sample_df.to_csv(out_per_sample, sep='\t', index=False)
    print("Per-sample defense system counts:")
    print(per_sample_df.to_markdown(index=False))

    # Write overall summary sorted by frequency
    overall_list = [{'system': sys, 'count': c} for sys, c in overall_counts.items()]
    overall_df = pd.DataFrame(overall_list).sort_values('count', ascending=False)
    overall_df.to_csv(out_overall, sep='\t', index=False)
    print("\nOverall defense system frequencies:")
    print(overall_df.to_markdown(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Summarize Padloc defense systems across samples'
    )
    parser.add_argument(
        '--padloc-results', required=True,
        help='Base directory containing Padloc results subfolders'
    )
    parser.add_argument(
        '--out-sample', required=True,
        help='Path to save per-sample summary TSV'
    )
    parser.add_argument(
        '--out-overall', required=True,
        help='Path to save overall summary TSV'
    )
    args = parser.parse_args()
    main(args.padloc_results, args.out_sample, args.out_overall)
