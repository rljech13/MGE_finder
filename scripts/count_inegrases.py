import os
import glob
import argparse
import pandas as pd

def collect_summary_files(results_dir):
    return sorted(glob.glob(os.path.join(results_dir, "*", "integrase_hits_summary.tsv")))

def merge_summaries(files):
    dfs = []
    for path in files:
        sample = os.path.basename(os.path.dirname(path))
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "sample", sample)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

def main(results_dir, output_file):
    summary_files = collect_summary_files(results_dir)
    if not summary_files:
        print("Нет файлов интеграз для объединения.")
        return

    merged_df = merge_summaries(summary_files)
    merged_df.to_csv(output_file, sep="\t", index=False)
    print("Объединённый файл записан в:", output_file)
    print(merged_df.to_markdown(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate integrase counts from all genomes")
    parser.add_argument("--results", required=True, help="Path to results/ directory")
    parser.add_argument("--out", required=True, help="Path to output .tsv file")
    args = parser.parse_args()

    main(args.results, args.out)