import os
import glob
import argparse
import pandas as pd

def collect_file(results_dir, sample, filename):
    """Collect a file path for a given sample.

    Constructs the file path from results_dir, sample, and filename. Returns the path if the file exists,
    otherwise returns None.

    Args:
        results_dir (str): Base directory where sample folders are located.
        sample (str): Sample identifier (folder name).
        filename (str): Name of the file to look for.

    Returns:
        str or None: The full file path if it exists; otherwise, None.
    """
    path = os.path.join(results_dir, sample, filename)
    return path if os.path.exists(path) else None

def get_sample_list(results_dir):
    """Retrieve a sorted list of sample folder names from the results directory.

    Args:
        results_dir (str): Directory containing sample subdirectories.

    Returns:
        list of str: Sorted list of sample identifiers.
    """
    samples = [os.path.basename(d) for d in glob.glob(os.path.join(results_dir, "*")) if os.path.isdir(d)]
    return sorted(samples)

def parse_summary_file(file_path):
    """Parse a summary TSV file and return the number of rows.

    Reads a TSV file into a DataFrame and returns the number of rows. If any error occurs, returns 0.

    Args:
        file_path (str): Path to the TSV file.

    Returns:
        int: Number of rows in the TSV file, or 0 if an error occurs.
    """
    try:
        df = pd.read_csv(file_path, sep="\t")
        return len(df)
    except Exception:
        return 0

def main(results_dir, output_file):
    """Collect summary information from each sample folder and save it as a TSV file.

    This function obtains a sorted list of sample folders from results_dir, then for each sample,
    it collects counts from integrase_hits_summary.tsv, integrase_trna.tsv, and mge_blast.tsv.
    It aggregates the counts for all samples and writes a summary TSV file.
    Finally, it prints a markdown summary and overall totals.

    Args:
        results_dir (str): Path to the results directory containing sample folders.
        output_file (str): Path to save the final summary TSV file.

    Returns:
        None
    """
    samples = get_sample_list(results_dir)
    summary = []
    
    total_integrases = 0
    total_trna = 0
    total_blast = 0
    
    for sample in samples:
        integrase_file = collect_file(results_dir, sample, "integrase_hits_summary.tsv")
        trna_file = collect_file(results_dir, sample, "integrase_trna.tsv")
        blast_file = collect_file(results_dir, sample, "mge_blast.tsv")
        
        integrase_count = parse_summary_file(integrase_file) if integrase_file else 0
        trna_count = parse_summary_file(trna_file) if trna_file else 0
        blast_count = parse_summary_file(blast_file) if blast_file else 0
        
        total_integrases += integrase_count
        total_trna += trna_count
        total_blast += blast_count
        
        summary.append({
            "sample": sample,
            "integrases": integrase_count,
            "trna_found": trna_count,
            "blast_hits": blast_count
        })
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_file, sep="\t", index=False)
    
    print("Summary information by sample:")
    print(summary_df.to_markdown(index=False))
    print("\nOverall totals:")
    print(f"Samples: {len(samples)}")
    print(f"Total integrases: {total_integrases}")
    print(f"Integrases with nearby tRNA: {total_trna}")
    print(f"Total BLAST hits: {total_blast}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect summary information on integrases, tRNA, and BLAST hits"
    )
    parser.add_argument("--results", required=True, help="Path to the results directory")
    parser.add_argument("--out", required=True, help="Path to save the final TSV summary file")
    args = parser.parse_args()
    
    main(args.results, args.out)