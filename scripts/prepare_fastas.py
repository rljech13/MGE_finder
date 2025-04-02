import os
import glob
import shutil
import argparse
import yaml
from Bio import SeqIO
from logger import Logger

log = Logger(name="prepare_fasta", draw_progress=True).get_logger()


def find_files(base_dir, pattern, recursive=True):
    """Recursively find files in base_dir matching the given pattern.

    Args:
        base_dir (str): The directory to search in.
        pattern (str): The glob pattern to match file names.
        recursive (bool): Whether to search subdirectories recursively.

    Returns:
        list: List of file paths that match the pattern.
    """
    if recursive:
        matches = []
        for root, _, filenames in os.walk(base_dir):
            for filename in filenames:
                if glob.fnmatch.fnmatch(filename, pattern):
                    matches.append(os.path.join(root, filename))
        return matches
    else:
        return glob.glob(os.path.join(base_dir, pattern))


def process_input_sources(sources, out_dir, ncbi_pattern="*_genomic.fna", hybracter_pattern="barcode*.fastq_final.fasta"):
    """Process input source directories and prepare genome FASTA files.

    Creates the output directory if it does not exist, then for each source directory,
    it searches for files matching the provided NCBI and Hybracter patterns. NCBI files
    are copied directly while Hybracter files are converted to FASTA format.

    Args:
        sources (list of str): List of directories containing input files.
        out_dir (str): Directory where processed FASTA files will be stored.
        ncbi_pattern (str): Glob pattern for NCBI FASTA files.
        hybracter_pattern (str): Glob pattern for Hybracter FASTQ files.

    Returns:
        list of str: List of sample identifiers that have been processed.
    """
    os.makedirs(out_dir, exist_ok=True)
    processed = []
    task = None

    total_files = 0
    for src in sources:
        total_files += len(find_files(src, ncbi_pattern))
        total_files += len(find_files(src, hybracter_pattern))

    if total_files > 0:
        task = Logger().progress_task("Copying genomes", total=total_files)

    for src in sources:
        # Process NCBI-style files
        fna_files = find_files(src, ncbi_pattern)
        for fna in fna_files:
            # Extract sample name from parent directory (e.g., "GCA_000008125.1")
            sample = os.path.basename(os.path.dirname(fna))
            dst = os.path.join(out_dir, f"{sample}.fna")
            try:
                shutil.copy(fna, dst)
                if sample not in processed:
                    processed.append(sample)
                log.info(f"[NCBI] {sample} → {dst}")
            except Exception as e:
                log.error(f"Error in copying {fna}: {e}")
            Logger().advance_progress(task)

        # Process Hybracter-style files
        fastq_files = find_files(src, hybracter_pattern)
        for fasta in fastq_files:
            # Extract sample name from file name (e.g., "barcode01.fastq_final.fasta" -> "barcode01")
            sample = os.path.basename(fasta).split(".")[0]
            dst = os.path.join(out_dir, f"{sample}.fna")
            try:
                with open(fasta) as fin, open(dst, "w") as fout:
                    count = SeqIO.write(SeqIO.parse(fin, "fasta"), fout, "fasta")
                if sample not in processed:
                    processed.append(sample)
                log.info(f"[HYBRACTER] {sample} ({count} seqs) → {dst}")
            except Exception as e:
                log.error(f"Error in conversion {fasta}: {e}")
            Logger().advance_progress(task)

    Logger().finish_progress(task)
    log.info(f"Totally processed: {len(processed)} samples.")
    log.info("Processed sample names: " + ", ".join(processed))
    return processed


def main(config_path, output_done=None):
    """Main function to prepare genome FASTA files.

    Loads the configuration from a YAML file, processes the input source directories
    using the provided custom file patterns (if any), and writes a completion file if specified.

    Args:
        config_path (str): Path to the configuration YAML file.
        output_done (str, optional): Path to the output completion file (used in Snakemake).

    Returns:
        None
    """
    log.info(f"Load config: {config_path}")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    sources = config["input_sources"]
    out_dir = config["genomes_dir"]

    # Get customizable patterns with defaults if not provided
    ncbi_pattern = config.get("ncbi_pattern", "*_genomic.fna")
    hybracter_pattern = config.get("hybracter_pattern", "barcode*.fastq_final.fasta")

    processed = process_input_sources(sources, out_dir, ncbi_pattern, hybracter_pattern)

    if output_done:
        with open(output_done, "w") as f:
            f.write("done\n")
        log.info(f"Completion file written: {output_done}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare genome FASTA files")
    parser.add_argument("--config", required=False, default="config.yaml", help="Path to config.yaml")
    parser.add_argument("--done", required=False, help="(Optional) Path to .complete output file (used in Snakemake)")
    args = parser.parse_args()
    main(args.config, args.done)
else:
    main("config.yaml", snakemake.output[0])