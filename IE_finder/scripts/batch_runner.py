"""Run the finder pipeline in genome batches from a manifest file.

Reads a manifest of genome FASTA paths, stages each batch into the working
genomes directory, invokes Snakemake, and clears the staging directory before
processing the next batch. Progress is persisted so interrupted runs can resume.
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

BATCH_STATE_PREFIX = ".batch_state"


def read_manifest(manifest_path):
    """Load genome file paths from a manifest.

    Blank lines and lines starting with ``#`` are ignored.

    Args:
        manifest_path: Path to a text file with one genome path per line.

    Returns:
        List of non-empty, non-comment file paths.
    """
    with open(manifest_path) as f:
        lines = [line.strip() for line in f]
    return [line for line in lines if line and not line.startswith("#")]


def infer_sample_name(path):
    """Derive a sample identifier from a genome file path.

    Strips common FASTA extensions (``.gz``, ``.fna``, ``.fa``, ``.fasta``)
    from the basename.

    Args:
        path: Path to a genome FASTA file.

    Returns:
        Sample name used as the output ``{sample}.fna`` basename.
    """
    name = os.path.basename(path)
    for ext in (".gz", ".fna", ".fa", ".fasta"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name


def ensure_clean_dir(path):
    """Remove all contents of a directory and recreate it.

    Args:
        path: Directory to empty. The directory itself is preserved.
    """
    Path(path).mkdir(parents=True, exist_ok=True)
    for entry in Path(path).iterdir():
        if entry.is_file() or entry.is_symlink():
            entry.unlink()
        elif entry.is_dir():
            shutil.rmtree(entry)


def place_genome(src, dst, mode):
    """Copy, symlink, or decompress a genome into the staging directory.

    Gzip-compressed inputs (``.gz`` suffix) are always decompressed to plain
    FASTA regardless of ``mode``.

    Args:
        src: Source genome path.
        dst: Destination path (typically ``{work_dir}/{sample}.fna``).
        mode: Placement mode for non-gzip inputs: ``"copy"`` or ``"symlink"``.

    Raises:
        OSError: When creating a symlink or writing the destination fails.
    """
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    elif mode == "symlink":
        try:
            os.symlink(src, dst)
        except FileExistsError:
            os.remove(dst)
            os.symlink(src, dst)
    else:
        shutil.copy2(src, dst)


def save_state(state_path, data):
    """Write batch progress to a JSON state file.

    Args:
        state_path: Path to the state file.
        data: Serializable dictionary with keys ``next_index`` and ``completed``.
    """
    Path(state_path).parent.mkdir(parents=True, exist_ok=True)
    with open(state_path, "w") as f:
        json.dump(data, f, indent=2)


def load_state(state_path):
    """Load batch progress from a JSON state file.

    Args:
        state_path: Path to the state file.

    Returns:
        Parsed state dictionary, or None when the file does not exist.
    """
    if os.path.exists(state_path):
        with open(state_path) as f:
            return json.load(f)
    return None


def compute_state_path(work_dir, manifest_path):
    """Build the default path for the batch state file.

    The state file is stored next to the pipeline root (parent of ``work_dir``)
    and named from the manifest basename.

    Args:
        work_dir: Staging directory used as ``GENOMES_DIR`` during batch runs.
        manifest_path: Path to the genome manifest file.

    Returns:
        Absolute path to the JSON state file.
    """
    work_path = Path(work_dir).resolve()
    manifest_name = Path(manifest_path).stem
    state_name = f"{BATCH_STATE_PREFIX}_{manifest_name}.json"
    return work_path.parent / state_name


def run_batch(snakemake_cmd, work_dir):
    """Execute Snakemake from the directory that contains the Snakefile.

    Walks up at most five parent directories from ``work_dir`` to locate
    ``Snakefile``.

    Args:
        snakemake_cmd: Full Snakemake shell command string.
        work_dir: Staging directory passed to the pipeline (used to find Snakefile).

    Raises:
        FileNotFoundError: When no Snakefile is found within the search limit.
        RuntimeError: When Snakemake exits with a non-zero status.
    """
    current = os.path.abspath(work_dir)
    snakefile_dir = None

    for _ in range(5):
        if os.path.exists(os.path.join(current, "Snakefile")):
            snakefile_dir = current
            break
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent

    if not snakefile_dir:
        raise FileNotFoundError(f"Snakefile not found starting from {work_dir}")

    result = subprocess.run(
        snakemake_cmd,
        shell=True,
        executable="/bin/bash",
        cwd=snakefile_dir,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Snakemake failed with exit code {result.returncode}")


def build_snakemake_cmd(args):
    """Construct a default Snakemake command from parsed CLI arguments.

    Args:
        args: Namespace returned by ``argparse`` with Snakemake-related fields.

    Returns:
        Shell command string, or ``args.snakemake_cmd`` when already provided.
    """
    if args.snakemake_cmd:
        return args.snakemake_cmd

    parts = [
        "snakemake",
        "--cores",
        str(args.snakemake_cores),
        "--configfile",
        args.snakemake_config,
    ]

    if args.use_conda:
        parts.append("--use-conda")
        if args.conda_prefix:
            parts.extend(["--conda-prefix", args.conda_prefix])

    return " ".join(parts)


def process_batches(args):
    """Process all manifest entries in batches of ``args.batch_size``.

    On failure, genomes remain in ``args.work_dir`` for inspection and the
    exception is re-raised. Successful batches update the persisted state file.

    Args:
        args: Namespace returned by ``argparse`` with manifest, work directory,
            batch size, copy mode, and Snakemake settings.

    Raises:
        FileNotFoundError: When a manifest path or genome file is missing.
        RuntimeError: When Snakemake fails for a batch.
    """
    manifest = read_manifest(args.manifest)
    if not manifest:
        print("Manifest is empty", file=sys.stderr)
        return

    state_path = compute_state_path(args.work_dir, args.manifest)
    if args.restart and os.path.exists(state_path):
        os.remove(state_path)

    state = load_state(state_path) or {"next_index": 0, "completed": []}

    while state["next_index"] < len(manifest):
        batch_id = len(state["completed"]) + 1
        batch_paths = manifest[
            state["next_index"] : state["next_index"] + args.batch_size
        ]
        print(f"=== Batch {batch_id} ({len(batch_paths)} genomes) ===")

        ensure_clean_dir(args.work_dir)

        for src in batch_paths:
            if not os.path.exists(src):
                raise FileNotFoundError(f"Genome file not found: {src}")
            sample = infer_sample_name(src)
            dst = os.path.join(args.work_dir, f"{sample}.fna")
            place_genome(src, dst, args.copy_mode)

        try:
            run_batch(args.snakemake_cmd, args.work_dir)
        except Exception as exc:
            print(f"Batch {batch_id} failed: {exc}", file=sys.stderr)
            print(f"Genomes preserved in {args.work_dir} for inspection.", file=sys.stderr)
            raise

        ensure_clean_dir(args.work_dir)
        state["next_index"] += len(batch_paths)
        state["completed"].append(batch_id)
        save_state(state_path, state)
        print(f"Batch {batch_id} completed.")

    print("All batches processed.")


def main():
    """Parse command-line arguments and run batched Snakemake execution."""
    parser = argparse.ArgumentParser(description="Run the finder pipeline in batches.")
    parser.add_argument(
        "--manifest",
        required=True,
        help="Text file with genome paths (one per line).",
    )
    parser.add_argument(
        "--work-dir",
        required=True,
        help="Directory used as GENOMES_DIR during each batch run.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Number of genomes per batch (default: 1000).",
    )
    parser.add_argument(
        "--copy-mode",
        choices=["copy", "symlink"],
        default="copy",
        help="How to place non-gzip genomes into the work directory (default: copy).",
    )
    parser.add_argument(
        "--snakemake-cmd",
        default=None,
        help="Custom Snakemake command. Built automatically when omitted.",
    )
    parser.add_argument(
        "--snakemake-cores",
        type=int,
        default=16,
        help="Core count when the Snakemake command is built automatically (default: 16).",
    )
    parser.add_argument(
        "--snakemake-config",
        default="ie_finder_config.yaml",
        help="Config file passed to Snakemake when the command is built automatically.",
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help="Add --use-conda when building the Snakemake command automatically.",
    )
    parser.add_argument(
        "--conda-prefix",
        default=None,
        help="Path for --conda-prefix (only used together with --use-conda).",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        help="Ignore existing batch state and start from the first manifest entry.",
    )
    args = parser.parse_args()

    args.work_dir = os.path.abspath(args.work_dir)
    args.manifest = os.path.abspath(args.manifest)
    if args.snakemake_config:
        args.snakemake_config = os.path.abspath(args.snakemake_config)
    if args.conda_prefix:
        args.conda_prefix = os.path.abspath(args.conda_prefix)

    if not args.snakemake_cmd:
        args.snakemake_cmd = build_snakemake_cmd(args)

    process_batches(args)


if __name__ == "__main__":
    main()
