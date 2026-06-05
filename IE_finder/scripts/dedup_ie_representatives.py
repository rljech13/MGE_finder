#!/usr/bin/env python3
"""Cluster confident integrative elements with MMseqs2 and write representative outputs."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from logger import Logger

logger = Logger(name="dedup_ie_representatives").get_logger()


def load_dedup_config(config_path: Path | None) -> dict:
    """Load MMseqs deduplication settings from a pipeline configuration file.

    Args:
        config_path: Path to ``ie_finder_config.yaml``, or None to use defaults.

    Returns:
        Dictionary with keys ``enabled``, ``min_seq_id``, ``min_coverage``,
        and ``mmseqs_bin``.
    """
    defaults = {
        "enabled": True,
        "min_seq_id": 0.9,
        "min_coverage": 0.8,
        "mmseqs_bin": None,
    }
    if config_path is None or not config_path.is_file():
        return defaults
    with config_path.open() as fh:
        cfg = yaml.safe_load(fh) or {}
    dedup = cfg.get("dedup", {})
    return {
        "enabled": bool(dedup.get("enabled", True)),
        "min_seq_id": float(dedup.get("min_seq_id", 0.9)),
        "min_coverage": float(dedup.get("min_coverage", 0.8)),
        "mmseqs_bin": dedup.get("mmseqs_bin"),
    }


def resolve_mmseqs(explicit: str | None) -> str:
    """Locate the MMseqs2 executable.

    Args:
        explicit: Optional path from configuration; searched first when set.

    Returns:
        Absolute path to the ``mmseqs`` binary.

    Raises:
        FileNotFoundError: When no executable can be resolved.
    """
    if explicit:
        p = Path(explicit)
        if p.is_file():
            return str(p)
        raise FileNotFoundError(f"MMseqs binary not found: {explicit}")
    found = shutil.which("mmseqs")
    if found:
        return found
    raise FileNotFoundError(
        "mmseqs not found on PATH; install mmseqs2 or set dedup.mmseqs_bin in config"
    )


def safe_filename(text: str, max_len: int = 80) -> str:
    """Convert arbitrary text into a filesystem-safe filename fragment.

    Args:
        text: Input string, typically a sample identifier.
        max_len: Maximum returned string length.

    Returns:
        Sanitized filename fragment.
    """
    out = re.sub(r'[<>:"/\\|?*\s]', "_", text)
    out = out.strip("._") or "unknown"
    return out[:max_len]


def cluster_id(ie_id: str, sample_id: str) -> str:
    """Build a cluster member identifier used in downstream metadata tables.

    Args:
        ie_id: Integrative element FASTA record identifier.
        sample_id: Sample (genome) identifier.

    Returns:
        String of the form ``ie_id|sample=sample_id``.
    """
    return f"{ie_id}|sample={sample_id}"


def discover_confident(results_dir: Path) -> list[tuple[str, Path, Path]]:
    """Find per-sample confident IE files under a Snakemake results directory.

    Args:
        results_dir: Root results directory (for example ``results/``).

    Returns:
        List of tuples ``(sample_id, ie_confident.fa path, ie_confident.gbk path)``.
    """
    found: list[tuple[str, Path, Path]] = []
    if not results_dir.is_dir():
        return found
    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name in ("combined", "dedup"):
            continue
        fa = sample_dir / "ie_confident.fa"
        gbk = sample_dir / "ie_confident.gbk"
        if fa.is_file() and fa.stat().st_size > 0:
            found.append((sample_dir.name, fa, gbk if gbk.is_file() else fa))
    return found


def build_input_fasta(entries: list[tuple[str, Path, Path]], out_fa: Path) -> pd.DataFrame:
    """Merge per-sample confident FASTA records into one MMseqs input file.

    Args:
        entries: Output of ``discover_confident``.
        out_fa: Path to the combined nucleotide FASTA written for MMseqs.

    Returns:
        Manifest DataFrame mapping numeric sequence identifiers to sample metadata.
    """
    rows: list[dict] = []
    seq_id = 0
    with out_fa.open("w") as w:
        for sample_id, fa_path, gbk_path in entries:
            for rec in SeqIO.parse(fa_path, "fasta"):
                ie_id = rec.id
                integrase_id = ie_id.split(":", 1)[0]
                seq_id += 1
                sid = str(seq_id)
                w.write(f">{sid}\n{str(rec.seq)}\n")
                rows.append(
                    {
                        "seq_db_id": sid,
                        "sample_id": sample_id,
                        "ie_id": ie_id,
                        "integrase_id": integrase_id,
                        "cluster_member": cluster_id(ie_id, sample_id),
                        "gbk_path": str(gbk_path),
                    }
                )
    return pd.DataFrame(rows)


def run_mmseqs(
    mmseqs_bin: str,
    fa_in: Path,
    cluster_prefix: Path,
    tmp_dir: Path,
    min_seq_id: float,
    min_coverage: float,
) -> None:
    """Run ``mmseqs easy-cluster`` on a nucleotide FASTA file.

    Args:
        mmseqs_bin: Path to the MMseqs2 executable.
        fa_in: Input FASTA path.
        cluster_prefix: Output prefix passed to MMseqs2.
        tmp_dir: Temporary directory for MMseqs2 intermediate files.
        min_seq_id: Minimum sequence identity (0–1).
        min_coverage: Minimum alignment coverage (0–1).

    Raises:
        subprocess.CalledProcessError: When MMseqs2 exits with a non-zero code.
    """
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)
    threads = max(1, min(16, os.cpu_count() or 4))
    cmd = [
        mmseqs_bin,
        "easy-cluster",
        str(fa_in),
        str(cluster_prefix),
        str(tmp_dir),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(min_coverage),
        "--shuffle",
        "0",
        "--threads",
        str(threads),
    ]
    logger.info("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def resolve_rep(
    sid: str,
    member_to_rep: dict[str, str],
    rep_seq_ids: set[str],
    cl: pd.DataFrame,
) -> str:
    """Map a cluster member sequence identifier to its representative identifier.

    Args:
        sid: Numeric sequence identifier from the MMseqs input FASTA.
        member_to_rep: Member-to-representative mapping from the cluster TSV.
        rep_seq_ids: Set of representative sequence identifiers.
        cl: Full cluster assignment table from MMseqs2.

    Returns:
        Representative sequence identifier.

    Raises:
        RuntimeError: When the representative cannot be resolved.
    """
    if sid in member_to_rep:
        return member_to_rep[sid]
    if sid in rep_seq_ids:
        return sid
    alt = cl.loc[cl["rep_id"] == sid, "member_id"]
    if not alt.empty:
        anchor = str(alt.iloc[0])
        return member_to_rep.get(anchor, anchor)
    raise RuntimeError(f"Cannot resolve MMseqs representative for seq_db_id={sid}")


def write_empty_outputs(out_dir: Path) -> None:
    """Create empty placeholder files when deduplication is skipped or has no input.

    Args:
        out_dir: Deduplication output directory (for example ``results/dedup``).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "ie_all_confident.fa",
        "ie_representatives.fa",
        "ie_representatives.gbk",
        "ie_cluster_manifest.tsv",
        "ie_dedup_mapping.tsv",
        "ie_cluster_metadata.tsv",
        "rep_index.tsv",
    ):
        (out_dir / name).write_text("")
    (out_dir / "representatives").mkdir(exist_ok=True)


def dedup_representatives(
    results_dir: Path,
    out_dir: Path,
    config_path: Path | None,
) -> int:
    """Cluster confident integrative elements and write representative outputs.

    Args:
        results_dir: Snakemake results root containing per-sample directories.
        out_dir: Output directory for deduplication artifacts.
        config_path: Path to ``ie_finder_config.yaml``.

    Returns:
        Number of non-redundant representative integrase identifiers.

    Raises:
        FileNotFoundError: When MMseqs2 or expected cluster output is missing.
        subprocess.CalledProcessError: When MMseqs2 fails.
    """
    cfg = load_dedup_config(config_path)
    if not cfg["enabled"]:
        logger.info("dedup.enabled=false; skipping")
        write_empty_outputs(out_dir)
        return 0

    entries = discover_confident(results_dir)
    if not entries:
        logger.warning("No confident IE sequences found; writing empty dedup outputs")
        write_empty_outputs(out_dir)
        return 0

    mmseqs_bin = resolve_mmseqs(cfg["mmseqs_bin"])
    out_dir.mkdir(parents=True, exist_ok=True)
    mmseqs_dir = out_dir / "mmseqs"
    mmseqs_dir.mkdir(exist_ok=True)

    fa_in = out_dir / "ie_all_confident.fa"
    manifest = build_input_fasta(entries, fa_in)
    n_pre = len(manifest)
    manifest.to_csv(out_dir / "ie_cluster_manifest.tsv", sep="\t", index=False)
    logger.info(f"Collected {n_pre} confident IE(s) from {len(entries)} sample(s)")

    cluster_prefix = mmseqs_dir / "ie_cluster"
    run_mmseqs(
        mmseqs_bin,
        fa_in,
        cluster_prefix,
        mmseqs_dir / "tmp",
        cfg["min_seq_id"],
        cfg["min_coverage"],
    )

    cluster_tsv = mmseqs_dir / "ie_cluster_cluster.tsv"
    if not cluster_tsv.is_file():
        raise FileNotFoundError(f"MMseqs cluster TSV missing: {cluster_tsv}")

    cl = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["member_id", "rep_id"])
    cl["member_id"] = cl["member_id"].astype(str)
    cl["rep_id"] = cl["rep_id"].astype(str)
    cl["_self"] = (cl["member_id"] == cl["rep_id"]).astype(int)
    cl_map = (
        cl.sort_values(["member_id", "_self"], ascending=[True, False])
        .drop_duplicates("member_id", keep="first")
        .drop(columns=["_self"])
    )
    member_to_rep = dict(zip(cl_map["member_id"], cl_map["rep_id"]))

    rep_seq_ids: set[str] = set()
    with (mmseqs_dir / "ie_cluster_rep_seq.fasta").open() as fh:
        for line in fh:
            if line.startswith(">"):
                rep_seq_ids.add(line[1:].strip().split()[0])

    id_to_row = manifest.set_index("seq_db_id").to_dict("index")

    mapping_rows = []
    for _, row in manifest.iterrows():
        sid = str(row["seq_db_id"])
        rep_sid = resolve_rep(sid, member_to_rep, rep_seq_ids, cl)
        rep_row = id_to_row.get(rep_sid, row)
        mapping_rows.append(
            {
                "integrase_id": row["integrase_id"],
                "ie_id": row["ie_id"],
                "sample_id": row["sample_id"],
                "seq_db_id": sid,
                "rep_seq_db_id": rep_sid,
                "rep_integrase_id": rep_row["integrase_id"],
                "rep_ie_id": rep_row["ie_id"],
                "rep_sample_id": rep_row["sample_id"],
                "is_representative": rep_sid == sid,
            }
        )
    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(out_dir / "ie_dedup_mapping.tsv", sep="\t", index=False)

    cluster_sizes = mapping.groupby("rep_seq_db_id").size().to_dict()
    metadata_rows = []
    for _, row in mapping.iterrows():
        metadata_rows.append(
            {
                "cluster_id": row["rep_ie_id"] + f"|sample={row['rep_sample_id']}",
                "member": row["ie_id"] + f"|sample={row['sample_id']}",
                "sample_id": row["sample_id"],
                "integrase_id": row["integrase_id"],
                "ie_id": row["ie_id"],
                "cluster_size": cluster_sizes.get(row["rep_seq_db_id"], 1),
                "is_representative": row["is_representative"],
            }
        )
    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(out_dir / "ie_cluster_metadata.tsv", sep="\t", index=False)

    rep_manifest = manifest[manifest["seq_db_id"].astype(str).isin(rep_seq_ids)].copy()
    rep_integrase_ids = set(mapping.loc[mapping["is_representative"], "integrase_id"].astype(str))
    n_post = len(rep_integrase_ids)
    logger.info(
        f"MMseqs representatives: {n_post}/{n_pre} "
        f"(removed {n_pre - n_post}, {100 * (n_pre - n_post) / n_pre:.2f}%)"
    )

    rep_records_fa = []
    rep_records_gbk = []
    gbk_cache: dict[tuple[str, str], dict[str, object]] = {}

    for _, row in rep_manifest.iterrows():
        sid = str(row["seq_db_id"])
        for rec in SeqIO.parse(fa_in, "fasta"):
            if rec.id == sid:
                rec.id = row["cluster_member"]
                rec.description = ""
                rep_records_fa.append(rec)
                break

        sample_id = str(row["sample_id"])
        ie_id = str(row["ie_id"])
        gbk_path = Path(str(row["gbk_path"]))
        cache_key = (sample_id, str(gbk_path))
        if cache_key not in gbk_cache:
            gbk_cache[cache_key] = {
                r.id.split(":")[0] if ":" in r.id else r.id: r
                for r in SeqIO.parse(gbk_path, "genbank")
            } if gbk_path.is_file() else {}
        integrase_key = ie_id.split(":")[0]
        gbk_rec = gbk_cache[cache_key].get(integrase_key)
        if gbk_rec is not None:
            gbk_copy = deepcopy(gbk_rec)
            gbk_copy.id = row["cluster_member"]
            gbk_copy.description = f"representative|sample={sample_id}"
            rep_records_gbk.append(gbk_copy)

    SeqIO.write(rep_records_fa, out_dir / "ie_representatives.fa", "fasta")
    SeqIO.write(rep_records_gbk, out_dir / "ie_representatives.gbk", "genbank")

    rep_dir = out_dir / "representatives"
    if rep_dir.exists():
        shutil.rmtree(rep_dir)
    rep_dir.mkdir(parents=True)
    index_rows = []
    rep_rows = metadata[metadata["is_representative"]].drop_duplicates("integrase_id")
    for idx, (_, row) in enumerate(rep_rows.iterrows(), start=1):
        integrase_id = str(row["integrase_id"])
        sample_id = str(row["sample_id"])
        gbk_match = next((r for r in rep_records_gbk if r.id.startswith(integrase_id)), None)
        if gbk_match is None:
            continue
        fname = f"rep_{idx:06d}_{safe_filename(sample_id, 60)}.gbk"
        SeqIO.write(gbk_match, rep_dir / fname, "genbank")
        index_rows.append(
            {
                "rep_id": idx,
                "file": f"representatives/{fname}",
                "cluster_id": row["member"],
                "sample_id": sample_id,
                "integrase_id": integrase_id,
                "ie_id": row["ie_id"],
                "cluster_size": int(row["cluster_size"]),
            }
        )
    pd.DataFrame(index_rows).to_csv(out_dir / "rep_index.tsv", sep="\t", index=False)
    return n_post


def main() -> None:
    """Command-line entry point for cohort-level IE deduplication."""
    parser = argparse.ArgumentParser(
        description="Cluster confident integrative elements with MMseqs2."
    )
    parser.add_argument("--results-dir", required=True, help="Snakemake results directory.")
    parser.add_argument("--out-dir", required=True, help="Deduplication output directory.")
    parser.add_argument("--config", default="ie_finder_config.yaml", help="Pipeline config.")
    args = parser.parse_args()
    n = dedup_representatives(Path(args.results_dir), Path(args.out_dir), Path(args.config))
    logger.info(f"Dedup complete: {n} representative IE(s)")


if __name__ == "__main__":
    main()
