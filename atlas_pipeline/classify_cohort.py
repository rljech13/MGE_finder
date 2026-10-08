#!/usr/bin/env python3
"""Classify a cohort the user chooses.

Every input path is an argument. The marker rules ship with this package;
the cohort, the locus table, and the GenBank files do not.

  python atlas_pipeline/classify_cohort.py \
      --cohort my_ids.tsv \
      --truth integrase_truth.tsv \
      --rep-index rep_index.tsv \
      --locus-table locus_eggnog_merged.tsv \
      --gbff-dir path/to/gbff \
      --out classification.tsv

``my_ids.tsv`` has a column ``integrase_id``. Rows absent from ``--truth``
are listed in the output with ``error=not_in_truth_table`` and are not labelled.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ie_classifier.classify import classify_truth_row
from ie_classifier.config import load_config
from ie_classifier.loaders import build_rep_locus_file_lookup, load_merge_table

DEFAULT_LABELS = HERE / "ie_classifier" / "ie_atlas_labels.yaml"


def classify_cohort(
    cohort_ids: list[str],
    truth: pd.DataFrame,
    cfg: dict,
    merge_df: pd.DataFrame,
    rep_lut: dict,
    gbff_dir: Path,
) -> pd.DataFrame:
    by_id = truth.drop_duplicates("integrase_id_lookup").set_index("integrase_id_lookup", drop=False)
    rows = []
    for iid in cohort_ids:
        if iid not in by_id.index:
            rows.append({
                "integrase_id": iid,
                "genome_id": "",
                "primary_label": "",
                "coarse_label": "",
                "confidence": "",
                "n_cds": "",
                "rep_gbff": "",
                "error": "not_in_truth_table",
            })
            continue
        rec = classify_truth_row(by_id.loc[iid], cfg, merge_df, rep_lut, gbff_dir)
        rows.append({
            "integrase_id": iid,
            "genome_id": rec.get("genome_id", ""),
            "primary_label": rec.get("primary_label", ""),
            "coarse_label": rec.get("coarse_label", ""),
            "confidence": rec.get("confidence", ""),
            "n_cds": rec.get("n_cds", ""),
            "rep_gbff": rec.get("rep_gbff", ""),
            "error": rec.get("error", ""),
        })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Classify integrative elements from a cohort table.")
    ap.add_argument("--cohort", required=True, type=Path, help="TSV with column integrase_id")
    ap.add_argument("--truth", required=True, type=Path, help="Truth table with integrase_id_lookup")
    ap.add_argument("--rep-index", required=True, type=Path, help="rep_index.tsv")
    ap.add_argument("--locus-table", required=True, type=Path, help="Per-locus eggNOG table")
    ap.add_argument("--gbff-dir", required=True, type=Path, help="Directory of representative GBFF files")
    ap.add_argument("--labels", type=Path, default=DEFAULT_LABELS, help="Marker yaml")
    ap.add_argument("--out", required=True, type=Path, help="Classification TSV to write")
    args = ap.parse_args(argv)

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    if "integrase_id" not in cohort.columns:
        raise SystemExit(f"{args.cohort} has no integrase_id column")
    ids = cohort["integrase_id"].astype(str).tolist()
    truth = pd.read_csv(args.truth, sep="\t", dtype=str, low_memory=False)
    cfg = load_config(args.labels)
    rep = pd.read_csv(args.rep_index, sep="\t", dtype=str)
    rep_lut = build_rep_locus_file_lookup(rep)
    merge_df = load_merge_table(args.locus_table)
    out = classify_cohort(ids, truth, cfg, merge_df, rep_lut, args.gbff_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    n_err = int((out["error"].astype(str).str.len() > 0).sum())
    print(f"wrote {args.out}  rows={len(out)}  errors={n_err}")
    labelled = out[out["coarse_label"].astype(str).str.len() > 0]
    if not labelled.empty:
        print(labelled["coarse_label"].value_counts().to_string())


if __name__ == "__main__":
    main()
