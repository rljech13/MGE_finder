#!/usr/bin/env python3
"""Build IE atlas for Stage 5 (strand-aware attL, no downstream filters)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import FF, REPO

ART = FF / "artifacts"
SCRIPTS = FF / "02_atlas"
NB7 = REPO / "results_organized/genome_clustering/project/notebook7_extended"

PATH_S5 = ART / "stage5_attL_lengths.tsv"
PATH_FULL = NB7 / "ie_instance_classification_full25591_strict.tsv"
_truth = os.environ.get("MGE_INTEGRAS_TRUTH")
if not _truth:
    raise SystemExit(
        "MGE_INTEGRAS_TRUTH is not set. Point it at the integrase truth table, "
        "or run atlas_pipeline/classify_cohort.py with --truth."
    )
PATH_TRUTH = Path(_truth)
PATH_REP = REPO / "results_organized/merged_deduplicated/rep_index.tsv"
PATH_LOCUS = NB7 / "locus_eggnog_merged.tsv"
PATH_CFG = FF / "02_atlas" / "ie_classifier" / "ie_atlas_labels.yaml"
PATH_GBFF = REPO / "results_organized/merged_deduplicated"

OUT_TRUTH = ART / "stage5_truth_for_classifier.tsv"
OUT_CLASSIF = ART / "stage5_atlas_classification.tsv"
OUT_JOIN = ART / "stage5_atlas_join.tsv"
OUT_MODE = ART / "stage5_attL_length_mode.tsv"
THRESH = 32


def _qseqid_to_element_id(qseqid: str) -> str:
    s = str(qseqid or "").strip()
    if not s:
        return ""
    parts = s.rsplit(":", 1)
    if len(parts) == 2 and parts[1] in {"+", "-"}:
        return parts[0]
    return s


def load_stage5() -> pd.DataFrame:
    df = pd.read_csv(PATH_S5, sep="\t", dtype=str)
    df = df[df["status"].astype(str) == "ok"].copy()
    df["integrase_id"] = df["integrase_id"].astype(str)
    df = df.drop_duplicates("integrase_id", keep="first")
    df["attL_len_bp"] = pd.to_numeric(df["best_attl_len_stage5"], errors="coerce")
    df["attL_length_mode"] = df["attL_len_bp"].apply(
        lambda x: "short" if pd.notna(x) and x <= THRESH else ("long" if pd.notna(x) else "")
    )
    return df


def build_truth(stage5: pd.DataFrame) -> pd.DataFrame:
    truth = pd.read_csv(PATH_TRUTH, sep="\t", dtype=str, low_memory=False)
    rep = pd.read_csv(PATH_REP, sep="\t", dtype=str)
    rep["cluster_key"] = rep["cluster_id"].astype(str).str.split("|sample=", regex=False).str[0]
    rep["integrase_id"] = rep["cluster_key"].str.split(":", n=1, regex=False).str[0]
    rep_lut = rep.drop_duplicates("integrase_id", keep="first").set_index("integrase_id")

    s5_ids = set(stage5["integrase_id"])
    via_truth = truth[truth["integrase_id_lookup"].astype(str).isin(s5_ids)].copy()
    have_truth = set(via_truth["integrase_id_lookup"].astype(str))

    synth_rows = []
    for iid in sorted(s5_ids - have_truth):
        if iid not in rep_lut.index:
            continue
        r = rep_lut.loc[iid]
        s5_row = stage5.loc[stage5["integrase_id"] == iid].iloc[0]
        synth_rows.append(
            {
                "matrix_row_genome_id": str(r["rep_id"]),
                "sample_id": str(r.get("sample_id", s5_row.get("sample", ""))),
                "leiden_hd_cluster": "",
                "broad_mge_type": "other",
                "mge_canonical_id": str(r["cluster_key"]),
                "integrase_id_lookup": iid,
                "integrase_model": "",
                "integrase_strand": "",
                "integrase_contig": iid.rsplit("_", 1)[0] if "_" in iid else iid,
                "integrase_rows_in_trna_file": "1",
                "tRNA_type_assigned": "",
                "finder_data_missing": "",
                "ie_len_nt": "",
            }
        )
    synth = pd.DataFrame(synth_rows)
    out_cols = list(truth.columns)
    for c in out_cols:
        if c not in via_truth.columns:
            via_truth[c] = ""
        if c not in synth.columns:
            synth[c] = ""
    out = pd.concat([via_truth[out_cols], synth[out_cols]], ignore_index=True)
    out = out.drop_duplicates(subset=["integrase_id_lookup"], keep="first")
    OUT_TRUTH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TRUTH, sep="\t", index=False)
    print(f"truth-backed: {len(via_truth):,}  synthetic: {len(synth):,}  -> {OUT_TRUTH}")
    return out


def classify_missing(truth_s5: pd.DataFrame, existing_ids: set[str]) -> pd.DataFrame:
    need = truth_s5[~truth_s5["integrase_id_lookup"].astype(str).isin(existing_ids)].copy()
    if need.empty:
        return pd.DataFrame()
    print(f"classifying {len(need):,} IE not in full25591 atlas...")

    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    from ie_classifier.classify import classify_truth_row
    from ie_classifier.config import load_config
    from ie_classifier.loaders import build_rep_locus_file_lookup, load_merge_table
    from ie_classifier.schema import TSV_COLUMNS, flatten_record

    cfg = load_config(PATH_CFG)
    merge_df = load_merge_table(PATH_LOCUS)
    rep_df = pd.read_csv(PATH_REP, sep="\t", dtype=str)
    rep_lut = build_rep_locus_file_lookup(rep_df)

    rows = []
    total = len(need)
    for i, (_, row) in enumerate(need.iterrows(), start=1):
        rec = classify_truth_row(
            row, cfg, merge_df, rep_lut, PATH_GBFF,
            txsscan_hits_lut=None, external_hits_lut=None,
        )
        rows.append(flatten_record(rec))
        if i == 1 or i % 250 == 0 or i == total:
            print(f"  progress {i}/{total}")

    return pd.DataFrame(rows, columns=TSV_COLUMNS)


def stub_no_rep(stage5: pd.DataFrame, classified_ids: set[str]) -> pd.DataFrame:
    missing = stage5[~stage5["integrase_id"].isin(classified_ids)].copy()
    if missing.empty:
        return pd.DataFrame()
    rows = []
    for _, r in missing.iterrows():
        eid = _qseqid_to_element_id(r.get("qseqid", "")) or str(r["integrase_id"])
        rows.append(
            {
                "element_id": eid,
                "genome_id": "",
                "sample_id": str(r.get("sample", "")),
                "primary_label": "other.unclassified_IE",
                "coarse_label": "unclassified",
                "confidence": "low",
                "error": "no_rep_index_match_for_cluster",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    stage5 = load_stage5()
    print(f"Stage 5 universe: {len(stage5):,} IE")

    truth_s5 = build_truth(stage5)

    full = pd.read_csv(PATH_FULL, sep="\t", dtype=str, low_memory=False)
    full["integrase_id"] = full["element_id"].astype(str).str.split(":").str[0]
    full_s5 = full[full["integrase_id"].isin(set(stage5["integrase_id"]))].copy()
    existing_ids = set(full_s5["integrase_id"])
    print(f"from full25591 atlas: {len(full_s5):,}")

    new_cls = classify_missing(truth_s5, existing_ids)
    if not new_cls.empty:
        print(f"newly classified: {len(new_cls):,}")

    classified_ids = set(existing_ids)
    parts = [full_s5]
    if not new_cls.empty:
        new_cls = new_cls.copy()
        new_cls["integrase_id"] = new_cls["element_id"].astype(str).str.split(":").str[0]
        classified_ids |= set(new_cls["integrase_id"])
        parts.append(new_cls)

    stubs = stub_no_rep(stage5, classified_ids)
    if not stubs.empty:
        stubs["integrase_id"] = stubs["element_id"].astype(str).str.split(":").str[0]
        print(f"no-rep stubs: {len(stubs):,}")
        parts.append(stubs)

    classif = pd.concat(parts, ignore_index=True, sort=False)
    classif["integrase_id"] = classif["integrase_id"].astype(str)
    classif = classif.drop_duplicates("integrase_id", keep="first")

    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    from ie_classifier.rules import coarse_label_from_primary

    classif["coarse_label"] = classif["primary_label"].astype(str).map(coarse_label_from_primary)

    classif.to_csv(OUT_CLASSIF, sep="\t", index=False)
    print(f"-> {OUT_CLASSIF} ({len(classif):,} rows)")

    join = stage5.merge(classif, on="integrase_id", how="left", suffixes=("", "_cls"))
    join.to_csv(OUT_JOIN, sep="\t", index=False)
    mode_cols = [
        "integrase_id", "sample", "best_attl_len_stage5", "attL_len_bp", "attL_length_mode",
        "coarse_label", "primary_label", "phage_score", "conjugation_score", "replication_score",
        "mobilome_score", "island_module_score", "gta_score", "pici_score", "error",
    ]
    mode_cols = [c for c in mode_cols if c in join.columns]
    join[mode_cols].to_csv(OUT_MODE, sep="\t", index=False)
    print(f"-> {OUT_JOIN} ({len(join):,} rows)")
    print(f"-> {OUT_MODE}")
    print(join["coarse_label"].value_counts().head(12).to_string())


if __name__ == "__main__":
    main()
