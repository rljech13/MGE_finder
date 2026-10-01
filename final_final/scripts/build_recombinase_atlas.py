#!/usr/bin/env python3
"""IE Atlas classification + join for serine recombinase cohort (489 reps).

Uses the same ie_classifier rules as tyrosine NB7, but:
  - CDS from Bakta annotation.gbff (full gene content)
  - attL/attR flags from finder dedup rep .gbk
  - locus Pfam merge from Bakta GFF3 Dbxref

Outputs under results_organized_recombinase/genome_clustering/project/notebook7_extended/
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import FF, REPO

SCRIPTS = FF / "02_atlas"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from ie_classifier.classify import classify_from_genes  # noqa: E402
from ie_classifier.config import load_config  # noqa: E402
from ie_classifier.loaders import (  # noqa: E402
    find_att_notes,
    hydrate_genes,
    parse_gbff_cds,
    truth_row_meta,
)
from ie_classifier.schema import TSV_COLUMNS, flatten_record  # noqa: E402

REP_INDEX = REPO / "finder_pipeline/results_recombinase/dedup/rep_index.tsv"
DEDUP_DIR = REPO / "finder_pipeline/results_recombinase/dedup"
BAKTA_ROOT = REPO / "bakta_pipeline/results_bakta_recombinase"
RESULTS_ROOT = REPO / "finder_pipeline/results_recombinase"
ATLAS_CFG = FF / "02_atlas" / "ie_classifier" / "ie_atlas_labels.yaml"
OUT_NB7 = REPO / "results_organized_recombinase/genome_clustering/project/notebook7_extended"
ART = REPO / "final_final/artifacts/recombinase_atlas"


def sanitize_ie_id(ie_id: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", str(ie_id))


def parse_pfam_from_dbxref(raw: str) -> str:
    pfams = []
    for part in str(raw).split(","):
        part = part.strip()
        if part.upper().startswith("PFAM:"):
            pfams.append(part.split(":", 1)[1])
        elif re.match(r"PF\d{5}", part, re.I):
            pfams.append(part)
    return ",".join(pfams)


def build_locus_merge(rep_index: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    rows = []
    for _, row in rep_index.iterrows():
        rep_id = str(row["rep_id"])
        gff3 = BAKTA_ROOT / sanitize_ie_id(row["ie_id"]) / "annotation.gff3"
        if not gff3.is_file():
            continue
        with gff3.open() as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                attrs = parts[8]
                locus = ""
                product = ""
                dbxref = ""
                for tok in attrs.split(";"):
                    tok = tok.strip()
                    if tok.startswith("locus_tag="):
                        locus = tok.split("=", 1)[1]
                    elif tok.startswith("product="):
                        product = tok.split("=", 1)[1]
                    elif tok.startswith("Dbxref="):
                        dbxref = tok.split("=", 1)[1]
                if not locus:
                    continue
                rows.append(
                    {
                        "genome_id": rep_id,
                        "locus_tag": locus,
                        "product": product,
                        "pfam_id": parse_pfam_from_dbxref(dbxref),
                        "emapper_PFAMs": parse_pfam_from_dbxref(dbxref),
                        "eggNOG_Description": "",
                        "seed_eggNOG_ortholog": "",
                        "KEGG_ko": "",
                        "family_id": "",
                        "leiden_hd_cluster": "",
                        "predicted_protein_name": product,
                    }
                )
    df = pd.DataFrame(rows)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"→ {out_tsv} ({len(df):,} locus rows, {df['genome_id'].nunique()} reps)")
    return df


def load_audit_row(sample_id: str, integrase_id: str) -> dict | None:
    audit = RESULTS_ROOT / sample_id / "ie_filter_audit.tsv"
    if not audit.is_file():
        return None
    df = pd.read_csv(audit, sep="\t", dtype=str, low_memory=False)
    sub = df[df["integrase_id"].astype(str) == integrase_id]
    conf = sub[sub["passed_confident"].astype(str).str.lower() == "true"]
    if conf.empty:
        conf = sub
    if conf.empty:
        return None
    return conf.iloc[0].to_dict()


def load_integrase_model(sample_id: str, integrase_id: str) -> str:
    hits = RESULTS_ROOT / sample_id / "integrase_hits_summary.tsv"
    if not hits.is_file():
        return ""
    df = pd.read_csv(hits, sep="\t", dtype=str, low_memory=False)
    id_col = "integrase_id" if "integrase_id" in df.columns else "orf_id"
    sub = df[df[id_col].astype(str) == integrase_id]
    if sub.empty:
        return ""
    row = sub.iloc[0]
    return str(row.get("model", "") or row.get("model_accession", "") or row.get("pfam", ""))


def load_trna_type(sample_id: str, integrase_id: str) -> str:
    trna = RESULTS_ROOT / sample_id / "integrase_trna.tsv"
    if not trna.is_file():
        return ""
    df = pd.read_csv(trna, sep="\t", dtype=str, low_memory=False)
    sub = df[df["integrase_id"].astype(str) == integrase_id].copy()
    if sub.empty:
        return ""
    sub["distance"] = pd.to_numeric(sub["distance"], errors="coerce")
    sub = sub.sort_values("distance")
    return str(sub.iloc[0].get("tRNA_type", ""))


def build_truth_table(rep_index: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    rows = []
    for _, row in rep_index.iterrows():
        rep_id = str(row["rep_id"])
        sample_id = str(row["sample_id"])
        integrase_id = str(row["integrase_id"])
        ie_id = str(row["ie_id"])
        audit = load_audit_row(sample_id, integrase_id) or {}
        rows.append(
            {
                "matrix_row_genome_id": rep_id,
                "sample_id": sample_id,
                "leiden_hd_cluster": "",
                "broad_mge_type": "other",
                "mge_canonical_id": ie_id,
                "integrase_id_lookup": integrase_id,
                "integrase_model": load_integrase_model(sample_id, integrase_id)
                or audit.get("integrase_model", ""),
                "integrase_strand": "",
                "integrase_contig": audit.get("contig", ""),
                "integrase_rows_in_trna_file": "1",
                "tRNA_type_assigned": load_trna_type(sample_id, integrase_id),
                "trna_start": audit.get("trna_start", ""),
                "trna_end": audit.get("trna_end", ""),
                "trna_strand": audit.get("trna_strand", ""),
                "finder_distance_reported": audit.get("best_attl_len_bp", ""),
                "finder_data_missing": "",
                "finder_integrase_trna_path": str(RESULTS_ROOT / sample_id / "integrase_trna.tsv"),
                "orfs_faa_path": "",
                "distance_rule_agrees_min_distance": "",
                "finder_matched_rel": "",
                "nb6_label_policy": "",
                "finder_searched_roots": "results_recombinase",
                "ie_len_nt": audit.get("ie_len_nt", ""),
                "integrase_len_aa": audit.get("integrase_len_aa", ""),
                "rep_file": str(row["file"]),
            }
        )
    df = pd.DataFrame(rows)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"→ {out_tsv} ({len(df):,} rows)")
    return df


def merge_lookup_by_genome(merge_df: pd.DataFrame) -> dict[str, dict[str, dict[str, str]]]:
    out: dict[str, dict[str, dict[str, str]]] = {}
    for _, row in merge_df.iterrows():
        gid = str(row["genome_id"])
        lt = str(row["locus_tag"])
        out.setdefault(gid, {})[lt] = {k: ("" if pd.isna(row[k]) else str(row[k])) for k in row.index}
    return out


def classify_cohort(
    truth: pd.DataFrame,
    merge_df: pd.DataFrame,
    cfg: dict,
    out_tsv: Path,
) -> pd.DataFrame:
    merge_lut = merge_lookup_by_genome(merge_df)
    rows_out = []
    t0 = time.time()
    for i, (_, row) in enumerate(truth.iterrows(), start=1):
        meta = truth_row_meta(row)
        meta["element_id"] = str(row["mge_canonical_id"])
        meta["genome_id"] = str(row["matrix_row_genome_id"])
        meta["sample_id"] = str(row["sample_id"])
        rep_id = str(row["matrix_row_genome_id"])
        bakta_gbff = BAKTA_ROOT / sanitize_ie_id(row["mge_canonical_id"]) / "annotation.gbff"
        finder_gbk = DEDUP_DIR / str(row["rep_file"])
        err = ""
        if not bakta_gbff.is_file():
            err = f"missing_bakta:{bakta_gbff.name}"
            genes = []
            att = {}
        else:
            cds = parse_gbff_cds(bakta_gbff)
            genes = hydrate_genes(cds, merge_lut.get(rep_id, {}))
            att = find_att_notes(finder_gbk) if finder_gbk.is_file() else {}
        if err:
            rec = classify_from_genes([], meta, cfg, att=att)
            rec["error"] = err
            rec["primary_label"] = "other.unclassified_IE"
            rec["coarse_label"] = "unclassified"
            rec["confidence"] = "low"
        else:
            rec = classify_from_genes(
                genes,
                meta,
                cfg,
                att=att,
                rep_gbff=str(finder_gbk.name),
            )
        rows_out.append(flatten_record(rec))
        if i == 1 or i % 100 == 0 or i == len(truth):
            print(f"  classified {i}/{len(truth)} ({time.time()-t0:.1f}s)")

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"→ {out_tsv}")
    return pd.DataFrame(rows_out)


def build_atlas_join(classif: pd.DataFrame, truth: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    cls = classif.copy()
    cls["integrase_id"] = cls["element_id"].astype(str).str.split(":").str[0]
    truth_slim = truth[
        [
            "mge_canonical_id",
            "matrix_row_genome_id",
            "sample_id",
            "integrase_id_lookup",
            "integrase_model",
            "tRNA_type_assigned",
            "ie_len_nt",
            "integrase_len_aa",
        ]
    ].rename(columns={"mge_canonical_id": "element_id", "matrix_row_genome_id": "genome_id"})
    join = cls.merge(
        truth_slim,
        on="element_id",
        how="left",
        suffixes=("", "_truth"),
    )
    gtdb = REPO / "results_organized/misc/gtdb_taxonomy.tsv"
    if gtdb.is_file():
        g = pd.read_csv(gtdb, sep="\t", low_memory=False)
        id_col = "sample_id" if "sample_id" in g.columns else g.columns[0]
        ph_col = next((c for c in g.columns if "phylum" in c.lower()), None)
        if ph_col:
            g = g[[id_col, ph_col]].drop_duplicates(id_col).rename(
                columns={id_col: "sample_id", ph_col: "gtdb_phylum"}
            )
            join = join.merge(g, on="sample_id", how="left")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    join.to_csv(out_tsv, sep="\t", index=False)
    print(f"→ {out_tsv} ({len(join):,} rows)")
    return join


def write_summary(classif: pd.DataFrame, join: pd.DataFrame, out_md: Path) -> None:
    coarse = classif["coarse_label"].value_counts()
    primary = classif["primary_label"].value_counts().head(15)
    novel = classif["novel_tier"].fillna("none").value_counts()
    lines = [
        "# Serine recombinase · IE Atlas classification",
        "",
        f"**Cohort:** {len(classif):,} confident IE reps (PF00239/PF07508 finder)",
        "",
        "## coarse_label",
        "",
        "| coarse_label | n |",
        "|---|---:|",
    ]
    for k, v in coarse.items():
        lines.append(f"| {k} | {v} |")
    lines += ["", "## top primary_label", "", "| primary_label | n |", "|---|---:|"]
    for k, v in primary.items():
        lines.append(f"| {k} | {v} |")
    lines += ["", "## novel_tier", "", "| tier | n |", "|---|---:|"]
    for k, v in novel.items():
        lines.append(f"| {k} | {v} |")
    if "integrase_model" in join.columns:
        lines += ["", "## integrase HMM model (finder)", "", "| model | n |", "|---|---:|"]
        for k, v in join["integrase_model"].value_counts().head(10).items():
            lines.append(f"| {k} | {v} |")
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"→ {out_md}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not REP_INDEX.is_file():
        print(f"ERROR: missing {REP_INDEX}", file=sys.stderr)
        return 1
    if not ATLAS_CFG.is_file():
        print(f"ERROR: missing {ATLAS_CFG}", file=sys.stderr)
        return 1

    OUT_NB7.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    rep_index = pd.read_csv(REP_INDEX, sep="\t")
    cfg = load_config(ATLAS_CFG)

    truth_tsv = OUT_NB7 / "integrase_truth_table_recombinase.tsv"
    merge_tsv = OUT_NB7 / "locus_eggnog_merged.tsv"
    classif_tsv = OUT_NB7 / "ie_instance_classification.tsv"
    join_tsv = ART / "recombinase_atlas_join.tsv"
    summary_md = ART / "REPORT.md"

    if args.force or not truth_tsv.exists():
        truth = build_truth_table(rep_index, truth_tsv)
    else:
        truth = pd.read_csv(truth_tsv, sep="\t", dtype=str, low_memory=False)
        print(f"[SKIP] {truth_tsv}")

    if args.force or not merge_tsv.exists():
        merge_df = build_locus_merge(rep_index, merge_tsv)
    else:
        merge_df = pd.read_csv(merge_tsv, sep="\t", dtype=str, low_memory=False)
        print(f"[SKIP] {merge_tsv}")

    if args.force or not classif_tsv.exists():
        classify_cohort(truth, merge_df, cfg, classif_tsv)
    else:
        print(f"[SKIP] {classif_tsv}")

    classif = pd.read_csv(classif_tsv, sep="\t", dtype=str, low_memory=False)
    join = build_atlas_join(classif, truth, join_tsv)
    write_summary(classif, join, summary_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
