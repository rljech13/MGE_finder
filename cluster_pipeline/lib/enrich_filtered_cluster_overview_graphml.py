#!/usr/bin/env python3
"""Attach per-cluster biology / protein / defense / CRISPR metadata to filtered_ie_cluster_overview*.graphml.

Nodes in the overview graphmls correspond to Leiden-HD clusters (``cluster_id``).
We merge the cluster-level TSVs sitting next to the graphml (produced by the
clustering pipeline) so that Cytoscape can colour / filter clusters by:

* basis for the cluster split — top discriminating Pfam markers
  (``cluster_markers_unclassified_leiden_hd.tsv``);
* protein composition — hallmark roles, top named products, biology summary
  (``cluster_hallmark_composition.tsv``, ``cluster_top_products.tsv``,
  ``cluster_biology_summary.tsv``);
* defence and CRISPR systems — Padloc summaries split into CRISPR-Cas vs
  non-Cas defence (``cluster_padloc_summary.tsv``, ``cluster_padloc_systems.tsv``);
* length / atlas / host summaries.

Each new attribute is stored as a GraphML string (NetworkX requires homogeneous
types per attribute; we serialise floats / ints to canonical decimal strings
and leave ``""`` for missing values).

Usage:
  conda run -n MGE_finder python enrich_filtered_cluster_overview_graphml.py
  conda run -n MGE_finder python enrich_filtered_cluster_overview_graphml.py \
      --cohort unclassified_only
  conda run -n MGE_finder python enrich_filtered_cluster_overview_graphml.py \
      --input /path/to/filtered_ie_cluster_overview_mixing_perm_merged.graphml \
      --metadata-dir /path/to/dir/with/cluster_*.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import networkx as nx
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent

DEFAULT_GRAPHMLS = (
    "filtered_ie_cluster_overview_mixing_perm_merged.graphml",
    "filtered_ie_cluster_overview_perm_filtered.graphml",
    "filtered_ie_cluster_overview_permutation.graphml",
    "filtered_ie_leiden_cluster_overview.graphml",
)

TEXT_CAP = 4000


def _cap(val: object) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return ""
    if len(s) > TEXT_CAP:
        return s[: TEXT_CAP - 1] + "\u2026"
    return s


def _float_str(val: object, ndigits: int | None = None) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return ""
    try:
        f = float(s)
    except (TypeError, ValueError):
        return _cap(val)
    if ndigits is not None:
        return f"{round(f, ndigits):.{ndigits}f}"
    return repr(f)


def _int_str(val: object) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return ""
    try:
        return str(int(float(s)))
    except (TypeError, ValueError):
        return _cap(val)


def _norm_cluster(val: object) -> str | None:
    """Cluster ids may be int / float / str — coerce to canonical str like "0"."""
    if val is None:
        return None
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return None
    try:
        return str(int(float(s)))
    except (TypeError, ValueError):
        return s


def _safe_attr_name(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    if not s:
        s = "x"
    if s[0].isdigit():
        s = "_" + s
    return s


# -----------------------------------------------------------------------------
# Loaders for the cluster_*.tsv tables
# -----------------------------------------------------------------------------

BIOLOGY_TEXT_COLS = (
    "group_composition",
    "top_host_phylum",
    "top_named_products",
)
BIOLOGY_INT_COLS = (
    "n_ie",
    "ie_len_median",
    "ie_len_q25",
    "ie_len_q75",
    "n_distinct_host_phyla",
)
BIOLOGY_FLOAT_COLS = (
    "top_host_phylum_frac",
    "median_IE_core_score",
    "median_replication_score",
    "frac_integrase_tyrosine",
    "frac_integrase_serine",
    "frac_transposase",
    "frac_replication",
    "frac_methyltransferase",
    "frac_restriction_T1",
    "frac_defense_STAND",
    "frac_toxin_antitoxin",
    "frac_regulator_HTH",
    "frac_aaa_atpase",
    "frac_tpr_repeat",
    "frac_helicase",
    "frac_zn_ribbon",
)


def load_biology(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = _norm_cluster(row.get("cluster"))
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for col in BIOLOGY_TEXT_COLS:
            rec[f"bio_{col}"] = _cap(row.get(col))
        for col in BIOLOGY_INT_COLS:
            rec[f"bio_{col}"] = _int_str(row.get(col))
        for col in BIOLOGY_FLOAT_COLS:
            rec[f"bio_{col}"] = _float_str(row.get(col), ndigits=4)
        out[cid] = rec
    return out


PADLOC_SUMMARY_TEXT_COLS = ("top_systems",)
PADLOC_SUMMARY_INT_COLS = ("n_ie", "n_ie_with_padloc", "n_distinct_systems")
PADLOC_SUMMARY_FLOAT_COLS = ("frac_ie_with_padloc",)


def load_padloc_summary(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = _norm_cluster(row.get("cluster"))
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for col in PADLOC_SUMMARY_TEXT_COLS:
            rec[f"padloc_{col}"] = _cap(row.get(col))
        for col in PADLOC_SUMMARY_INT_COLS:
            rec[f"padloc_{col}"] = _int_str(row.get(col))
        for col in PADLOC_SUMMARY_FLOAT_COLS:
            rec[f"padloc_{col}"] = _float_str(row.get(col), ndigits=4)
        out[cid] = rec
    return out


CRISPR_RE = re.compile(r"^(cas_|crispr)", re.IGNORECASE)


def _is_crispr(system: str) -> bool:
    return bool(CRISPR_RE.match(str(system or "")))


def _fmt_system_token(name: str, n_ie: str | float, frac: str | float) -> str:
    n_ie_s = _int_str(n_ie)
    frac_s = _float_str(frac, ndigits=3)
    if n_ie_s == "" and frac_s == "":
        return str(name)
    parts = [str(name)]
    extras = []
    if n_ie_s:
        extras.append(f"{n_ie_s} IE")
    if frac_s:
        extras.append(frac_s)
    return f"{parts[0]} ({', '.join(extras)})"


def load_padloc_systems(path: Path, top_n: int = 10) -> dict[str, dict[str, str]]:
    """Split per-cluster Padloc hits into CRISPR-Cas vs non-Cas defence."""
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df["system"] = df["system"].astype(str)
    df["nie_num"] = pd.to_numeric(df["n_ie_with_system"], errors="coerce").fillna(0.0)
    df["frac_num"] = pd.to_numeric(df["frac_ie_with_system"], errors="coerce").fillna(0.0)
    df["hits_num"] = pd.to_numeric(df.get("n_total_hits"), errors="coerce").fillna(0.0)

    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby("cluster", sort=False):
        cid = _norm_cluster(cid_raw)
        if cid is None:
            continue
        crispr = sub[sub["system"].map(_is_crispr)].sort_values(
            ["nie_num", "frac_num"], ascending=[False, False]
        )
        defence = sub[~sub["system"].map(_is_crispr)].sort_values(
            ["nie_num", "frac_num"], ascending=[False, False]
        )

        rec: dict[str, str] = {}
        rec["padloc_n_distinct_crispr_systems"] = str(len(crispr))
        rec["padloc_n_distinct_defense_systems"] = str(len(defence))
        rec["padloc_crispr_n_ie_total"] = _int_str(crispr["nie_num"].sum())
        rec["padloc_defense_n_ie_total"] = _int_str(defence["nie_num"].sum())
        rec["padloc_crispr_n_hits_total"] = _int_str(crispr["hits_num"].sum())
        rec["padloc_defense_n_hits_total"] = _int_str(defence["hits_num"].sum())

        rec["padloc_top_crispr_systems"] = "; ".join(
            _fmt_system_token(r.system, r.nie_num, r.frac_num)
            for r in crispr.head(top_n).itertuples()
        )
        rec["padloc_top_defense_systems"] = "; ".join(
            _fmt_system_token(r.system, r.nie_num, r.frac_num)
            for r in defence.head(top_n).itertuples()
        )
        rec["padloc_crispr_systems_list"] = "; ".join(crispr["system"].tolist())
        out[cid] = rec
    return out


def load_hallmark(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Pivot cluster_hallmark_composition.tsv into per-role columns.

    For each role we emit two attributes:
      - hallmark_<role>_frac          (fraction of IE in the cluster with the role)
      - hallmark_<role>_mean_hits     (mean pfam hits/IE for that role)

    Plus a single summary attribute ``hallmark_top_roles`` with the top 5 roles
    by fraction (excluding ``other``).
    """
    if not path.is_file():
        return {}, []
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df["frac_num"] = pd.to_numeric(df.get("frac_ie_with_role"), errors="coerce").fillna(0.0)
    df["mean_num"] = pd.to_numeric(df.get("mean_pfam_hits_per_ie"), errors="coerce").fillna(0.0)
    df["role"] = df["role"].astype(str)

    roles = sorted(df["role"].dropna().unique().tolist())
    role_attrs: list[str] = []
    out: dict[str, dict[str, str]] = {}

    for cid_raw, sub in df.groupby("cluster", sort=False):
        cid = _norm_cluster(cid_raw)
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for _, row in sub.iterrows():
            role = _safe_attr_name(row["role"])
            rec[f"hallmark_{role}_frac"] = _float_str(row["frac_num"], ndigits=4)
            rec[f"hallmark_{role}_mean_hits"] = _float_str(row["mean_num"], ndigits=3)
        non_other = sub[sub["role"].str.lower() != "other"].sort_values(
            "frac_num", ascending=False
        )
        top_tokens = []
        for r in non_other.head(5).itertuples():
            label_raw = getattr(r, "role_label", None)
            if label_raw is None or pd.isna(label_raw) or str(label_raw).strip() == "":
                label = r.role
            else:
                label = str(label_raw)
            top_tokens.append(f"{label} ({r.frac_num:.2f})")
        rec["hallmark_top_roles"] = "; ".join(top_tokens)
        out[cid] = rec

    for role in roles:
        safe = _safe_attr_name(role)
        role_attrs.append(f"hallmark_{safe}_frac")
        role_attrs.append(f"hallmark_{safe}_mean_hits")
    role_attrs.append("hallmark_top_roles")
    return out, role_attrs


def load_top_products(path: Path, top_n: int = 10) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cluster_col = "leiden_hd_cluster" if "leiden_hd_cluster" in df.columns else "cluster"
    df["per_ie_num"] = pd.to_numeric(df.get("per_ie"), errors="coerce").fillna(0.0)
    df["nprot_num"] = pd.to_numeric(df.get("n_proteins"), errors="coerce").fillna(0.0)
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby(cluster_col, sort=False):
        cid = _norm_cluster(cid_raw)
        if cid is None:
            continue
        sub = sub.sort_values(["per_ie_num", "nprot_num"], ascending=[False, False]).head(top_n)
        tokens = []
        for r in sub.itertuples():
            prod_raw = getattr(r, "product_norm", "")
            prod = "" if pd.isna(prod_raw) else str(prod_raw)
            if not prod:
                continue
            tokens.append(
                f"{prod} ({r.per_ie_num:.2f}/IE, n={int(r.nprot_num)})"
            )
        out[cid] = {"top_products": _cap("; ".join(tokens))}
    return out


def load_markers(path: Path, top_n: int = 10) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df["score_num"] = pd.to_numeric(df.get("marker_score"), errors="coerce").fillna(0.0)
    df["enrich_num"] = pd.to_numeric(df.get("enrichment"), errors="coerce").fillna(0.0)
    df["pin_num"] = pd.to_numeric(df.get("p_in"), errors="coerce").fillna(0.0)
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby("cluster", sort=False):
        cid = _norm_cluster(cid_raw)
        if cid is None:
            continue
        sub = sub.sort_values("score_num", ascending=False).head(top_n)
        tokens = []
        roles_seen: list[str] = []
        for r in sub.itertuples():
            short_raw = getattr(r, "short_name", "")
            short = "" if pd.isna(short_raw) else str(short_raw)
            feat = getattr(r, "feature_id", "")
            name = short.strip() if short.strip() else str(feat)
            tokens.append(
                f"{name} ({r.enrich_num:.1f}x, p_in={r.pin_num:.2f}, score={r.score_num:.2f})"
            )
            role_raw = getattr(r, "role", "")
            role = "" if pd.isna(role_raw) else str(role_raw)
            if role and role.lower() != "other" and role not in roles_seen:
                roles_seen.append(role)
        rec = {
            "markers_top": _cap("; ".join(tokens)),
            "markers_top_roles": "; ".join(roles_seen[:5]),
            "markers_n_features": str(len(sub)),
        }
        out[cid] = rec
    return out


IE_LEN_INT_COLS = ("n",)
IE_LEN_FLOAT_COLS = ("median", "mean", "q25", "q75")


def load_ie_length(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = _norm_cluster(row.get("cluster"))
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for col in IE_LEN_INT_COLS:
            rec[f"ie_len_{col}"] = _int_str(row.get(col))
        for col in IE_LEN_FLOAT_COLS:
            rec[f"ie_len_{col}"] = _float_str(row.get(col), ndigits=1)
        out[cid] = rec
    return out


ATLAS_FLOAT_COLS = (
    "IE_core_score",
    "phage_score",
    "inovirus_score",
    "conjugation_score",
    "replication_score",
)


def load_atlas_medians(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = _norm_cluster(row.get("cluster"))
        if cid is None:
            continue
        rec: dict[str, str] = {}
        rec["atlas_median_n_ie"] = _int_str(row.get("n_ie"))
        for col in ATLAS_FLOAT_COLS:
            rec[f"atlas_median_{col}"] = _float_str(row.get(col), ndigits=4)
        out[cid] = rec
    return out


ARCHETYPE_COLS = ("archetype_short", "archetype_family", "archetype_color", "archetype_long")


def load_archetype_labels(path: Path) -> dict[str, dict[str, str]]:
    """Hand-curated biological label per Leiden cluster (cluster_archetype_labels.tsv)."""
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = _norm_cluster(row.get("cluster"))
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for col in ARCHETYPE_COLS:
            rec[col] = _cap(row.get(col))
        out[cid] = rec
    return out


# -----------------------------------------------------------------------------
# Enrichment driver
# -----------------------------------------------------------------------------

def _ensure_attr(graph: nx.Graph, attr: str) -> None:
    """Make sure every node carries the attribute (string default '')."""
    for node in graph.nodes():
        if attr not in graph.nodes[node]:
            graph.nodes[node][attr] = ""


def _collect_attr_names(lookups: Iterable[dict[str, dict[str, str]]]) -> list[str]:
    attrs: set[str] = set()
    for lut in lookups:
        for rec in lut.values():
            attrs.update(rec.keys())
    return sorted(attrs)


def enrich_graph(
    graph: nx.Graph,
    bio: dict[str, dict[str, str]],
    padloc_summary: dict[str, dict[str, str]],
    padloc_systems: dict[str, dict[str, str]],
    hallmark: dict[str, dict[str, str]],
    hallmark_attrs: list[str],
    products: dict[str, dict[str, str]],
    markers: dict[str, dict[str, str]],
    ie_len: dict[str, dict[str, str]],
    atlas: dict[str, dict[str, str]],
    archetype: dict[str, dict[str, str]],
) -> dict[str, int]:
    luts = (bio, padloc_summary, padloc_systems, hallmark, products, markers,
            ie_len, atlas, archetype)
    attr_names = set(_collect_attr_names(luts)) | set(hallmark_attrs) | set(ARCHETYPE_COLS)

    for attr in attr_names:
        _ensure_attr(graph, attr)

    stats = {
        "nodes": 0,
        "bio": 0,
        "padloc_summary": 0,
        "padloc_systems": 0,
        "hallmark": 0,
        "products": 0,
        "markers": 0,
        "ie_len": 0,
        "atlas": 0,
        "archetype": 0,
    }

    for node, attr in list(graph.nodes(data=True)):
        stats["nodes"] += 1
        cid_raw = attr.get("cluster_id", attr.get("cluster_id_str", ""))
        cid = _norm_cluster(cid_raw)
        if cid is None:
            continue
        for key, lut in (
            ("bio", bio),
            ("padloc_summary", padloc_summary),
            ("padloc_systems", padloc_systems),
            ("hallmark", hallmark),
            ("products", products),
            ("markers", markers),
            ("ie_len", ie_len),
            ("atlas", atlas),
            ("archetype", archetype),
        ):
            rec = lut.get(cid)
            if not rec:
                continue
            stats[key] += 1
            for k, v in rec.items():
                graph.nodes[node][k] = v
    return stats


def resolve_metadata_dir(graphml_path: Path, override: Path | None) -> Path:
    if override is not None:
        return override
    return graphml_path.parent


def enrich_one(graphml_path: Path, meta_dir: Path, out_path: Path) -> None:
    graph = nx.read_graphml(graphml_path)

    bio = load_biology(meta_dir / "cluster_biology_summary.tsv")
    padloc_summary = load_padloc_summary(meta_dir / "cluster_padloc_summary.tsv")
    padloc_systems = load_padloc_systems(meta_dir / "cluster_padloc_systems.tsv")
    hallmark, hallmark_attrs = load_hallmark(meta_dir / "cluster_hallmark_composition.tsv")
    products = load_top_products(meta_dir / "cluster_top_products.tsv")
    markers_path = meta_dir / "cluster_markers_unclassified_leiden_hd.tsv"
    if not markers_path.is_file():
        for cand in meta_dir.glob("cluster_markers_*.tsv"):
            markers_path = cand
            break
    markers = load_markers(markers_path)
    ie_len = load_ie_length(meta_dir / "cluster_ie_length_stats.tsv")
    atlas = load_atlas_medians(meta_dir / "cluster_atlas_score_medians.tsv")
    archetype = load_archetype_labels(meta_dir / "cluster_archetype_labels.tsv")

    stats = enrich_graph(
        graph,
        bio=bio,
        padloc_summary=padloc_summary,
        padloc_systems=padloc_systems,
        hallmark=hallmark,
        hallmark_attrs=hallmark_attrs,
        products=products,
        markers=markers,
        ie_len=ie_len,
        atlas=atlas,
        archetype=archetype,
    )

    existing = graph.graph.get("description", "")
    note = (
        "Per-cluster enrichment: cluster_biology_summary.tsv, "
        "cluster_padloc_{summary,systems}.tsv (CRISPR=cas_*/CRISPR*), "
        "cluster_hallmark_composition.tsv, cluster_top_products.tsv, "
        "cluster_markers_*.tsv, cluster_ie_length_stats.tsv, "
        "cluster_atlas_score_medians.tsv, "
        "cluster_archetype_labels.tsv (hand-curated biological grouping: "
        "archetype_short/family/color/long)."
    )
    graph.graph["cluster_enrichment_note"] = note
    if existing and "cluster_enrichment_note" not in existing:
        graph.graph["description"] = f"{existing} | {note}"
    elif not existing:
        graph.graph["description"] = note

    out_path.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(graph, str(out_path))

    print(f"[enrich] {graphml_path.name} -> {out_path}")
    print(
        f"  nodes={stats['nodes']} "
        f"bio={stats['bio']} padloc_summary={stats['padloc_summary']} "
        f"padloc_systems={stats['padloc_systems']} hallmark={stats['hallmark']} "
        f"products={stats['products']} markers={stats['markers']} "
        f"ie_len={stats['ie_len']} atlas={stats['atlas']} "
        f"archetype={stats['archetype']}"
    )


def iter_default_inputs(cohort_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for name in DEFAULT_GRAPHMLS:
        p = cohort_dir / name
        if p.is_file():
            paths.append(p)
    return paths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cohort",
        choices=("unclassified_only", "all"),
        default="unclassified_only",
        help="Which cohort folder under final/clustering_filtered_ie/ to enrich.",
    )
    ap.add_argument(
        "--input",
        type=Path,
        action="append",
        default=None,
        help="Override input graphml (repeatable). Defaults: filtered_ie_cluster_overview_*.graphml",
    )
    ap.add_argument(
        "--metadata-dir",
        type=Path,
        default=None,
        help="Override directory containing cluster_*.tsv (defaults to graphml's directory).",
    )
    ap.add_argument(
        "--out-suffix",
        default="_enriched",
        help="Output suffix appended before .graphml (default: _enriched).",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Write enriched graphml in-place (suffix is ignored).",
    )
    args = ap.parse_args()

    if args.cohort == "all":
        cohort_dir = BASE_DIR
    else:
        cohort_dir = BASE_DIR / "unclassified_only"

    inputs = list(args.input) if args.input else iter_default_inputs(cohort_dir)
    if not inputs:
        raise SystemExit(f"No graphml inputs found under {cohort_dir}")

    for graphml_path in inputs:
        graphml_path = Path(graphml_path).resolve()
        if not graphml_path.is_file():
            print(f"[skip] missing: {graphml_path}")
            continue
        meta_dir = resolve_metadata_dir(graphml_path, args.metadata_dir)
        if args.overwrite:
            out_path = graphml_path
        else:
            out_path = graphml_path.with_name(
                graphml_path.stem + args.out_suffix + graphml_path.suffix
            )
        enrich_one(graphml_path, meta_dir, out_path)


if __name__ == "__main__":
    main()
