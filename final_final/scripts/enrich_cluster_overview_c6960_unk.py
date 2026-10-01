#!/usr/bin/env python3
"""Attach full per-cluster biology / marker / defense metadata to Stage 9 unc graphmls.

Mirrors ``final/clustering_filtered_ie/enrich_filtered_cluster_overview_graphml.py``
for the Stage 9 unclassified Leiden rebuild (14 clusters, n=2,589 IE).

Adds Cytoscape-friendly marker colouring fields (compatible with the old archetype_* names):
  marker_label_short, marker_family, marker_color, marker_label_long
  archetype_short / archetype_family / archetype_color / archetype_long  (aliases)

Usage:
  conda run -n MGE_finder python enrich_cluster_overview_c6960_unk.py
"""
from __future__ import annotations

import math
import re
import sys
from collections import Counter
from pathlib import Path

import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import REPO
OLD_ENRICH = REPO / "final_final" / "03_clustering" / "lib" / "enrich_filtered_cluster_overview_graphml.py"
sys.path.insert(0, str(OLD_ENRICH.parent))
import enrich_filtered_cluster_overview_graphml as ef  # noqa: E402

ART = REPO / "final_final/artifacts/nb04_cluster_unknown_6960"

GRAPHMLS = (
    "cluster_overview_mixing_null.graphml",
    "cluster_overview_permutation.graphml",
    "cluster_overview_perm_filtered.graphml",
    "cluster_overview_mixing_perm_merged.graphml",
    "cluster_network_cytoscape.graphml",
    "cluster_network_strong.graphml",
)

SKIP_HALLMARK_FOR_COLOR = {"integrase_tyrosine", "regulator_HTH", "other"}

MARKER_FAMILY_PALETTE = {
    "HTH_regulatory": "#d62728",
    "replication": "#ff7f0e",
    "integrase_SAM": "#1f77b4",
    "integrase_serine": "#e377c2",
    "defense_RM": "#2ca02c",
    "toxin_antitoxin": "#9467bd",
    "transposase": "#8c564b",
    "AAA_motor": "#bcbd22",
    "TPR_scaffold": "#17becf",
    "coiled_coil": "#aec7e8",
    "zn_ribbon": "#ff9896",
    "minimal_integrase": "#5fa9d6",
    "other": "#7f7f7f",
}

# Human-readable Cytoscape labels (variant 1: colour = marker_family, label = C# · family).
MARKER_FAMILY_HUMAN = {
    "HTH_regulatory": "HTH regulators",
    "replication": "replication",
    "integrase_SAM": "integrase (SAM)",
    "integrase_serine": "serine integrase",
    "minimal_integrase": "phage integrase",
    "defense_RM": "RM / methyltransferase",
    "toxin_antitoxin": "toxin–antitoxin",
    "transposase": "transposase",
    "AAA_motor": "AAA-ATPase",
    "TPR_scaffold": "TPR repeats",
    "coiled_coil": "coiled-coil",
    "zn_ribbon": "zinc ribbon",
    "other": "other",
}


def _cluster_col(df: pd.DataFrame) -> str:
    for c in ("leiden_cluster", "cluster", "cluster_id"):
        if c in df.columns:
            return c
    raise KeyError(f"No cluster column in {list(df.columns)}")


def _infer_marker_family(short_name: str, top_role: str) -> str:
    s = (short_name or "").upper()
    role = (top_role or "").lower()
    if any(x in s for x in ("HTH", "MER", "XRE", "CRO", "TRMB")):
        return "HTH_regulatory"
    if any(x in s for x in ("RM_", "MTASE", "METHYL", "RESTRICTION", "HSD", "N6_")):
        return "defense_RM"
    if any(x in s for x in ("REP", "HELIC", "PRI", "DNAB", "REPLITRON")):
        return "replication"
    if any(x in s for x in ("PHAGE_INT", "INT_N", "INT_SAM", "INTEGRASE")):
        return "integrase_SAM" if "SAM" in s else "minimal_integrase"
    if "TPR" in s:
        return "TPR_scaffold"
    if "AAA" in s:
        return "AAA_motor"
    if "ZN" in s and "RIBBON" in s:
        return "zn_ribbon"
    if "CCDC" in s or "COIL" in s:
        return "coiled_coil"
    if any(x in s for x in ("TRANSP", "IS_", "INS")):
        return "transposase"
    if role == "integrase_serine" or "INT_M" in s or "RESOLVASE" in s:
        return "integrase_serine"
    if role == "toxin_antitoxin" or any(x in s for x in ("ABI", "ANTITOX", "TOXIN")):
        return "toxin_antitoxin"
    if role == "replication":
        return "replication"
    if role == "tpr_repeat":
        return "TPR_scaffold"
    if role == "aaa_atpase":
        return "AAA_motor"
    if role == "zn_ribbon":
        return "zn_ribbon"
    if role == "coiled_coil":
        return "coiled_coil"
    if role == "methyltransferase":
        return "defense_RM"
    if role == "transposase":
        return "transposase"
    return "other"


def load_stage8_biology(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = ef._norm_cluster(row.get(cc))
        if cid is None:
            continue
        rec = {
            "bio_n_ie": ef._int_str(row.get("n_IE")),
            "bio_ie_len_median": ef._int_str(row.get("ie_len_med")),
            "bio_ie_len_q25": ef._int_str(row.get("ie_len_q25")),
            "bio_ie_len_q75": ef._int_str(row.get("ie_len_q75")),
            "bio_integrase_aa_med": ef._int_str(row.get("integrase_aa_med")),
            "bio_top_host_phylum": ef._cap(row.get("top_phylum")),
            "bio_top_host_phylum_frac": ef._float_str(row.get("top_phylum_frac"), ndigits=4),
            "bio_n_distinct_host_phyla": ef._int_str(row.get("n_phyla")),
            "bio_frac_high_interest_novel": ef._float_str(row.get("frac_high_interest_novel"), ndigits=4),
            "bio_median_IE_core_score": ef._float_str(row.get("IE_core_med"), ndigits=4),
        }
        out[cid] = rec
    return out


def load_stage8_padloc_summary(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = ef._norm_cluster(row.get(cc))
        if cid is None:
            continue
        out[cid] = {
            "padloc_n_ie": ef._int_str(row.get("n_ie")),
            "padloc_n_ie_with_padloc": ef._int_str(row.get("n_with_padloc")),
            "padloc_frac_ie_with_padloc": ef._float_str(row.get("frac_with_padloc"), ndigits=4),
            "padloc_mean_systems": ef._float_str(row.get("mean_systems"), ndigits=3),
        }
    return out


def build_padloc_systems_from_per_ie(path: Path, top_n: int = 10) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cc = _cluster_col(df)
    rows: list[dict] = []
    for _, r in df.iterrows():
        cid = ef._norm_cluster(r.get(cc))
        if cid is None:
            continue
        ie = str(r.get("ie_id", r.get("integrase_id", "")))
        for sys_name in str(r.get("padloc_systems", "")).split("|"):
            sys_name = sys_name.strip()
            if sys_name:
                rows.append({"cluster": cid, "system": sys_name, "ie_id": ie})
    if not rows:
        return {}
    sys_df = pd.DataFrame(rows)
    cluster_sizes = df.groupby(cc).size()
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in sys_df.groupby("cluster", sort=False):
        cid = ef._norm_cluster(cid_raw)
        if cid is None:
            continue
        n_cluster = int(cluster_sizes.get(int(cid) if cid.isdigit() else cid_raw, len(sub)))
        agg = (
            sub.groupby("system")
            .agg(n_ie_with_system=("ie_id", "nunique"), n_total_hits=("ie_id", "count"))
            .reset_index()
        )
        agg["frac_ie_with_system"] = agg["n_ie_with_system"] / max(1, n_cluster)
        agg = agg.sort_values(["n_ie_with_system", "frac_ie_with_system"], ascending=False)
        crispr = agg[agg["system"].map(ef._is_crispr)]
        defence = agg[~agg["system"].map(ef._is_crispr)]
        tokens = [
            ef._fmt_system_token(r.system, r.n_ie_with_system, r.frac_ie_with_system)
            for r in agg.head(top_n).itertuples()
        ]
        out[cid] = {
            "padloc_top_systems": ef._cap("; ".join(tokens)),
            "padloc_top_crispr_systems": "; ".join(
                ef._fmt_system_token(r.system, r.n_ie_with_system, r.frac_ie_with_system)
                for r in crispr.head(top_n).itertuples()
            ),
            "padloc_top_defense_systems": "; ".join(
                ef._fmt_system_token(r.system, r.n_ie_with_system, r.frac_ie_with_system)
                for r in defence.head(top_n).itertuples()
            ),
            "padloc_n_distinct_crispr_systems": str(len(crispr)),
            "padloc_n_distinct_defense_systems": str(len(defence)),
            "padloc_crispr_n_ie_total": ef._int_str(crispr["n_ie_with_system"].sum()),
            "padloc_defense_n_ie_total": ef._int_str(defence["n_ie_with_system"].sum()),
        }
    return out


def load_stage8_hallmark(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    if not path.is_file():
        return {}, []
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    df["frac_num"] = pd.to_numeric(df.get("frac_ie_with_role"), errors="coerce").fillna(0.0)
    df["mean_num"] = pd.to_numeric(df.get("mean_pfam_hits_per_ie"), errors="coerce").fillna(0.0)
    df["role"] = df["role"].astype(str)
    roles = sorted(df["role"].dropna().unique().tolist())
    role_attrs: list[str] = []
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby(cc, sort=False):
        cid = ef._norm_cluster(cid_raw)
        if cid is None:
            continue
        rec: dict[str, str] = {}
        for _, row in sub.iterrows():
            role = ef._safe_attr_name(row["role"])
            rec[f"hallmark_{role}_frac"] = ef._float_str(row["frac_num"], ndigits=4)
            rec[f"hallmark_{role}_mean_hits"] = ef._float_str(row["mean_num"], ndigits=3)
        non_other = sub[sub["role"].str.lower() != "other"].sort_values("frac_num", ascending=False)
        top_tokens = [f"{r.role} ({r.frac_num:.2f})" for r in non_other.head(5).itertuples()]
        rec["hallmark_top_roles"] = "; ".join(top_tokens)
        out[cid] = rec
    for role in roles:
        safe = ef._safe_attr_name(role)
        role_attrs.extend([f"hallmark_{safe}_frac", f"hallmark_{safe}_mean_hits"])
    role_attrs.append("hallmark_top_roles")
    return out, role_attrs


def load_stage8_markers(path: Path, top_n: int = 10) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    df["score_num"] = pd.to_numeric(df.get("marker_score"), errors="coerce").fillna(0.0)
    df["enrich_num"] = pd.to_numeric(df.get("enrichment"), errors="coerce").fillna(0.0)
    df["pin_num"] = pd.to_numeric(df.get("p_in"), errors="coerce").fillna(0.0)
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby(cc, sort=False):
        cid = ef._norm_cluster(cid_raw)
        if cid is None:
            continue
        sub = sub.sort_values("score_num", ascending=False).head(top_n)
        tokens = []
        for r in sub.itertuples():
            short = "" if pd.isna(getattr(r, "short_name", "")) else str(getattr(r, "short_name", "")).strip()
            feat = getattr(r, "feature_id", "")
            pfam = getattr(r, "pfam_acc", "")
            name = short if short else str(feat)
            tokens.append(
                f"{name} ({r.enrich_num:.1f}x, p_in={r.pin_num:.2f}, score={r.score_num:.2f})"
            )
        top = sub.iloc[0]
        short1 = str(top.get("short_name", "") or "").strip()
        if not short1 or short1.lower() == "nan":
            pfam = str(top.get("pfam_acc", "") or top.get("feature_id", "")).strip()
            short1 = pfam.replace("pfam_", "").split(".")[0] if pfam else ""
        rec = {
            "markers_top": ef._cap("; ".join(tokens)),
            "markers_n_features": str(len(sub)),
            "markers_top1_short": ef._cap(short1),
            "markers_top1_pfam": ef._cap(top.get("pfam_acc", "")),
            "markers_top1_full": ef._cap(top.get("full_name", "")),
            "markers_top1_score": ef._float_str(top.get("marker_score"), ndigits=3),
            "markers_top1_enrichment": ef._float_str(top.get("enrichment"), ndigits=2),
            "markers_top1_p_in": ef._float_str(top.get("p_in"), ndigits=3),
        }
        out[cid] = rec
    return out


def load_stage8_host_purity(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = ef._norm_cluster(row.get(cc))
        if cid is None:
            continue
        out[cid] = {
            "host_purity_phylum": ef._float_str(row.get("purity_phylum"), ndigits=4),
            "host_shannon_phylum_bits": ef._float_str(row.get("shannon_phylum_bits"), ndigits=3),
            "host_effective_phyla": ef._float_str(row.get("effective_phyla"), ndigits=2),
            "host_n_with_phylum": ef._int_str(row.get("n_with_phylum")),
        }
    return out


def load_stage8_diversity(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = ef._norm_cluster(row.get(cc))
        if cid is None:
            continue
        out[cid] = {
            "diversity_shannon_phylum": ef._float_str(row.get("shannon_phylum"), ndigits=3),
            "diversity_shannon_genus": ef._float_str(row.get("shannon_genus"), ndigits=3),
            "diversity_eff_phyla": ef._float_str(row.get("eff_phyla"), ndigits=2),
            "diversity_eff_genera": ef._float_str(row.get("eff_genera"), ndigits=1),
            "diversity_top_phylum_frac": ef._float_str(row.get("top_phylum_frac"), ndigits=4),
            "diversity_top_genus_frac": ef._float_str(row.get("top_genus_frac"), ndigits=4),
        }
    return out


def load_stage8_integrase_div(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    cc = _cluster_col(df)
    out: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cid = ef._norm_cluster(row.get(cc))
        if cid is None:
            continue
        out[cid] = {
            "integrase_frac_PF00589": ef._float_str(row.get("frac_PF00589"), ndigits=3),
            "integrase_frac_PF22022": ef._float_str(row.get("frac_PF22022"), ndigits=3),
            "integrase_frac_serine_pfam": ef._float_str(row.get("frac_serine_pfam"), ndigits=3),
            "integrase_n_tyr_only": ef._int_str(row.get("n_tyr_only")),
            "integrase_n_tyr_plus_ser": ef._int_str(row.get("n_tyr_plus_ser")),
        }
    return out


def load_stage8_atlas_medians(
    rebuilt_path: Path,
    atlas_join_path: Path | None = None,
    leiden_path: Path | None = None,
) -> dict[str, dict[str, str]]:
    score_cols = ["IE_core_score", "phage_score", "inovirus_score",
                  "conjugation_score", "replication_score"]
    df = None
    if atlas_join_path and leiden_path and atlas_join_path.is_file() and leiden_path.is_file():
        atlas = pd.read_csv(atlas_join_path, sep="\t", low_memory=False)
        leid = pd.read_csv(leiden_path, sep="\t")
        atlas["integrase_id"] = atlas["element_id"].astype(str).str.split(":").str[0]
        cols = ["integrase_id"] + [c for c in score_cols if c in atlas.columns]
        df = leid.merge(atlas[cols], on="integrase_id", how="inner")
        cc = "leiden_cluster"
    elif rebuilt_path.is_file():
        df = pd.read_csv(rebuilt_path, sep="\t", low_memory=False)
        cc = "leiden_cluster"
        score_cols = [c for c in score_cols if c in df.columns]
    else:
        return {}
    if df is None or cc not in df.columns or not score_cols:
        return {}
    out: dict[str, dict[str, str]] = {}
    for cid_raw, sub in df.groupby(cc, sort=False):
        cid = ef._norm_cluster(cid_raw)
        if cid is None:
            continue
        rec = {"atlas_median_n_ie": ef._int_str(len(sub))}
        for col in score_cols:
            rec[f"atlas_median_{col}"] = ef._float_str(sub[col].median(), ndigits=4)
        out[cid] = rec
    return out


def _hallmark_bio_frac_aliases(hallmark: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    """Map hallmark role fractions to legacy bio_frac_* attribute names."""
    role_map = {
        "integrase_tyrosine": "bio_frac_integrase_tyrosine",
        "integrase_serine": "bio_frac_integrase_serine",
        "transposase": "bio_frac_transposase",
        "replication": "bio_frac_replication",
        "helicase": "bio_frac_helicase",
        "aaa_atpase": "bio_frac_aaa_atpase",
        "methyltransferase": "bio_frac_methyltransferase",
        "restriction_T1": "bio_frac_restriction_T1",
        "defense_STAND": "bio_frac_defense_STAND",
        "toxin_antitoxin": "bio_frac_toxin_antitoxin",
        "regulator_HTH": "bio_frac_regulator_HTH",
        "tpr_repeat": "bio_frac_tpr_repeat",
        "zn_ribbon": "bio_frac_zn_ribbon",
    }
    out: dict[str, dict[str, str]] = {}
    for cid, rec in hallmark.items():
        aliases: dict[str, str] = {}
        for role, bio_key in role_map.items():
            safe = ef._safe_attr_name(role)
            val = rec.get(f"hallmark_{safe}_frac", "")
            if val != "":
                aliases[bio_key] = val
        if aliases:
            out[cid] = aliases
    return out


def _top_hallmark_role(hallmark: dict[str, dict[str, str]], cid: str) -> str:
    rec = hallmark.get(cid, {})
    best_role, best_frac = "", -1.0
    for k, v in rec.items():
        if not k.startswith("hallmark_") or not k.endswith("_frac"):
            continue
        role = k[len("hallmark_") : -len("_frac")]
        if role in SKIP_HALLMARK_FOR_COLOR:
            continue
        try:
            frac = float(v)
        except (TypeError, ValueError):
            continue
        if frac > best_frac:
            best_frac, best_role = frac, role
    return best_role


def _cluster_subtype_human(
    family: str,
    mrec: dict[str, str],
    bio_rec: dict[str, str],
    host_rec: dict[str, str],
    intdiv_rec: dict[str, str],
) -> str:
    """Short tag explaining why clusters with the same marker_family still differ."""
    markers_top = mrec.get("markers_top", "")
    short1 = mrec.get("markers_top1_short", "")
    tags: list[str] = []

    if family == "HTH_regulatory":
        if "MarR" in markers_top:
            tags.append("MarR")
        elif "MerR" in markers_top:
            tags.append("MerR")
        if any(x in markers_top for x in ("HTH_3", "HTH_26", "HTH_19", "HTH_31", "HTH_17")):
            if "MarR" not in tags and "MerR" not in tags:
                tags.append("Cro-type HTH")
            elif any(x in markers_top for x in ("HTH_3", "HTH_26")):
                tags.append("Cro-type")
    elif family == "replication":
        if "REP_" in short1 or "REP_" in markers_top:
            tags.append("REP ORF2")
        elif "Rep_trans" in short1:
            tags.append("rep-init")
    elif family == "minimal_integrase":
        tags.append("phage INT-N")
    elif family == "defense_RM":
        tags.append("type I RM")

    phylum = bio_rec.get("bio_top_host_phylum", "")
    try:
        purity = float(host_rec.get("host_purity_phylum", "") or 0)
    except (TypeError, ValueError):
        purity = 0.0
    if phylum:
        if purity >= 0.6:
            tags.append(f"{phylum}↑")
        elif purity <= 0.28:
            tags.append("mixed hosts")
        else:
            tags.append(phylum)

    try:
        ser_frac = float(intdiv_rec.get("integrase_frac_serine_pfam", "") or 0)
    except (TypeError, ValueError):
        ser_frac = 0.0
    if ser_frac >= 0.15:
        tags.append("ser-INT↑")

    if not tags and short1:
        tags.append(short1)
    # dedupe preserving order
    seen: set[str] = set()
    uniq = []
    for t in tags:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    return " · ".join(uniq[:3])


def build_marker_coloring(
    markers: dict[str, dict[str, str]],
    hallmark: dict[str, dict[str, str]],
    bio: dict[str, dict[str, str]],
    host: dict[str, dict[str, str]] | None = None,
    integrase: dict[str, dict[str, str]] | None = None,
) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    all_cids = set(markers) | set(hallmark) | set(bio)
    for cid in all_cids:
        mrec = markers.get(cid, {})
        short1 = mrec.get("markers_top1_short", "")
        top_role = _top_hallmark_role(hallmark, cid)
        family = _infer_marker_family(short1, top_role)
        color = MARKER_FAMILY_PALETTE.get(family, MARKER_FAMILY_PALETTE["other"])
        n_ie = bio.get(cid, {}).get("bio_n_ie", "")
        phylum = bio.get(cid, {}).get("bio_top_host_phylum", "")
        long_parts = [
            f"C{cid} ({n_ie} IE)",
            f"top marker: {short1 or 'n/a'}",
            f"hallmark: {top_role or 'n/a'}",
            f"host: {phylum or 'n/a'}",
        ]
        if mrec.get("markers_top"):
            long_parts.append(f"markers: {mrec['markers_top'][:200]}")
        family_human = MARKER_FAMILY_HUMAN.get(family, MARKER_FAMILY_HUMAN["other"])
        subtype = _cluster_subtype_human(
            family,
            mrec,
            bio.get(cid, {}),
            (host or {}).get(cid, {}),
            (integrase or {}).get(cid, {}),
        )
        label_human = f"C{cid} · {family_human}"
        if n_ie:
            label_human = f"C{cid} ({n_ie}) · {family_human}"
        if subtype:
            label_human = f"{label_human} · {subtype}"
        rec = {
            "marker_label_short": short1 or top_role or f"C{cid}",
            "marker_family": family,
            "marker_family_human": family_human,
            "marker_subtype_human": subtype,
            "marker_label_human": label_human,
            "marker_color": color,
            "marker_label_long": ef._cap(" — ".join(long_parts)),
            "markers_top_roles": top_role,
            # Cytoscape aliases (same convention as filtered_ie pipeline)
            "archetype_short": label_human,
            "archetype_family": family,
            "archetype_color": color,
            "archetype_long": ef._cap(" — ".join(long_parts)),
            "node_color_by": "marker_family",
        }
        out[cid] = rec
    return out


def parse_n_ie(attr: dict) -> float:
    """Cluster size (# IE) from graph node attributes."""
    for key in ("node_weight", "bio_n_ie", "bio_n_IE", "n_reps", "n_reps_in_graph", "ie_len_n"):
        val = attr.get(key)
        if val is None or val == "":
            continue
        try:
            f = float(val)
            if f > 0:
                return f
        except (TypeError, ValueError):
            continue
    return 0.0


def derive_node_visual_attrs(graph: nx.Graph) -> None:
    """Derived node sizing fields (area ∝ n_IE), mirroring edge_width for edges."""
    n_map = {n: parse_n_ie(a) for n, a in graph.nodes(data=True)}
    if not n_map:
        return
    max_n = max(n_map.values()) or 1.0
    for n, n_ie in n_map.items():
        frac = n_ie / max_n
        sqrt_frac = math.sqrt(frac) if frac > 0 else 0.0
        graph.nodes[n]["node_weight"] = float(n_ie)
        graph.nodes[n]["node_area_fraction"] = float(round(frac, 6))
        # Cytoscape: map Size → node_size (typical visual range ~24–220)
        graph.nodes[n]["node_size"] = float(round(24 + 196 * sqrt_frac, 2))
        # matplotlib scatter s (points²): area proportional to n_IE
        graph.nodes[n]["node_size_matplotlib"] = float(round(250 + 1750 * frac, 1))


def enrich_one(graphml_path: Path, out_path: Path | None = None) -> dict[str, int]:
    graph = nx.read_graphml(graphml_path)
    meta = ART

    bio = load_stage8_biology(meta / "cluster_biology_c6960_unk.tsv")
    padloc_summary = load_stage8_padloc_summary(meta / "cluster_padloc_summary_c6960_unk.tsv")
    padloc_systems = build_padloc_systems_from_per_ie(meta / "padloc_per_ie_c6960_unk.tsv")
    hallmark, hallmark_attrs = load_stage8_hallmark(meta / "cluster_hallmark_composition_c6960_unk.tsv")
    markers = load_stage8_markers(meta / "cluster_pfam_markers_c6960_unk.tsv")
    host = load_stage8_host_purity(meta / "cluster_host_purity_c6960_unk.tsv")
    diversity = load_stage8_diversity(meta / "cluster_diversity_c6960_unk.tsv")
    integrase = load_stage8_integrase_div(meta / "integrase_diversity_by_cluster_c6960_unk.tsv")
    atlas_full = load_stage8_atlas_medians(
        REPO / "final_final/artifacts/c6960_unk_clusters_rebuilt.tsv",
        atlas_join_path=REPO / "results_organized/genome_clustering/project/notebook7_extended/ie_instance_classification.tsv",
        leiden_path=meta / "leiden_clusters_c6960_unk.tsv",
    )
    bio_frac = _hallmark_bio_frac_aliases(hallmark)
    coloring = build_marker_coloring(markers, hallmark, bio, host=host, integrase=integrase)

    # ie_len + atlas aliases from biology; merge full atlas medians when available
    ie_len: dict[str, dict[str, str]] = {}
    atlas: dict[str, dict[str, str]] = {}
    for cid, rec in bio.items():
        ie_len[cid] = {
            "ie_len_n": rec.get("bio_n_ie", ""),
            "ie_len_median": rec.get("bio_ie_len_median", ""),
            "ie_len_q25": rec.get("bio_ie_len_q25", ""),
            "ie_len_q75": rec.get("bio_ie_len_q75", ""),
        }
        atlas[cid] = dict(atlas_full.get(cid, {}))
        if not atlas[cid]:
            atlas[cid] = {
                "atlas_median_n_ie": rec.get("bio_n_ie", ""),
                "atlas_median_IE_core_score": rec.get("bio_median_IE_core_score", ""),
            }
        elif "atlas_median_IE_core_score" not in atlas[cid]:
            atlas[cid]["atlas_median_IE_core_score"] = rec.get("bio_median_IE_core_score", "")

    # Merge bio_frac into bio records for legacy compatibility
    for cid, frac_rec in bio_frac.items():
        bio.setdefault(cid, {}).update(frac_rec)

    extra_luts = (host, diversity, integrase, bio_frac, coloring, atlas_full)
    stats = ef.enrich_graph(
        graph,
        bio=bio,
        padloc_summary=padloc_summary,
        padloc_systems=padloc_systems,
        hallmark=hallmark,
        hallmark_attrs=hallmark_attrs,
        products={},
        markers=markers,
        ie_len=ie_len,
        atlas=atlas,
        archetype=coloring,
    )

    for attr in sorted({k for lut in extra_luts for rec in lut.values() for k in rec}):
        ef._ensure_attr(graph, attr)
    for node, ndata in graph.nodes(data=True):
        cid = ef._norm_cluster(ndata.get("cluster_id", ndata.get("cluster_id_str", "")))
        if cid is None:
            continue
        for lut in extra_luts:
            rec = lut.get(cid)
            if rec:
                ndata.update(rec)

    note = (
        "Stage 9 unclassified Leiden rebuild enrichment: cluster_biology_c6960_unk.tsv, "
        "cluster_pfam_markers_c6960_unk.tsv, cluster_hallmark_composition_c6960_unk.tsv, "
        "cluster_padloc_summary_c6960_unk.tsv + padloc_per_ie_c6960_unk.tsv, "
        "cluster_host_purity_c6960_unk.tsv, cluster_diversity_c6960_unk.tsv, "
        "integrase_diversity_by_cluster_c6960_unk.tsv. "
        "Node colour: marker_family / archetype_color (top PFAM marker + hallmark role). "
        "Node label: marker_label_human (C# · marker_family_human). "
        "Node size: node_size / node_size_matplotlib (area ∝ node_weight = n_IE)."
    )
    derive_node_visual_attrs(graph)
    graph.graph["cluster_enrichment_note"] = note
    existing = graph.graph.get("description", "")
    if existing and "cluster_enrichment_note" not in existing:
        graph.graph["description"] = f"{existing} | {note}"
    elif not existing:
        graph.graph["description"] = note

    dest = out_path or graphml_path
    nx.write_graphml(graph, str(dest))
    n_attrs = len(next(iter(graph.nodes(data=True)))[1]) if graph.number_of_nodes() else 0
    print(f"[enrich-stage8] {graphml_path.name} -> {dest.name}  nodes={graph.number_of_nodes()} node_attrs≈{n_attrs}")
    stats["node_attrs"] = n_attrs
    return stats


def main() -> None:
    for name in GRAPHMLS:
        path = ART / name
        if path.is_file():
            enrich_one(path)
        else:
            print(f"[skip] missing {path}")
    print("[done] enriched graphmls in", ART)


if __name__ == "__main__":
    main()
