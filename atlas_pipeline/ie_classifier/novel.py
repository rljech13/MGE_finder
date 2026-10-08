from typing import Any, Dict, List, Tuple


def assess_novel(
    primary: str,
    scores: Dict[str, float],
    cfg: Dict[str, Any],
    conflicts: List[str],
    ann_cov: float,
    genes: List[Dict[str, Any]],
) -> Tuple[bool, str, str]:
    """Return (nonstandard_IE_candidate, novel_tier, novel_rationale)."""
    _ = conflicts
    _ = genes
    crit = cfg.get("novel_outlier_criteria") or {}
    spec_keys: List[str] = list(cfg.get("specialized_modules_for_novel_filter") or [])
    iec = float(scores.get("IE_core_score", 0.0))
    host = float(scores.get("host_contam_score", 0.0))
    max_host = float(crit.get("max_host_contam_for_novel", 1.0))

    reasons: List[str] = []
    candidate = False

    min_ie_uncl = float(crit.get("min_ie_core_unclassified", 0.35))
    min_ann = float(crit.get("min_annotation_coverage", 0.4))
    if primary == "other.unclassified_IE" and iec >= min_ie_uncl and ann_cov >= min_ann and host <= max_host:
        candidate = True
        reasons.append("unclassified_with_strong_IE_core_and_coverage")

    thr_c = float(crit.get("min_module_score_conflict", 0.35))
    n_strong = sum(1 for k in spec_keys if float(scores.get(k, 0.0)) >= thr_c)
    need_labels = int(crit.get("conflict_labels_min", 2))
    if n_strong >= need_labels and host <= max_host:
        candidate = True
        reasons.append(f"multi_module_conflict_n={n_strong}>={need_labels}")

    all_below = float(crit.get("all_specialized_below", 0.22))
    min_ann_weak = float(crit.get("min_annotation_coverage_for_weak_module_novel", min_ann * 0.75))
    if primary == "other.unclassified_IE" and spec_keys and all(float(scores.get(k, 0.0)) < all_below for k in spec_keys):
        if iec >= min_ie_uncl * 0.85 and ann_cov >= min_ann_weak and host <= max_host:
            candidate = True
            reasons.append("all_specialized_modules_weak")

    # Do not mark novel from routine notes (e.g. packaging_only anti_phage) or duplicate
    # strong_modules — those are already covered by multi_module_conflict / scores.

    tier = "none"
    hi = crit.get("high_interest_tier_if") or {}
    if candidate:
        top_spec = max((float(scores.get(k, 0.0)) for k in spec_keys), default=0.0)
        if iec >= float(hi.get("min_ie_core", 0.45)) and top_spec <= float(hi.get("max_top_specialized", 0.35)):
            tier = "high_interest"
        elif primary == "other.unclassified_IE" and iec >= min_ie_uncl:
            tier = "medium"
        else:
            tier = "watchlist"

    rationale = "; ".join(reasons) if reasons else ""
    return candidate, tier, rationale


def module_conflict_note(scores: Dict[str, float], cfg: Dict[str, Any]) -> str:
    crit = cfg.get("novel_outlier_criteria") or {}
    spec_keys: List[str] = list(cfg.get("specialized_modules_for_novel_filter") or [])
    thr_c = float(crit.get("min_module_score_conflict", 0.35))
    strong = [k for k in spec_keys if float(scores.get(k, 0.0)) >= thr_c]
    if len(strong) >= int(crit.get("conflict_labels_min", 2)):
        return "strong_modules:" + ",".join(strong)
    return ""
