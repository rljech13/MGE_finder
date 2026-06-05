"""Filter rules for attachment-site validation and integrative element confidence scoring."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_SHIFT = 3
DEFAULT_ATTL_MIN_BP = 15
DEFAULT_INTEGRASE_MIN_AA = 300
DEFAULT_IE_MIN_NT = 2000


@dataclass(frozen=True)
class FilterThresholds:
    """Numeric thresholds applied when selecting confident integrative elements.

    Attributes:
        shift: Maximum allowed distance (bp) between a BLAST hit and the tRNA 3-prime end.
        attl_min_bp: Minimum attachment-site (attL) length in base pairs.
        integrase_min_aa: Minimum integrase length in amino acids (exclusive lower bound).
        ie_min_nt: Minimum integrative element length in nucleotides.
        reject_ambiguous_n_ie: When True, reject elements containing ambiguous N bases.
    """

    shift: int = DEFAULT_SHIFT
    attl_min_bp: int = DEFAULT_ATTL_MIN_BP
    integrase_min_aa: int = DEFAULT_INTEGRASE_MIN_AA
    ie_min_nt: int = DEFAULT_IE_MIN_NT
    reject_ambiguous_n_ie: bool = True


def parse_orfs_gff(path: Path | str) -> dict[str, list[tuple[int, int, str, str]]]:
    """Parse Prodigal CDS features from a GFF file.

    Args:
        path: Path to a Prodigal GFF file produced by ``predict_orfs.py``.

    Returns:
        Dictionary mapping contig identifier to a list of CDS tuples
        ``(start_1based, end_1based, strand, cds_id)``.
    """
    cds_by_contig: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    path = Path(path)
    if not path.is_file():
        return cds_by_contig
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            cds_id = ""
            for attr in parts[8].split(";"):
                if attr.startswith("ID="):
                    cds_id = attr[3:]
                    break
            cds_by_contig[parts[0]].append((start, end, strand, cds_id))
    return cds_by_contig


def _best_overlap(
    attl_lo: int,
    attl_hi: int,
    cds_list: list[tuple[int, int, str, str]],
) -> tuple[int, tuple[int, int, str, str] | None]:
    """Return the largest overlap between an interval and a list of CDS intervals.

    Args:
        attl_lo: Start coordinate of the attachment site (1-based inclusive).
        attl_hi: End coordinate of the attachment site (1-based inclusive).
        cds_list: CDS intervals on the same contig.

    Returns:
        Tuple ``(overlap_bp, best_cds)`` where ``best_cds`` is the CDS with the
        largest overlap, or None when ``cds_list`` is empty.
    """
    best_ov = 0
    best_cds = None
    for c0, c1, cstrand, cid in cds_list:
        ov = max(0, min(attl_hi, c1) - max(attl_lo, c0) + 1)
        if ov > best_ov:
            best_ov = ov
            best_cds = (c0, c1, cstrand, cid)
    return best_ov, best_cds


def classify_prodigal_overlap(
    attl_lo: int,
    attl_hi: int,
    cds_list: list[tuple[int, int, str, str]],
) -> dict[str, Any]:
    """Classify whether an attL interval overlaps Prodigal-predicted CDS features.

    Args:
        attl_lo: Start coordinate of attL on the genome (1-based inclusive).
        attl_hi: End coordinate of attL on the genome (1-based inclusive).
        cds_list: Prodigal CDS intervals on the attL contig.

    Returns:
        Dictionary with overlap length, overlap fraction, hit class
        (``intergenic``, ``partial_CDS``, ``fully_inside_CDS``, or
        ``no_cds_on_contig``), and metadata for the best-overlapping CDS.
    """
    attl_len = attl_hi - attl_lo + 1
    if not cds_list:
        return {
            "prodigal_overlap_bp": 0,
            "prodigal_overlap_frac": 0.0,
            "prodigal_hit_class": "no_cds_on_contig",
            "prodigal_hit_cds_id": "",
            "prodigal_hit_cds_strand": "",
        }
    best_ov, best_cds = _best_overlap(attl_lo, attl_hi, cds_list)
    if best_ov == 0:
        hit_class = "intergenic"
    elif best_ov >= attl_len:
        hit_class = "fully_inside_CDS"
    else:
        hit_class = "partial_CDS"
    return {
        "prodigal_overlap_bp": int(best_ov),
        "prodigal_overlap_frac": float(best_ov) / attl_len if attl_len > 0 else 0.0,
        "prodigal_hit_class": hit_class,
        "prodigal_hit_cds_id": best_cds[3] if best_cds else "",
        "prodigal_hit_cds_strand": best_cds[2] if best_cds else "",
    }


def select_v3ps_strict_hit(
    raw_hits: pd.DataFrame,
    trna_len: int,
    trna_strand: str,
    shift: int = DEFAULT_SHIFT,
) -> dict[str, Any] | None:
    """Select the longest BLAST hit that passes strict attachment-site criteria.

    A hit must anchor at the tRNA 3-prime end within ``shift`` base pairs and
    match the expected subject-strand orientation relative to tRNA strand.

    Args:
        raw_hits: Raw BLAST table for one tRNA query (``mge_blast_raw.tsv`` subset).
        trna_len: tRNA query length in nucleotides.
        trna_strand: tRNA strand (``+`` or ``-``).
        shift: Allowed distance from the tRNA 3-prime end (base pairs).

    Returns:
        Dictionary with attL coordinates and hit metadata for the selected hit,
        or None when no hit passes the strict filters.
    """
    if raw_hits.empty:
        return None
    anchor_thr = trna_len - shift + 1
    sub = raw_hits.copy()
    sub["qend"] = sub["qend"].astype(int)
    sub["sstart"] = sub["sstart"].astype(int)
    sub["send"] = sub["send"].astype(int)
    sub["length"] = sub["length"].astype(int)
    sub["wstart"] = sub["wstart"].astype(int)

    kept = sub[sub["qend"] >= anchor_thr]
    if trna_strand == "+":
        kept = kept[kept["sstart"] < kept["send"]]
    else:
        kept = kept[kept["sstart"] > kept["send"]]
    if kept.empty:
        return None

    kept = kept.copy()
    kept["s_lo"] = kept[["sstart", "send"]].min(axis=1)
    kept["s_hi"] = kept[["sstart", "send"]].max(axis=1)
    kept["abs_lo"] = kept["wstart"] + kept["s_lo"] - 1
    kept["abs_hi"] = kept["wstart"] + kept["s_hi"] - 1

    best = kept.loc[kept["length"].idxmax()]
    return {
        "best_attl_len_bp": int(best["length"]),
        "attL_abs_lo": int(best["abs_lo"]),
        "attL_abs_hi": int(best["abs_hi"]),
        "attL_strand": "+" if int(best["sstart"]) < int(best["send"]) else "-",
        "n_v3ps_hits": int(len(kept)),
        "qstart": int(best["qstart"]),
        "qend": int(best["qend"]),
    }


def integrase_aa_length(start: int, end: int) -> int:
    """Estimate integrase length in amino acids from nucleotide coordinates.

    Args:
        start: Integrase start coordinate (1-based inclusive).
        end: Integrase end coordinate (1-based inclusive).

    Returns:
        Integrase length in amino acids, computed as nucleotide span divided by three.
    """
    return (int(end) - int(start) + 1) // 3


def build_qseqid(
    integrase_id: str,
    contig: str,
    trna_start: int,
    trna_end: int,
    trna_strand: str,
) -> str:
    """Build the BLAST query identifier used in ``mge_query.fa`` headers.

    Args:
        integrase_id: Integrase ORF identifier.
        contig: Contig name.
        trna_start: tRNA start coordinate (1-based inclusive).
        trna_end: tRNA end coordinate (1-based inclusive).
        trna_strand: tRNA strand (``+`` or ``-``).

    Returns:
        Query identifier string in the form
        ``integrase_id:contig:start-end:strand``.
    """
    return f"{integrase_id}:{contig}:{trna_start}-{trna_end}:{trna_strand}"


def evaluate_ie_candidate(
    *,
    integrase_id: str,
    contig: str,
    trna_start: int,
    trna_end: int,
    trna_strand: str,
    trna_len: int,
    integrase_start: int,
    integrase_end: int,
    ie_id: str | None,
    ie_len_nt: int | None,
    ie_has_n: bool,
    raw_hits: pd.DataFrame | None,
    cds_by_contig: dict[str, list[tuple[int, int, str, str]]],
    thresholds: FilterThresholds,
) -> dict[str, Any]:
    """Evaluate one integrative element candidate against all confidence filters.

    Args:
        integrase_id: Integrase ORF identifier.
        contig: Contig carrying the integrase and tRNA pair.
        trna_start: tRNA start coordinate (1-based inclusive).
        trna_end: tRNA end coordinate (1-based inclusive).
        trna_strand: tRNA strand (``+`` or ``-``).
        trna_len: tRNA length in nucleotides.
        integrase_start: Integrase start coordinate (1-based inclusive).
        integrase_end: Integrase end coordinate (1-based inclusive).
        ie_id: Extracted element FASTA record identifier, if present.
        ie_len_nt: Extracted element length in nucleotides, if present.
        ie_has_n: True when the extracted element sequence contains ambiguous N.
        raw_hits: Raw BLAST hits for the tRNA query, or None when unavailable.
        cds_by_contig: Prodigal CDS intervals keyed by contig identifier.
        thresholds: Filter thresholds loaded from pipeline configuration.

    Returns:
        Audit dictionary with per-filter pass flags, measured values, and
        ``reject_reason`` when the candidate fails (empty string when it passes).
    """
    aa_len = integrase_aa_length(integrase_start, integrase_end)
    row: dict[str, Any] = {
        "integrase_id": integrase_id,
        "ie_id": ie_id or "",
        "contig": contig,
        "trna_start": trna_start,
        "trna_end": trna_end,
        "trna_strand": trna_strand,
        "trna_len": trna_len,
        "integrase_len_aa": aa_len,
        "ie_len_nt": ie_len_nt if ie_len_nt is not None else 0,
        "ie_has_ambiguous_n": ie_has_n,
        "best_attl_len_bp": 0,
        "attL_abs_lo": 0,
        "attL_abs_hi": 0,
        "attL_strand": "",
        "n_v3ps_hits": 0,
        "prodigal_hit_class": "",
        "prodigal_overlap_bp": 0,
        "reject_reason": "",
        "passed_v3ps_strict": False,
        "passed_attl_min": False,
        "passed_intergenic": False,
        "passed_integrase_len": False,
        "passed_ie_len": False,
        "passed_no_ambiguous_n": False,
        "passed_confident": False,
    }

    if raw_hits is None or raw_hits.empty:
        row["reject_reason"] = "v3ps_no_raw_blast"
        return row

    v3ps = select_v3ps_strict_hit(raw_hits, trna_len, trna_strand, thresholds.shift)
    if v3ps is None:
        row["reject_reason"] = "v3ps_no_strict_hit"
        return row

    row.update(v3ps)
    row["passed_v3ps_strict"] = True

    if v3ps["best_attl_len_bp"] < thresholds.attl_min_bp:
        row["reject_reason"] = "attl_too_short"
        return row
    row["passed_attl_min"] = True

    cds_list = cds_by_contig.get(contig, [])
    if not cds_list:
        for key, val in cds_by_contig.items():
            if contig in key or key in contig:
                cds_list = val
                break

    prodigal = classify_prodigal_overlap(v3ps["attL_abs_lo"], v3ps["attL_abs_hi"], cds_list)
    row.update(prodigal)
    if prodigal["prodigal_hit_class"] != "intergenic":
        row["reject_reason"] = f"cds_{prodigal['prodigal_hit_class']}"
        return row
    row["passed_intergenic"] = True

    if aa_len <= thresholds.integrase_min_aa:
        row["reject_reason"] = "integrase_too_short"
        return row
    row["passed_integrase_len"] = True

    if thresholds.reject_ambiguous_n_ie and ie_has_n:
        row["reject_reason"] = "ambiguous_n"
        return row
    row["passed_no_ambiguous_n"] = True

    if ie_len_nt is None or ie_len_nt < thresholds.ie_min_nt:
        row["reject_reason"] = "ie_too_short"
        return row
    row["passed_ie_len"] = True

    row["passed_confident"] = True
    row["reject_reason"] = ""
    return row
