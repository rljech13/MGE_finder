#!/usr/bin/env python3
"""Find T/AT pIEs that are still 'small' if replication genes are free.

Rule: n_ta + n_other ≤ 3. Integrase and Rep do not count.
Unknown atlas class only.
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ART = FF / "artifacts/unknown_le3_ta"
ATLAS = FF / "artifacts/attL_complement_blast_window/cohort_6960_atlas_join.tsv"
FUNC = FF / "artifacts/attL_complement_blast_window/cohort_6960_four_func_classes.tsv"
IDX = REPO / "results_organized_rep_6960/atlas_rep_index.tsv"
PAD = FF / "artifacts/nb04_cluster_unclassified/padloc_per_ie_stage8_unc.tsv"
REP_ROOT = REPO / "results_organized_rep_6960"
CAT = FF / "artifacts/unknown_le3_extra/catalog.tsv"

INT_RE = re.compile(
    r"integrase|recombinase|xerc|xerd|tyrosine-type|tyrosine recombinase",
    re.I,
)
REP_RE = re.compile(
    r"\brep\b|replication initiation|replication initiator|initiator rep|"
    r"rep protein|repa\b|repc\b|repabc|primase|primpol|dnag|"
    r"bacteriophage replication gene a|gpa\b|"
    r"plasmid replication",
    re.I,
)
HEL_REP_RE = re.compile(r"helicase rep|rep helicase|trw[c]|tr[ai] helicase", re.I)
TA_RE = re.compile(
    r"toxin[- ]antitoxin|antitoxin|\btoxin\b|addiction module|"
    r"\bhic[ab]\b|\brel[eb]\b|\bpar[ed]\b|\bmaz[ef]\b|\bpem[ik]\b|"
    r"\bvap[bc]\b|\bccd[ab]\b|\bhig[ab]\b|\bmqs[ar]\b|"
    r"abiei|abieii|abigii|\bphd\b|\bdoc\b|\byoeb\b|\byefm\b|"
    r"\bbrn[ta]\b|\bhip[ab]\b|zeta toxin",
    re.I,
)
PAD_TA_RE = re.compile(r"AbiE", re.I)


def find_gbff(row) -> Path | None:
    rel = str(getattr(row, "gbff_path", "") or "")
    realm = str(getattr(row, "source_realm", "") or "")
    cands = []
    if rel and rel != "nan":
        cands.append(REP_ROOT / rel)
        stem = Path(rel).name.replace(".annotation.gbff.merged", "")
        if realm:
            cands.append(REP_ROOT / realm / stem / "annotation.gbff.merged")
    for p in cands:
        if p.is_file() and p.stat().st_size > 0:
            return p
    return None


def kind_of(prod: str) -> str:
    p = prod or ""
    if INT_RE.search(p):
        return "int"
    if TA_RE.search(p):
        return "ta"
    if REP_RE.search(p) or HEL_REP_RE.search(p):
        return "rep"
    return "other"


def ta_role(prod: str) -> str:
    p = (prod or "").lower()
    if "antitoxin" in p or "addiction module" in p or "abiei" in p:
        return "antitoxin"
    if "toxin" in p or "abieii" in p or "abigii" in p:
        return "toxin"
    if TA_RE.search(prod or ""):
        return "ta"
    return ""


def parse_gbff(path: Path) -> list[dict]:
    rec = SeqIO.read(path, "genbank")
    rows = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        prod = (f.qualifiers.get("product") or [""])[0]
        rows.append({
            "start": int(f.location.start) + 1,
            "end": int(f.location.end),
            "strand": f.location.strand or 1,
            "product": prod,
            "kind": kind_of(prod),
            "aa": (int(f.location.end) - int(f.location.start)) // 3,
            "role": ta_role(prod),
        })
    return rows


def padloc_ta_spans(pad_csv: Path) -> list[tuple[int, int, str]]:
    if not pad_csv.is_file() or pad_csv.stat().st_size == 0:
        return []
    try:
        df = pd.read_csv(pad_csv)
    except Exception:
        return []
    if "system" not in df.columns:
        return []
    out = []
    for r in df.itertuples(index=False):
        sys = str(getattr(r, "system", "") or "")
        if not PAD_TA_RE.search(sys):
            continue
        try:
            a, b = int(r.start), int(r.end)
        except Exception:
            continue
        if a > b:
            a, b = b, a
        out.append((a, b, sys))
    return out


def overlay_padloc(genes: list[dict], spans: list[tuple[int, int, str]]) -> None:
    for g in genes:
        if g["kind"] == "ta":
            continue
        for a, b, _sys in spans:
            ov = min(g["end"], b) - max(g["start"], a)
            if ov > 0.5 * min(g["end"] - g["start"] + 1, b - a + 1):
                g["kind"] = "ta"
                if not g["role"]:
                    g["role"] = "ta"
                break


def main() -> None:
    ART.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(ATLAS, sep="\t")
    func = pd.read_csv(FUNC, sep="\t")
    idx = pd.read_csv(IDX, sep="\t")
    pad = pd.read_csv(PAD, sep="\t")
    m = atlas.merge(func[["integrase_id", "func_class"]], on="integrase_id")
    m = m.merge(idx[["integrase_id", "gbff_path", "source_realm"]], on="integrase_id", how="left")
    m = m.merge(pad[["integrase_id", "padloc_systems"]], on="integrase_id", how="left")
    in_le3 = set()
    if CAT.is_file():
        in_le3 = set(pd.read_csv(CAT, sep="\t").integrase_id)

    m["has_padloc_ta"] = m["padloc_systems"].fillna("").str.contains("AbiE", regex=False)
    m = m[m.func_class == "unknown"].copy()
    cand = m[(m.has_padloc_ta) | (m.n_cds.fillna(99) <= 10)].copy()
    print("candidates to parse", len(cand), "padloc_AbiE", int(m.has_padloc_ta.sum()))

    hits = []
    genes_out = []
    n_ok = n_miss = 0
    for r in cand.itertuples():
        p = find_gbff(r)
        if p is None:
            n_miss += 1
            continue
        n_ok += 1
        genes = parse_gbff(p)
        pad_cands = [p.parent / "padloc.csv"]
        rel = str(getattr(r, "gbff_path", "") or "")
        realm = str(getattr(r, "source_realm", "") or "")
        if rel and realm:
            stem = Path(rel).name.replace(".annotation.gbff.merged", "")
            pad_cands.append(REP_ROOT / realm / stem / "padloc.csv")
        spans = []
        for pc in pad_cands:
            spans = padloc_ta_spans(pc)
            if spans:
                break
        overlay_padloc(genes, spans)
        n_int = sum(g["kind"] == "int" for g in genes)
        n_rep = sum(g["kind"] == "rep" for g in genes)
        n_ta = sum(g["kind"] == "ta" for g in genes)
        n_other = sum(g["kind"] == "other" for g in genes)
        if n_int == 0 and genes:
            n_int = 1
            n_other = max(0, n_other - 1)
        n_budget = n_ta + n_other
        if n_ta < 1:
            continue
        if n_budget > 3:
            continue
        ta_prod = "; ".join(
            f"{g['product']}({g['aa']}aa,{g['start']}-{g['end']})"
            for g in genes if g["kind"] == "ta"
        )
        rep_prod = "; ".join(
            f"{g['product']}({g['aa']}aa)" for g in genes if g["kind"] == "rep"
        )
        oth_prod = "; ".join(
            f"{g['product']}({g['aa']}aa)" for g in genes if g["kind"] == "other"
        )
        hits.append({
            "integrase_id": r.integrase_id,
            "func_class": r.func_class,
            "ie_id": r.ie_id,
            "n_cds_gbff": len(genes),
            "n_int": n_int,
            "n_rep": n_rep,
            "n_ta": n_ta,
            "n_other": n_other,
            "n_budget": n_budget,
            "padloc": r.padloc_systems if pd.notna(r.padloc_systems) else "",
            "in_unknown_le3": r.integrase_id in in_le3,
            "ta_products": ta_prod,
            "rep_products": rep_prod,
            "other_products": oth_prod,
            "n_cds_atlas": r.n_cds,
        })
        for g in genes:
            genes_out.append({
                "integrase_id": r.integrase_id,
                **g,
            })

    out = pd.DataFrame(hits)
    if out.empty:
        print("NO HITS")
        return
    lens = pd.read_csv(
        FF / "artifacts/attL_complement_blast_window/cohort_6960_ie_lengths.tsv",
        sep="\t",
    )
    if "integrase_id" in lens.columns and "ie_len_nt" in lens.columns:
        out = out.merge(lens[["integrase_id", "ie_len_nt"]], on="integrase_id", how="left")
    out = out.sort_values(["n_budget", "n_rep", "n_ta", "n_cds_gbff", "ie_len_nt"])
    dest = ART / "ta_allow_rep.tsv"
    out.to_csv(dest, sep="\t", index=False)
    pd.DataFrame(genes_out).to_csv(ART / "ta_allow_rep_genes.tsv", sep="\t", index=False)
    print("parsed", n_ok, "missing_gbff", n_miss)
    print("hits", len(out), "in_le3", int(out.in_unknown_le3.sum()),
          "new", int((~out.in_unknown_le3).sum()))
    print(out.groupby(["in_unknown_le3", "n_rep"]).size().to_string())
    pd.set_option("display.max_colwidth", 100)
    pd.set_option("display.width", 240)
    cols = [
        "integrase_id", "n_cds_gbff", "n_int", "n_rep", "n_ta", "n_other",
        "n_budget", "ie_len_nt", "padloc", "in_unknown_le3",
        "ta_products", "rep_products", "other_products",
    ]
    print(out[cols].to_string(index=False))
    print("\nwrote", dest)


if __name__ == "__main__":
    main()
