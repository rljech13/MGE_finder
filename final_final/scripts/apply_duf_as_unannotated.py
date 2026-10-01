#!/usr/bin/env python3
"""Treat non-obligatory DUF as unannotated.

A DUF family is obligatory for a class if it occurs on ≥70% of unknown IEs
that carry that class. Class-hallmark DUFs (DUF1016 = PDC, DUF4238 = Rst)
are the class itself, not unannotated.

Minimal (locked + this rule):
  one class, no leftover, short <70 aa free,
  ≤1 unannotated CDS ≥70 aa,
  cannot have a non-obligatory DUF and another unannotated ≥70 aa together.
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import re
from pathlib import Path

import pandas as pd

OUT = FF / "artifacts/unknown_ime_classes/after_hypothesis"
IE = OUT / "unknown_minimal_le1_long.tsv"
ANN = FF / "artifacts/unknown_reannotate/unknown_all_annotated.tsv"

DUF_RE = re.compile(r"DUF(\d+)", re.I)
HALLMARK = {
    "DUF1016": "PDC",
    "DUF4238": "Rst_DUF4238",
}
CLASS_ORDER = [
    "replication", "TA", "RM", "PDC", "RT-defense", "AMR", "Par", "metals",
]


def duf_family(pfam: object, bakta: object) -> str:
    for text in (pfam, bakta):
        if text is None or (isinstance(text, float) and pd.isna(text)):
            continue
        m = DUF_RE.search(str(text))
        if m:
            return f"DUF{m.group(1)}"
    return ""


def parse_classes(s: object) -> set[str]:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return set()
    return {x for x in str(s).split("|") if x}


def main() -> None:
    ie = pd.read_csv(IE, sep="\t")
    ann = pd.read_csv(
        ANN, sep="\t",
        usecols=["protein_id", "integrase_id", "aa_len", "bakta_product", "pfam_hmm"],
    )
    ann["duf_family"] = [
        duf_family(r.pfam_hmm, r.bakta_product) for r in ann.itertuples()
    ]
    duf = ann[ann.duf_family != ""].copy()
    duf["hallmark"] = duf.duf_family.map(HALLMARK).fillna("")

    class_ie: dict[str, set[str]] = {}
    for r in ie.itertuples():
        for c in parse_classes(r.classes):
            class_ie.setdefault(c, set()).add(r.integrase_id)

    # family × class frequency among IEs that carry the class
    freq_rows = []
    oblig: set[tuple[str, str]] = set()
    for fam, sub in duf.groupby("duf_family"):
        ies_with_fam = set(sub.integrase_id)
        for cls, members in class_ie.items():
            n_cls = len(members)
            n_both = len(ies_with_fam & members)
            if n_cls == 0 or n_both == 0:
                continue
            frac = n_both / n_cls
            is_hall = HALLMARK.get(fam, "") == cls or (fam == "DUF1016" and cls == "PDC")
            is_ob = (not is_hall) and frac >= 0.70
            if is_ob:
                oblig.add((cls, fam))
            freq_rows.append({
                "duf_family": fam,
                "class": cls,
                "n_ie_with_class": n_cls,
                "n_ie_with_family_and_class": n_both,
                "frac": round(frac, 4),
                "obligatory_ge70pct": int(is_ob),
                "class_hallmark": int(bool(is_hall)),
            })
    freq = pd.DataFrame(freq_rows)
    if len(freq):
        freq = freq.sort_values(["obligatory_ge70pct", "frac", "class"], ascending=[False, False, True])
    freq.to_csv(OUT / "duf_family_by_class.tsv", sep="\t", index=False)

    duf_by_ie: dict[str, list[dict]] = {}
    for r in duf.itertuples():
        duf_by_ie.setdefault(r.integrase_id, []).append({
            "family": r.duf_family,
            "aa_len": int(r.aa_len),
            "hallmark": r.hallmark,
        })

    n_duf_ge70 = []
    n_duf_lt70 = []
    n_duf_any = []
    duf_ex = []
    n_unann = []
    both = []
    keep = []
    for r in ie.itertuples():
        cls = r.char_class or ""
        recs = duf_by_ie.get(r.integrase_id, [])
        nonob = []
        for x in recs:
            if x["hallmark"]:
                continue
            if (cls, x["family"]) in oblig:
                continue
            nonob.append(x)
        n_ge = sum(1 for x in nonob if x["aa_len"] >= 70)
        n_lt = sum(1 for x in nonob if x["aa_len"] < 70)
        n_unres70 = int(r.n_unresolved_ge70)
        n_un = n_unres70 + n_ge
        has_both = int(n_unres70 >= 1 and len(nonob) >= 1)
        n_duf_ge70.append(n_ge)
        n_duf_lt70.append(n_lt)
        n_duf_any.append(len(nonob))
        duf_ex.append(";".join(f"{x['family']}({x['aa_len']}aa)" for x in nonob[:6]))
        n_unann.append(n_un)
        both.append(has_both)
        keep.append(int(r.composition_minimal == 1 and n_un <= 1 and has_both == 0))

    ie["n_duf_nonoblig_ge70"] = n_duf_ge70
    ie["n_duf_nonoblig_lt70"] = n_duf_lt70
    ie["n_duf_nonoblig"] = n_duf_any
    ie["duf_nonoblig_examples"] = duf_ex
    ie["n_unannotated_ge70"] = n_unann
    ie["duf_plus_unannotated"] = both
    ie["minimal_le1_long"] = keep
    ie["minimal_no_long_dark"] = (ie.composition_minimal.eq(1) & ie.n_unannotated_ge70.eq(0)).astype(int)
    ie.to_csv(IE, sep="\t", index=False)

    m = ie[ie.minimal_le1_long == 1]
    n0 = int((m.n_unannotated_ge70 == 0).sum())
    n1 = int((m.n_unannotated_ge70 == 1).sum())
    dropped = ie[(ie.composition_minimal == 1) & (ie.minimal_le1_long == 0)]
    summary = pd.DataFrame([
        {"metric": "unknown_n", "n": len(ie)},
        {"metric": "composition_minimal_one_class_no_leftover", "n": int(ie.composition_minimal.sum())},
        {"metric": "minimal_le1_long_unannotated", "n": len(m)},
        {"metric": "of_them_zero_unannotated_ge70", "n": n0},
        {"metric": "of_them_one_unannotated_ge70", "n": n1},
        {"metric": "dropped_from_462_by_duf_rule", "n": len(dropped)},
        {"metric": "dropped_duf_plus_unannotated", "n": int((dropped.duf_plus_unannotated == 1).sum())},
        {"metric": "dropped_unann_ge70_ge2", "n": int((dropped.n_unannotated_ge70 >= 2).sum())},
        {"metric": "duf_families_obligatory", "n": len(oblig)},
        {"metric": "rep_only", "n": int((m.char_class == "replication").sum())},
        {"metric": "module_plus_rep", "n": int(((m.char_class != "replication") & (m.has_rep == 1)).sum())},
        {"metric": "module_no_rep", "n": int(((m.char_class != "replication") & (m.has_rep == 0)).sum())},
    ])
    for c in CLASS_ORDER:
        summary = pd.concat([summary, pd.DataFrame([{
            "metric": c, "n": int((m.char_class == c).sum()),
        }])], ignore_index=True)
    summary.to_csv(OUT / "minimal_le1_long_summary.tsv", sep="\t", index=False)

    print("obligatory (class, family):", sorted(oblig) or "none")
    print(summary.to_string(index=False))
    print("\nminimal by class / 0 vs 1 unann≥70:")
    for c in CLASS_ORDER:
        sub = m[m.char_class == c]
        if sub.empty:
            continue
        print(f"  {c:12s} n={len(sub):3d}  0={int((sub.n_unannotated_ge70==0).sum())}  1={int((sub.n_unannotated_ge70==1).sum())}")
    print("\ndropped composition-min reasons (may overlap):")
    print("  duf+unann", int((dropped.duf_plus_unannotated == 1).sum()))
    print("  unann≥2", int((dropped.n_unannotated_ge70 >= 2).sum()))
    print("\nfreq rows with frac≥0.10:")
    if len(freq):
        print(freq[freq.frac >= 0.10][["duf_family", "class", "n_ie_with_class", "n_ie_with_family_and_class", "frac", "obligatory_ge70pct", "class_hallmark"]].to_string(index=False))
    print("wrote", IE)


if __name__ == "__main__":
    main()
