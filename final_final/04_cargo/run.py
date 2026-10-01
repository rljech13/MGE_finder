#!/usr/bin/env python3
"""Check inputs, then run one module experiment.

From this directory:

  python run.py check
  python run.py check 03_class_inventory
  python run.py run 01_rm_ta_cohort

``check`` only looks for files. ``run`` executes the steps in order and stops
at the first missing input. A step whose inputs are the outputs of an earlier
step in the same experiment is still listed, so the chain is visible.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "03_class_inventory"))
from mge_paths import REPO  # noqa: E402

ART = "final_final/artifacts"
FIG = "final_final/figures"
COHORT = f"{ART}/attL_complement_blast_window"
UNK = f"{ART}/unknown_ime_classes"
ANN = f"{ART}/unknown_reannotate"

SHARED = [
    f"{COHORT}/cohort_6960_atlas_join.tsv",
    f"{COHORT}/cohort_6960_four_func_classes.tsv",
    "results_organized_rep_6960/atlas_rep_index.tsv",
    "results_organized/misc/gtdb_taxonomy.tsv",
    f"{ART}/nb04_cluster_unclassified/padloc_per_ie_stage8_unc.tsv",
]

EXPERIMENTS: dict[str, list[tuple[str, list[str]]]] = {
    "01_rm_ta_cohort": [
        ("scan_unknown_rm_ta.py", SHARED),
        ("plot_unknown_rm_ta_diversity.py", [
            f"{ART}/unknown_rm/unknown_rm_hits.tsv",
            f"{ART}/unknown_ta/unknown_ta_hits.tsv",
        ]),
    ],
    "02_small_rm_ta": [
        ("scan_small_rm_allow_rep.py", SHARED),
        ("scan_small_ta_allow_rep.py", SHARED),
        ("plot_small_rm_histograms.py", [
            f"{ART}/unknown_le3_rm/rm_allow_rep.tsv",
            f"{COHORT}/cohort_6960_atlas_join.tsv",
            "results_organized/misc/gtdb_taxonomy.tsv",
        ]),
        ("plot_small_ta_histograms.py", [
            f"{ART}/unknown_le3_ta/ta_allow_rep.tsv",
            f"{ART}/unknown_le3_ta/ta_allow_rep_genes.tsv",
            f"{COHORT}/cohort_6960_atlas_join.tsv",
            "results_organized/misc/gtdb_taxonomy.tsv",
        ]),
    ],
    "03_class_inventory": [
        ("count_unknown_minimal_ie.py", SHARED + [
            f"{COHORT}/cohort_6960_four_class_v2.tsv",
            f"{ART}/unknown_rm/unknown_rm_ta_scan.tsv",
        ]),
        ("score_unknown_module_minimal.py", [
            f"{UNK}/unknown_minimal_ie.tsv",
            f"{ANN}/proteins.tsv",
            f"{ART}/unknown_le3_rm/rm_allow_rep.tsv",
            f"{ART}/unknown_le3_ta/ta_allow_rep.tsv",
        ]),
        ("recount_unknown_minimal_annotated.py", [
            f"{UNK}/unknown_minimal_ie.tsv",
            f"{ANN}/unknown_all_annotated.tsv",
            f"{ART}/nb04_cluster_unclassified/padloc_per_ie_stage8_unc.tsv",
        ]),
        ("recount_unknown_minimal_after_hypothesis.py", [
            f"{UNK}/unknown_minimal_ie.tsv",
            f"{ANN}/unknown_all_annotated.tsv",
            f"{ANN}/still_dark_all/still_dark_calls.tsv",
            f"{ANN}/still_dark_all/leftover/leftover_calls.tsv",
            f"{ANN}/still_dark_all/hypothesis/hypothesis_calls.tsv",
            f"{ANN}/still_dark_all/hypothesis/cis_elements.tsv",
        ]),
        ("apply_duf_as_unannotated.py", [
            f"{UNK}/after_hypothesis/unknown_minimal_le1_long.tsv",
            f"{ANN}/unknown_all_annotated.tsv",
        ]),
        ("inventory_unknown_funnel_classes.py", [
            f"{UNK}/after_hypothesis/unknown_minimal_ie.tsv",
            f"{UNK}/after_hypothesis/by_class.tsv",
        ]),
        ("plot_unknown_module_minimal.py", [
            f"{UNK}/after_hypothesis/unknown_minimal_ie.tsv",
            f"{COHORT}/cohort_6960_atlas_join.tsv",
        ]),
        ("plot_unknown_minimal_diversity.py", [
            f"{UNK}/after_hypothesis/unknown_minimal_ie.tsv",
            f"{COHORT}/cohort_6960_atlas_join.tsv",
        ]),
        ("plot_unknown_minimal_histograms.py", [
            f"{UNK}/after_hypothesis/unknown_minimal_le1_long.tsv",
        ]),
    ],
}


def missing(paths: list[str]) -> list[str]:
    return [rel for rel in paths if not (REPO / rel).exists()]


def check(name: str) -> bool:
    ok = True
    print(name)
    for script, needs in EXPERIMENTS[name]:
        gone = missing(needs)
        if gone:
            ok = False
            print(f"  {script}")
            for rel in gone:
                print(f"    missing  {rel}")
        else:
            print(f"  {script}  inputs ok")
    return ok


def run(name: str) -> int:
    folder = Path(__file__).resolve().parent / name
    if not check(name):
        print("\nnot started: an input above is missing")
        return 1
    for script, _needs in EXPERIMENTS[name]:
        print(f"\n=== {name}/{script}")
        result = subprocess.run([sys.executable, str(folder / script)], cwd=folder)
        if result.returncode != 0:
            return result.returncode
    return 0


def main() -> int:
    args = sys.argv[1:]
    if not args or args[0] not in {"check", "run"}:
        print(__doc__)
        return 2
    chosen = args[1:] or list(EXPERIMENTS)
    unknown = [name for name in chosen if name not in EXPERIMENTS]
    if unknown:
        print("unknown experiment:", ", ".join(unknown))
        print("known:", ", ".join(EXPERIMENTS))
        return 2
    if args[0] == "check":
        failed = [name for name in chosen if not check(name)]
        return 1 if failed else 0
    for name in chosen:
        code = run(name)
        if code != 0:
            return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
