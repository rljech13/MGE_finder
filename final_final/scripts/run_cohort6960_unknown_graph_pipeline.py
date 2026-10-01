#!/usr/bin/env python3
"""Full diploma-style graph clustering for cohort-6960 Figure-4 Unknown.

Cohort: partial→unknown; ICE=T4SS only; Unknown n=4,774 from
  artifacts/attL_complement_blast_window/cohort_6960_atlas_classes_conj_split.tsv

Steps:
  1. build_cohort6960_unknown_clustering.py
  2. recompute_hallmark_hmm_imputed_c6960_unk.py
  3. plot_integrase_diversity_c6960_unk.py
  4. render_cluster_markers_table_c6960_unk.py
  5. render_cluster_markers_heatmap_c6960_unk.py
  6. build_cluster_network_c6960_unk.py
  7. enrich_cluster_overview_c6960_unk.py
  8. enrich_c6960_unk_ie_knn_graphml.py

Usage:
  conda run -n MGE_finder python run_cohort6960_unknown_graph_pipeline.py
  conda run -n MGE_finder python run_cohort6960_unknown_graph_pipeline.py --from-step 3
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent

STEPS = [
    ("clustering", "build_cohort6960_unknown_clustering.py"),
    ("hallmark", "recompute_hallmark_hmm_imputed_c6960_unk.py"),
    ("integrase_div", "plot_integrase_diversity_c6960_unk.py"),
    ("markers_table", "render_cluster_markers_table_c6960_unk.py"),
    ("markers_heatmap", "render_cluster_markers_heatmap_c6960_unk.py"),
    ("cluster_network", "build_cluster_network_c6960_unk.py"),
    ("enrich_overview", "enrich_cluster_overview_c6960_unk.py"),
    ("ie_knn_overlay", "enrich_c6960_unk_ie_knn_graphml.py"),
]


def run_step(script: str) -> None:
    path = SCRIPTS / script
    cmd = [sys.executable, str(path)]
    print("\n" + "=" * 72)
    print(" ".join(cmd))
    print("=" * 72)
    subprocess.run(cmd, check=True, cwd=str(SCRIPTS))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-step", type=int, default=1, help="1-based step index")
    ap.add_argument("--to-step", type=int, default=len(STEPS), help="1-based step index (inclusive)")
    args = ap.parse_args()

    for i, (name, script) in enumerate(STEPS, start=1):
        if i < args.from_step or i > args.to_step:
            continue
        print(f"\n>>> Step {i}/{len(STEPS)}: {name}")
        run_step(script)

    print("\n✓ Cohort-6960 Unknown graph pipeline complete.")
    print("  Artifacts: final_final/artifacts/nb04_cluster_unknown_6960/")
    print("  Figures:   final_final/figures/04_cluster_unknown_6960/")


if __name__ == "__main__":
    main()
