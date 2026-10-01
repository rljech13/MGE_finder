"""Shared paths for cohort-6960 atlas Unknown graph clustering (diploma-style).

Cohort: Figure-4 class scheme on locked 6,960 IE —
  func_class == unknown (partial ICE → unknown; ICE = T4SS only).
Source table:
  artifacts/attL_complement_blast_window/cohort_6960_atlas_classes_conj_split.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import FF, REPO
ART = FF / "artifacts/nb04_cluster_unknown_6960"
FIG = FF / "figures/04_cluster_unknown_6960"
IDS = FF / "id_lists"
TAG = "c6960_unk"

PROJ_GENO = REPO / "results_organized/genome_clustering/project/genome_matrix"
PATH_X = PROJ_GENO / "genome_feature_matrix.npz"
PATH_LAB = PROJ_GENO / "genome_feature_labels.npz"

PATH_CLASSIF = (
    REPO / "results_organized/genome_clustering/project/notebook7_extended"
    / "ie_instance_classification.tsv"
)
PATH_TRUTH = (
    REPO / "results_organized/genome_clustering/project/notebook7_extended"
    / "integrase_truth_table_stage9.tsv"
)
# Not used as a hard filter (6960 unknown is the cohort); kept for API parity.
PATH_STAGE9_FINAL = FF / "artifacts/stage9_FINAL.tsv"
PATH_COHORT_UNKNOWN = (
    FF / "artifacts/attL_complement_blast_window"
    / "cohort_6960_atlas_classes_conj_split.tsv"
)
PATH_COHORT_JOIN = (
    FF / "artifacts/attL_complement_blast_window" / "cohort_6960_atlas_join.tsv"
)
# attL–attR-only fragments (~100–300 bp, integrase outside window) — not valid IE
PATH_EXCLUDE_ATT_ONLY = ART / "exclude_att_only_short_ie.txt"
PATH_GTDB = REPO / "results_organized/misc/gtdb_taxonomy.tsv"
PATH_SAMPLES_META = REPO / "results_organized/summary/samples_metadata.tsv"

UN_BASE = FF / "03_clustering" / "reference"
PATH_PFAM_DESC = UN_BASE / "pfam_descriptions.tsv"
PATH_LEIDEN_OLD = UN_BASE / "filtered_ie_leiden_hd.tsv"
PATH_LEIDEN_S8 = FF / "artifacts/nb04_cluster_unclassified/leiden_clusters_stage8_unc.tsv"
PATH_LEIDEN_S9 = FF / "artifacts/nb04_cluster_unclassified_stage9/leiden_clusters_stage9_unc.tsv"

GRAPH_K = 15
HDBSCAN_MIN_CLUSTER_SIZE = 5
HDBSCAN_MIN_SAMPLES = 3
LEIDEN_SEED = 42
# Leiden communities smaller than this are collapsed to noise (-1)
MIN_LEIDEN_CLUSTER_SIZE = 5


def ensure_dirs() -> None:
    ART.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    IDS.mkdir(parents=True, exist_ok=True)


def art(name: str) -> Path:
    return ART / name.replace("stage8_unc", TAG).replace("stage9_unc", TAG).replace(
        "c6960_unk", TAG
    )


def fig(name: str) -> Path:
    return FIG / name
