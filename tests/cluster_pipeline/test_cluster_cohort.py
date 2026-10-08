"""Clustering parameters and inputs come from the caller, not from a cohort name."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cluster_pipeline"))

from cluster_cohort import cluster_label, leiden_knn


def _two_groups() -> csr_matrix:
    # Six elements in group A share features 0-1. Six in group B share features 2-3.
    rows = []
    for _ in range(6):
        rows.append([1, 1, 0, 0])
    for _ in range(6):
        rows.append([0, 0, 1, 1])
    return csr_matrix(np.array(rows, dtype=np.float32))


class TestLeidenOnSyntheticMatrix(unittest.TestCase):
    def test_two_disjoint_groups_stay_apart(self):
        from cluster_cohort import jaccard
        labels = leiden_knn(jaccard(_two_groups()), k=3, min_cluster_size=5, seed=1)
        self.assertEqual(set(np.unique(labels[:6])), set(np.unique(labels[:6])))
        self.assertEqual(len(set(labels[:6])), 1)
        self.assertEqual(len(set(labels[6:])), 1)
        self.assertNotEqual(labels[0], labels[6])
        self.assertTrue((labels >= 0).all())

    def test_tiny_community_becomes_noise(self):
        sim = np.eye(3, dtype=np.float32)
        labels = leiden_knn(sim, k=1, min_cluster_size=5, seed=1)
        self.assertTrue((labels < 0).all())

    def test_other_labels_are_left_out(self):
        classification = pd.DataFrame({
            "integrase_id": [f"i{i}" for i in range(12)],
            "genome_id": [f"g{i}" for i in range(12)],
            "coarse_label": ["unclassified"] * 6 + ["conj"] * 6,
        })
        genomes = np.array([f"g{i}" for i in range(12)])
        out = cluster_label(
            classification, _two_groups(), genomes,
            keep_label="unclassified",
            label_column="coarse_label",
            id_column="integrase_id",
            genome_column="genome_id",
            k=3, min_cluster_size=5, seed=1,
        )
        self.assertEqual(len(out), 6)
        self.assertTrue(out["integrase_id"].isin([f"i{i}" for i in range(6)]).all())

    def test_genome_missing_from_matrix_is_an_error(self):
        classification = pd.DataFrame({
            "integrase_id": ["i0"],
            "genome_id": ["missing"],
            "coarse_label": ["unclassified"],
        })
        out = cluster_label(
            classification, _two_groups(), np.array(["g0"]),
            keep_label="unclassified",
            label_column="coarse_label",
            id_column="integrase_id",
            genome_column="genome_id",
            k=3, min_cluster_size=5, seed=1,
        )
        self.assertEqual(out.loc[0, "error"], "genome_not_in_matrix")

    def test_script_has_no_manuscript_cohort(self):
        text = (ROOT / "cluster_pipeline" / "cluster_cohort.py").read_text()
        self.assertNotIn("/home/", text)
        self.assertNotIn("6960", text)
        self.assertNotIn("c6960", text)


class TestSaveRoundtrip(unittest.TestCase):
    def test_cli(self):
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            X = _two_groups()
            save_npz(tmp / "X.npz", X)
            np.savez(tmp / "genomes.npz", genomes=np.array([f"g{i}" for i in range(12)]))
            pd.DataFrame({
                "integrase_id": [f"i{i}" for i in range(12)],
                "genome_id": [f"g{i}" for i in range(12)],
                "coarse_label": ["unclassified"] * 12,
            }).to_csv(tmp / "cls.tsv", sep="\t", index=False)
            out = tmp / "leiden.tsv"
            proc = subprocess.run(
                [
                    sys.executable, str(ROOT / "cluster_pipeline" / "cluster_cohort.py"),
                    "--classification", str(tmp / "cls.tsv"),
                    "--matrix", str(tmp / "X.npz"),
                    "--genomes", str(tmp / "genomes.npz"),
                    "--k", "3",
                    "--min-cluster-size", "5",
                    "--out", str(out),
                ],
                capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            got = pd.read_csv(out, sep="\t")
            self.assertEqual(got["leiden_cluster"].nunique(), 2)


if __name__ == "__main__":
    unittest.main()
