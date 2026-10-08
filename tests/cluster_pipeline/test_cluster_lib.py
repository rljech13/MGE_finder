"""Pure pieces of the unknown-element clustering, without the feature matrix."""
import sys
import unittest
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cluster_pipeline" / "lib"))

from cluster_biology import assign_role
from clustering_ie_graph_analysis import _bh_fdr, _compute_pair_stats, label_map_from_graph
from enrich_filtered_ie_knn_graphml import normalize_cluster


class TestDomainRole(unittest.TestCase):
    def test_manual_map_beats_the_name(self):
        self.assertEqual(assign_role("PF00589", full_name="something else"), "integrase_tyrosine")

    def test_clan_heuristic(self):
        self.assertEqual(assign_role("PF00000", clan="TPR"), "tpr_repeat")

    def test_name_heuristic_then_other(self):
        self.assertEqual(assign_role("PF00000", full_name="DEAD-box helicase"), "helicase")
        self.assertEqual(assign_role("PF00000", full_name="hypothetical protein"), "other")


class TestPairStats(unittest.TestCase):
    def test_bh_fdr_is_monotone_in_rank(self):
        adj = _bh_fdr(np.array([0.01, 0.04, 0.03]))
        np.testing.assert_allclose(adj, [0.03, 0.04, 0.04])

    def test_empty_pvalues(self):
        self.assertEqual(len(_bh_fdr(np.array([]))), 0)

    def test_observed_above_mixing_null(self):
        stats = _compute_pair_stats(O=10, wsum=1.0, cap=100, E=2.0)
        self.assertAlmostEqual(stats["cross_cluster_edge_density"], 0.1)
        self.assertAlmostEqual(stats["relatedness_z"], (10 - 2) / np.sqrt(2))
        self.assertGreater(stats["fold_enrichment_vs_mixing_null"], 1)


class TestGraphLabels(unittest.TestCase):
    def test_negative_and_blank_clusters_are_dropped(self):
        g = nx.Graph()
        g.add_node("a", leiden_hd_cluster=3)
        g.add_node("b", leiden_hd_cluster=-1)
        g.add_node("c", leiden_hd_cluster="")
        self.assertEqual(label_map_from_graph(g), {"a": 3})

    def test_normalize_cluster_strips_sample_suffix(self):
        self.assertEqual(normalize_cluster("12|sample=GCA_1"), "12")
        self.assertEqual(normalize_cluster(" 7 "), "7")


if __name__ == "__main__":
    unittest.main()
