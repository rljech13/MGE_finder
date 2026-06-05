"""Unit tests for IE deduplication helpers (no MMseqs required)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "finder_pipeline" / "scripts"))

from dedup_ie_representatives import cluster_id, resolve_rep


class TestDedupHelpers(unittest.TestCase):
    def test_cluster_id_format(self):
        cid = cluster_id("orf1:ctg:1-100", "GCA_test")
        self.assertEqual(cid, "orf1:ctg:1-100|sample=GCA_test")

    def test_resolve_rep_via_member_map(self):
        cl = pd.DataFrame([["1", "1"], ["2", "1"]], columns=["member_id", "rep_id"])
        member_to_rep = {"1": "1", "2": "1"}
        self.assertEqual(resolve_rep("2", member_to_rep, {"1"}, cl), "1")


if __name__ == "__main__":
    unittest.main()
