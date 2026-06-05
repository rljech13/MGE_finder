"""Unit tests for V_3pS strict and confident IE filter logic."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "IE_finder" / "scripts"))

from v3ps_filters import (
    FilterThresholds,
    classify_prodigal_overlap,
    evaluate_ie_candidate,
    select_v3ps_strict_hit,
)


class TestV3psStrict(unittest.TestCase):
    def test_select_longest_strict_hit(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 68, "sstart": 100, "send": 200,
             "length": 68, "wstart": 1000},
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 69, "sstart": 100, "send": 220,
             "length": 85, "wstart": 1000},
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 50, "sstart": 100, "send": 150,
             "length": 50, "wstart": 1000},
        ])
        hit = select_v3ps_strict_hit(raw, trna_len=70, trna_strand="+", shift=3)
        self.assertIsNotNone(hit)
        self.assertEqual(hit["best_attl_len_bp"], 85)

    def test_rejects_short_3p_anchor(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 60, "sstart": 100, "send": 200,
             "length": 60, "wstart": 1000},
        ])
        self.assertIsNone(select_v3ps_strict_hit(raw, trna_len=70, trna_strand="+", shift=3))

    def test_prodigal_intergenic(self):
        info = classify_prodigal_overlap(500, 520, [(600, 700, "+", "cds1")])
        self.assertEqual(info["prodigal_hit_class"], "intergenic")

    def test_prodigal_partial_overlap(self):
        info = classify_prodigal_overlap(500, 520, [(510, 700, "+", "cds1")])
        self.assertEqual(info["prodigal_hit_class"], "partial_CDS")

    def test_confident_passes_all_filters(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 69, "sstart": 100, "send": 220,
             "length": 20, "wstart": 1000},
        ])
        cds = {"c1": [(600, 700, "+", "cds1")]}
        audit = evaluate_ie_candidate(
            integrase_id="i1",
            contig="c1",
            trna_start=1,
            trna_end=70,
            trna_strand="+",
            trna_len=70,
            integrase_start=1000,
            integrase_end=1902,
            ie_id="i1:c1:100-3500",
            ie_len_nt=2500,
            ie_has_n=False,
            raw_hits=raw,
            cds_by_contig=cds,
            thresholds=FilterThresholds(),
        )
        self.assertTrue(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "")

    def test_confident_rejects_short_integrase(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 69, "sstart": 100, "send": 220,
             "length": 20, "wstart": 1000},
        ])
        audit = evaluate_ie_candidate(
            integrase_id="i1",
            contig="c1",
            trna_start=1,
            trna_end=70,
            trna_strand="+",
            trna_len=70,
            integrase_start=1000,
            integrase_end=1500,
            ie_id="i1:c1:100-3500",
            ie_len_nt=2500,
            ie_has_n=False,
            raw_hits=raw,
            cds_by_contig={"c1": [(2000, 3000, "+", "cds1")]},
            thresholds=FilterThresholds(),
        )
        self.assertFalse(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "integrase_too_short")


if __name__ == "__main__":
    unittest.main()
