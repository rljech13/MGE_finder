"""Classifier behaviour on a locus, without GBFF or eggNOG tables."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "atlas_pipeline"))

from ie_classifier.classify import classify_from_genes
from ie_classifier.config import load_config
from ie_classifier.rules import coarse_label_from_primary

LABELS = ROOT / "atlas_pipeline" / "ie_classifier" / "ie_atlas_labels.yaml"


def _gene(i, pfam, product):
    return {
        "locus_tag": f"g{i}",
        "start": i * 100,
        "end": i * 100 + 90,
        "pfam_id": pfam,
        "product": product,
    }


class TestCoarseLabel(unittest.TestCase):
    def test_published_buckets(self):
        self.assertEqual(coarse_label_from_primary("phage.caudoviral.complete"), "phage_caudoviral")
        self.assertEqual(coarse_label_from_primary("conj.T4SS_ICE"), "conj")
        self.assertEqual(coarse_label_from_primary("rep.helicase_only"), "rep")
        self.assertEqual(coarse_label_from_primary("other.unclassified_IE"), "unclassified")
        self.assertEqual(coarse_label_from_primary("island.RM_module"), "island_module")


class TestClassifyFromGenes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(LABELS)
        cls.meta = {
            "element_id": "e1",
            "genome_id": "g",
            "tRNA_type_assigned": "tRNA-Leu",
            "integrase_id_lookup": "i1",
        }
        cls.att = {"attL_attR_present": True}

    def test_integrase_alone_stays_unclassified(self):
        rec = classify_from_genes([_gene(1, "PF00589", "integrase")], self.meta, self.cfg, self.att)
        self.assertEqual(rec["primary_label"], "other.unclassified_IE")
        self.assertEqual(rec["coarse_label"], "unclassified")
        self.assertGreater(rec["IE_core_score"], 0)

    def test_integrase_plus_t4ss_is_conjugative(self):
        genes = [
            _gene(1, "PF00589", "integrase"),
            _gene(2, "PF03432", "relaxase"),
            _gene(3, "PF12696", "VirD4"),
            _gene(4, "PF00437", "VirB4"),
        ]
        rec = classify_from_genes(genes, self.meta, self.cfg, self.att)
        self.assertEqual(rec["primary_label"], "conj.T4SS_ICE")
        self.assertEqual(rec["coarse_label"], "conj")
        self.assertEqual(rec["conjugation_score"], 1.0)

    def test_empty_locus_is_unclassified(self):
        rec = classify_from_genes([], self.meta, self.cfg, self.att)
        self.assertEqual(rec["primary_label"], "other.unclassified_IE")
        self.assertEqual(rec["n_cds"], 0)


if __name__ == "__main__":
    unittest.main()
