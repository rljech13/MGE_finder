"""The cohort classifier takes the user's tables. It has no manuscript path."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "atlas_pipeline"))

from classify_cohort import classify_cohort


class TestClassifyCohort(unittest.TestCase):
    def test_missing_truth_row_is_reported(self):
        truth = pd.DataFrame({"integrase_id_lookup": ["known"]})
        out = classify_cohort(["absent"], truth, {}, pd.DataFrame(), {}, Path("."))
        self.assertEqual(out.loc[0, "error"], "not_in_truth_table")
        self.assertEqual(out.loc[0, "coarse_label"], "")

    def test_cli_requires_the_user_paths(self):
        script = ROOT / "atlas_pipeline" / "classify_cohort.py"
        proc = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--cohort", proc.stderr)

    def test_script_has_no_machine_path(self):
        text = (ROOT / "atlas_pipeline" / "classify_cohort.py").read_text()
        self.assertNotIn("/home/", text)
        self.assertNotIn("6960", text)


class TestCliRoundtripOfMissingId(unittest.TestCase):
    def test_writes_the_missing_row(self):
        script = ROOT / "atlas_pipeline" / "classify_cohort.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "cohort.tsv").write_text("integrase_id\nnope\n")
            (tmp / "truth.tsv").write_text("integrase_id_lookup\n")
            (tmp / "rep.tsv").write_text("cluster_id\tfile\trep_id\n")
            (tmp / "locus.tsv").write_text("genome_id\tlocus_tag\n")
            (tmp / "gbff").mkdir()
            out = tmp / "out.tsv"
            proc = subprocess.run(
                [
                    sys.executable, str(script),
                    "--cohort", str(tmp / "cohort.tsv"),
                    "--truth", str(tmp / "truth.tsv"),
                    "--rep-index", str(tmp / "rep.tsv"),
                    "--locus-table", str(tmp / "locus.tsv"),
                    "--gbff-dir", str(tmp / "gbff"),
                    "--out", str(out),
                ],
                capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            got = pd.read_csv(out, sep="\t", dtype=str)
            self.assertEqual(got.loc[0, "error"], "not_in_truth_table")


if __name__ == "__main__":
    unittest.main()
