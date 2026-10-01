# 2. Atlas

After the cohort, each element receives a class from the markers in its locus. The classifier is `ie_classifier/` in this directory. Thresholds and the marker list are `ie_classifier/ie_atlas_labels.yaml`. The rule logic in that file is unchanged. The agreed four primary classes (phage, conjugative, phage_like, unknown) are described in `docs/corpus/IE_ATLAS_FOUR_CLASS.md`.

The old paths `results_organized/scripts/ie_classifier` and `ie_atlas_labels.yaml` under the notebook directory are links to this package.

| script | what it labels | where it writes |
|---|---|---|
| `build_stage5_atlas.py` | Stage 5 elements missing from the frozen 25,591-row atlas | `artifacts/stage5_atlas_classification.tsv`, `stage5_atlas_join.tsv` |
| `build_blast_window_40kb_atlas.py` | The same rules on a window of up to 40 kb from the tRNA | `artifacts/blast_window_40kb_atlas_classification.tsv` |
| `build_recombinase_atlas.py` | Serine recombinase cohort, CDS from Bakta | `artifacts/recombinase_atlas/` |

```bash
cd final_final/02_atlas
python build_stage5_atlas.py
```

Inputs are GBFF and eggNOG tables from `results_organized/`. The stage 5 truth table lives outside this repository. Point `MGE_INTEGRAS_TRUTH` at that TSV.
