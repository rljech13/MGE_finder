# Atlas

Assigns each integrative element a class from the markers in its locus. The rules and the marker list are `ie_classifier/ie_atlas_labels.yaml`.

The cohort is yours. Paths to the truth table, the representative index, the locus annotations, and the GenBank files are arguments. Nothing in this command points at the manuscript tables.

```bash
python atlas_pipeline/classify_cohort.py \
    --cohort my_ids.tsv \
    --truth integrase_truth.tsv \
    --rep-index rep_index.tsv \
    --locus-table locus_eggnog_merged.tsv \
    --gbff-dir path/to/gbff \
    --out classification.tsv
```

`my_ids.tsv` has a column `integrase_id`. To try another marker list, pass `--labels`.

`final_final/02_atlas/build_stage5_atlas.py` is the manuscript run. It reads the Stage 5 table already on disk and writes back into `final_final/artifacts/`. It is not the way to classify a new cohort. Its truth table comes only from `MGE_INTEGRAS_TRUTH`.
