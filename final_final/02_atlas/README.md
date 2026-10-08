# 2. Atlas

Each element is labelled from the markers in its locus. The rule package is `atlas_pipeline/ie_classifier/`. `ie_classifier/` in this directory is a link to that package.

A new cohort is classified with `atlas_pipeline/classify_cohort.py`. The scripts in this directory are the manuscript runs: they read `results_organized/` and write `final_final/artifacts/`.

| script | what it labels |
|---|---|
| `build_stage5_atlas.py` | Stage 5 elements missing from the frozen atlas. Truth table: `MGE_INTEGRAS_TRUTH`. |
| `build_blast_window_40kb_atlas.py` | The same rules on a window of up to 40 kb from the tRNA |
| `build_recombinase_atlas.py` | Serine recombinase cohort, CDS from Bakta |
