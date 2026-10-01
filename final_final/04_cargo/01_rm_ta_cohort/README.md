# 01. RM and TA on the unknown cohort

Which unknown elements carry restriction-modification or a toxin-antitoxin pair. No size cut. The 951 unknown elements added later have no PADLOC rows, so this experiment stays on the 3,738 that do.

```bash
python ../run.py run 01_rm_ta_cohort
```

| step | script | reads | writes |
|---|---|---|---|
| 1 | `scan_unknown_rm_ta.py` | `cohort_6960_atlas_join.tsv`, `cohort_6960_four_func_classes.tsv`, `atlas_rep_index.tsv`, `gtdb_taxonomy.tsv`, `padloc_per_ie_stage8_unc.tsv`, representative GBFF | `artifacts/unknown_rm/unknown_rm_ta_scan.tsv`, `unknown_rm_hits.tsv`, `artifacts/unknown_ta/unknown_ta_hits.tsv` |
| 2 | `plot_unknown_rm_ta_diversity.py` | the two hit tables from step 1 | `figures/unknown_rm/`, `figures/unknown_ta/` |

The cohort tables are in `final_final/artifacts/attL_complement_blast_window/`. The representative sequences are in `results_organized_rep_6960/`.
