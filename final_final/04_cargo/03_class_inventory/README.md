# 03. Cargo classes and minimal modules

Which literature classes sit on the 4,689 unknown elements, and which elements are the integrase plus exactly one class.

This experiment reads the RM/TA scan from 01 and the small-RM / small-TA tables from 02.

```bash
python ../run.py run 03_class_inventory
```

| step | script | reads | writes |
|---|---|---|---|
| 1 | `count_unknown_minimal_ie.py` | cohort class table, `unknown_rm_ta_scan.tsv` from 01, PADLOC, representative GBFF. `unknown_metals_markers.py` is the metals vocabulary | `artifacts/unknown_ime_classes/unknown_minimal_ie.tsv`, `ime_class_funnel.tsv` |
| 2 | `score_unknown_module_minimal.py` | step 1, `proteins.tsv`, `rm_allow_rep.tsv` and `ta_allow_rep.tsv` from 02 | module-minimal table in the same directory |
| 3 | `recount_unknown_minimal_annotated.py` | step 1, `unknown_all_annotated.tsv`, PADLOC | recount after Pfam names on dark CDS |
| 4 | `recount_unknown_minimal_after_hypothesis.py` | step 1 plus the still-dark stack (`still_dark_calls.tsv`, `leftover_calls.tsv`, `hypothesis_calls.tsv`, `cis_elements.tsv`) | `after_hypothesis/unknown_minimal_ie.tsv`, `by_class.tsv` |
| 5 | `apply_duf_as_unannotated.py` | `after_hypothesis/unknown_minimal_le1_long.tsv`, `unknown_all_annotated.tsv` | updates that table in place. A DUF is obligatory only at ≥ 70% of elements that carry the class |
| 6 | `inventory_unknown_funnel_classes.py` | `after_hypothesis/unknown_minimal_ie.tsv`, `by_class.tsv` | `class_inventory_summary.tsv` |
| 7–9 | `plot_unknown_module_minimal.py`, `plot_unknown_minimal_diversity.py`, `plot_unknown_minimal_histograms.py` | the tables from steps 4 and 5, plus the atlas join | `figures/unknown_ime_classes/` |

Annotation tables are under `final_final/artifacts/unknown_reannotate/`.
