# 4. Cargo on unknown elements

What already occurs as cargo next to the integrase on unknown elements. This chapter does not decide whether modules must combine. It is an inventory.

Three steps. Step 2 reads the tables from step 1. Step 3 reads both.

```bash
cd final_final/04_cargo
python run.py check
python run.py run 01_rm_ta_cohort
```

| step | directory | question | writes |
|---|---|---|---|
| 1 | `01_rm_ta_cohort/` | Which unknown elements carry RM or TA. Only the 3,738 that have PADLOC rows. | `artifacts/unknown_rm/`, `artifacts/unknown_ta/` |
| 2 | `02_small_rm_ta/` | Which of those stay small once replication is not counted as extra cargo. | `artifacts/unknown_le3_rm/`, `artifacts/unknown_le3_ta/` |
| 3 | `03_class_inventory/` | Which literature classes occur on all 4,689 unknown elements, and where the class is exactly one. | `artifacts/unknown_ime_classes/` |

Each step has its own `README.md`: script, files read, files written. The shared cohort is the tables in `artifacts/attL_complement_blast_window/` and the representative GBFF in `results_organized_rep_6960/`.

Chapter 5 reads `artifacts/unknown_ime_classes/after_hypothesis/unknown_minimal_ie.tsv`, which appears at step 3.
