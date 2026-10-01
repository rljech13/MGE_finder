# 02. Small RM and small TA

Same cargo as experiment 01, kept only when the element stays small once replication genes are not counted as extra cargo.

RM (`scan_small_rm_allow_rep.py`): CDS that are neither integrase nor replication ≤ 3. The RM genes themselves count toward the three.

TA (`scan_small_ta_allow_rep.py`): TA genes plus other non-integrase, non-replication CDS ≤ 3. Unknown atlas class only.

```bash
python ../run.py run 02_small_rm_ta
```

| step | script | reads | writes |
|---|---|---|---|
| 1 | `scan_small_rm_allow_rep.py` | same cohort tables as experiment 01 | `artifacts/unknown_le3_rm/rm_allow_rep.tsv` |
| 2 | `scan_small_ta_allow_rep.py` | same | `artifacts/unknown_le3_ta/ta_allow_rep.tsv`, `ta_allow_rep_genes.tsv` |
| 3 | `plot_small_rm_histograms.py` | `rm_allow_rep.tsv`, atlas join, GTDB taxonomy | `figures/unknown_le3_rm/` |
| 4 | `plot_small_ta_histograms.py` | `ta_allow_rep.tsv`, `ta_allow_rep_genes.tsv`, atlas join, taxonomy | `figures/unknown_le3_ta/` |

`small_rm_loci.tsv` and `unknown_le3_extra/catalog.tsv`, if present, only narrow the candidate list. The scans run without them.
