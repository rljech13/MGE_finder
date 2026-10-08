# 3. Clusters of unknown elements

The atlas class `unknown` in the manuscript is compared by protein content. The network and domain-role code is `cluster_pipeline/lib/`. `lib/` here is a link to that package.

A cohort you choose is clustered with `cluster_pipeline/cluster_cohort.py`. The script below rebuilds the manuscript figures for the frozen unknown set. It is not the general command.

```bash
python final_final/03_clustering/run_cohort6960_unknown_graph_pipeline.py
```

That run reads `artifacts/attL_complement_blast_window/cohort_6960_atlas_classes_conj_split.tsv` and `results_organized/genome_clustering/project/genome_matrix/genome_feature_matrix.npz`.
