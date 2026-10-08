# Clustering

Groups one class of a classification table by shared proteins. Jaccard similarity, a k-nearest-neighbour graph, then Leiden. The matrix, the genome names, the label to keep, `k`, the minimum cluster size, and the seed are arguments.

```bash
python cluster_pipeline/cluster_cohort.py \
    --classification classification.tsv \
    --matrix genome_feature_matrix.npz \
    --genomes genome_feature_labels.npz \
    --keep-label unclassified \
    --k 15 \
    --min-cluster-size 5 \
    --seed 42 \
    --out leiden.tsv
```

`--genomes` is an npz with an array named `genomes`, in the same order as the matrix rows. `genome_id` in the classification table has to match those names. Elements with another label are left out.

`final_final/03_clustering/run_cohort6960_unknown_graph_pipeline.py` rebuilds the manuscript figures for that one cohort. It is not a general clustering command.

`protein_clusterization/` and `classification_investigation/clusterMGE/` are different clusterings.
