# 3. Clusters of unknown elements

The atlas leaves a large unknown class: those elements do not pass the phage, conjugative, or phage_like rules. Here they are compared by protein content, not by marker weights.

The matrix is the same one used for the whole MGE collection: `results_organized/genome_clustering/project/genome_matrix/genome_feature_matrix.npz`. Jaccard is computed on the cohort rows, a 15-nearest-neighbour graph is built, then Leiden and, for comparison, HDBSCAN. Pfam markers, the network between clusters, and labels on that graph follow.

The element list is `artifacts/attL_complement_blast_window/cohort_6960_atlas_classes_conj_split.tsv`, class unknown. Shared paths are in `cohort6960_unknown_paths.py`.

```bash
cd final_final/03_clustering
python run_cohort6960_unknown_graph_pipeline.py
python run_cohort6960_unknown_graph_pipeline.py --from-step 3
```

| step | script | what it does |
|---|---|---|
| 1 | `build_cohort6960_unknown_clustering.py` | Submatrix, kNN, Leiden, HDBSCAN |
| 2 | `recompute_hallmark_hmm_imputed_c6960_unk.py` | Domain roles; an integrase that passed the discovery HMM is kept |
| 3 | `plot_integrase_diversity_c6960_unk.py` | Integrase subtypes by cluster |
| 4 | `render_cluster_markers_table_c6960_unk.py` | Marker table |
| 5 | `render_cluster_markers_heatmap_c6960_unk.py` | Domain heatmap |
| 6 | `build_cluster_network_c6960_unk.py` | Network of links between clusters |
| 7 | `enrich_cluster_overview_c6960_unk.py` | Biology labels on the network nodes |
| 8 | `enrich_c6960_unk_ie_knn_graphml.py` | Host and atlas class on the element graph |

`lib/` holds the network and domain-role functions that used to live in the old `final/` tree. `reference/` holds Pfam descriptions and the previous Leiden labelling, which the new run is compared against. The `MGE_finder` environment is required (`igraph`, `leidenalg`, `hdbscan`).
