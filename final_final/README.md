# Analysis order in the paper

The numbers follow the order in which the analysis is told. Each chapter's code lives in the directory with the same number. Tables and figures are still written under `artifacts/` and `figures/` in this tree. Those outputs are not steps.

| chapter | directory | what it answers |
|---|---|---|
| 1 | `01_discovery/` | How genomes become the confident IE cohort. The pipeline itself is `finder_pipeline/`. |
| 2 | `02_atlas/` | How an element is labelled. The package is `atlas_pipeline/ie_classifier/`. |
| 3 | `03_clustering/` | How unknown elements group. The library is `cluster_pipeline/lib/`. |
| 4 | `04_cargo/` | Which functions sit on unknown elements: RM, TA, the other classes, and a minimal module. |
| 5 | `05_abundance/` | Whether a class is more common inside the element than in the rest of the same genome. |

Chapter 2 reads loci and annotations from `results_organized/`. Chapter 3 reads the chapter 2 labels and the feature matrix in `results_organized/genome_clustering/project/genome_matrix/`. Chapter 4 reads the cohort and, for the class inventory, the RM/TA tables it writes in its own first step. Chapter 5 reads the class table that chapter 4 writes to `artifacts/unknown_ime_classes/after_hypothesis/`.
