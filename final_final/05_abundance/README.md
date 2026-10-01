# 5. Enrichment against the rest of the genome

The last step of the story. For a class from chapter 4, compare how often it occurs inside the unknown element with how often it occurs in the same genome outside that element.

Reads `artifacts/unknown_ime_classes/after_hypothesis/unknown_minimal_ie.tsv` and the representative index `results_organized_rep_6960/atlas_rep_index.tsv`. Proteins outside the element are Prodigal ORFs from `finder_pipeline/results`. Profiles: Pfam, PADLOC, DefenseFinder, CasFinder. HMMER comes from `PATH`, or from `MGE_HMMSEARCH` and the matching `MGE_HMMFETCH`, `MGE_HMMPRESS`, `MGE_HMMSTAT`. Libraries can be overridden with `MGE_PFAM`, `MGE_PADLOC_HMM`, `MGE_DEFENSE_FINDER`, and `MGE_CASFINDER`.

```bash
cd final_final/05_abundance
python enrich_unknown_ie_vs_genome.py
python enrich_unknown_ie_per_genome.py
```

The first script accepts `--stage` (`hmm`, `intervals`, `faa`, `search`, `summarize`, `plot`) and `--max-genomes`. With no arguments it runs every stage and writes `after_hypothesis/genome_enrichment/`. The second script reads `enrichment_by_class.tsv` from that directory and reports a per-genome fold. Fold is `(n_IE / L_IE) / (n_rest / L_rest)`. This is a long run over genomes.
