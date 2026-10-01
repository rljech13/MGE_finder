# 1. Confident IE discovery

This is the first part of the paper: from a genome to the cohort used by the atlas, the clusters, and the modules. The code is the Snakemake workflow in `finder_pipeline/`, not in this directory.

The order in `finder_pipeline/scripts/v3ps_filters.py` and `finder_pipeline/README.md` is:

1. Phage-type tyrosine integrase, profiles PF00589 and PF22022.
2. A tRNA on the opposite strand, 3′ ends no more than 500 bp apart.
3. An exact direct repeat of that tRNA 3′ end, at least 8 bp, in a window of up to 300 kb upstream of the integrase.
4. Integrase longer than 300 aa. Length is `(end − start + 1) // 3`. The cut is exclusive: greater than 300.
5. attL does not overlap a CDS. The pipeline uses Prodigal and keeps the class `intergenic`.
6. Site length: an exact repeat of at least 14 bp, or a BLAST hit of at least 17 bp. Mismatches are allowed in that second case. There is no separate 10% cap.
7. Alignment without a coordinate gap: the hit length equals both the query span and the subject span.
8. MMseqs2 deduplication at identity 0.9 and coverage 0.8. The frozen cohort after this step is 6,960 elements.

How to run the pipeline is in `finder_pipeline/README.md`. The attL search window defaults to 300 kb (`filters.attl_window_bp`).
