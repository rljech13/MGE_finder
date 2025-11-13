TODO: complete readme

Overview

This dir contains an end-to-end Snakemake workflow (finder_pipeline) that screens bacterial / archaeal assemblies for mobile genetic elements (MGEs) integrating into tRNA genes.
This workflow starts with raw genome FASTA files, predicts ORFs, screens them with custom Pfam HMMs for integrases, finds the nearest tRNA within 500 bp, carves out the putative MGE region, and finally produces an annotated GenBank file with attL/attR sites.

---

finder_pipeline layout

finder_pipeline/
├── Snakefile                 # workflow definition
├── scripts/                  # helper scripts 
│   ├── prepare_fastas.py
│   ├── predict_orfs.py
│   ├── hmm_search.py
│   ├── annotate_trna_proximity.py
│   ├── extract_trna_region.py
│   ├── annotate_mge_region.py
│   ├── extract_mge_regions.py
│   └── annotate_and_orient_mge.py
├── finder_config.yaml               
└── finder_results/                  # output will be created here (git-ignored)

logger.py lives at finder_pipeline/scripts/logger.py and is imported by every script; 
⸻

Quick start

COMMING SOON! 

Output files: 

results/⟨sample⟩/ will now hold:
	•	integrase_hits_summary.tsv – full table of candidate integrases
	•	integrase_trna.tsv – integrase ↔ closest tRNA pairs
	•	mge_region.fa – extracted DNA of each putative element
	•	mge_annotated.gbk – GenBank with attL/attR + integrase CDS
	•	attachment_sites.tsv – table of computed att coordinates

Config:
TODO: carefully recreate config, change the way MGE_finder is intiated, provide not absolute, but relative links

Configuration (finder_config.yaml)

paths:
  input_sources:                # list of dirs searched recursively
    - "/data/(provide your genomes here!)"

  genomes_dir: "/scratch/genomes_fna"
  results_dir: "/scratch/mge_results"
pfam_profiles:
  - "/db/pfam/PF00589.hmm"      # example
  - "/db/pfam/PF13408.hmm"
execution:
  conda_env: "MGE_finder"       # the env created above


⸻

Snakemake workflow

Rule	What it does	Outputs
prepare_fasta	Converts/collects raw assemblies into canonical *.fna (NCBI & Hybracter formats supported).	genomes_dir/*.fna + .complete flag
predict_orfs	Runs Prodigal in single-genome mode, produces ORF GFF + nucleotide/protein FASTA.	orfs.gff, orfs.ffn, orfs.faa
build_combined_hmm	Concatenates user-listed Pfam HMMs and runs hmmpress.	pfam_combined.hmm (+ .h3*)
hmm_search	hmmscan each proteome against combined HMM; summarises integrase coordinates.	integrase_hits_summary.tsv, integrase_orfs.tsv
predict_trna	Calls Aragorn to find genomic tRNAs.	trna.tsv
trna_proximity	Matches integrase ↔ nearest opposite-strand tRNA ≤ 500 bp.	integrase_trna.tsv
extract_trna_region	Takes only the closest tRNA per integrase and writes its FASTA sequence (needed later as BLAST query).	mge_query.fa
blast_mge	BLASTs each extracted tRNA against a ±100 kb window around its integrase; retains partial matches starting at tRNA 5′ end (attL candidates).	mge_blast.tsv
extract_mge_region	From BLAST + tRNA table derives final region coordinates and saves DNA.	mge_region.fa
annotate_mge	Orientates region, reverse-complements when needed, writes GenBank with attL/attR features and integrase CDS.	mge_annotated.gbk, attachment_sites.tsv




Script docs: 

prepare_fastas.py

Function	Purpose
find_files(base_dir, pattern)	Recursive glob helper supporting *_genomic.fna (NCBI) and barcode*.fastq_final.fasta (Hybracter).
process_input_sources(sources, out_dir, …)	Converts / copies all matching files into out_dir/⟨sample⟩.fna, using BioPython for FASTQ→FASTA conversion.
main(config_path, output_done)	Loads YAML, calls the above and writes an optional .complete sentinel for Snakemake.

predict_orfs.py

Thin wrapper around Prodigal single-genome mode.

Function	Purpose
predict_with_prodigal()	Builds and runs the CLI command; logs all paths.
main()	CLI entry-point / Snakemake shim.

hmm_search.py

Function	Purpose
run_hmmscan()	Executes hmmscan --tblout against combined HMM.
parse_tblout()	Extracts ORF ID → Pfam accession mapping.
parse_faa()	Reads Prodigal FAA headers (>contig_id_# # start # end # strand) and stores coordinates.
parse_gff()	Pulls contig lengths from the Prodigal GFF meta-lines.
write_outputs()	Creates integrase_hits_summary.tsv (full coordinates + contig length) and integrase_orfs.tsv (short list).

annotate_trna_proximity.py

Function	Purpose
parse_integrases()	Reads summary TSV, casts numeric columns.
parse_trna()	Parses raw Aragorn text output; returns contig / coords / strand / type.
effective_coord_integrase() / effective_coord_trna()	Define where to measure distance on each strand (integrase 3′ vs tRNA 3′).
find_nearby_trnas()	Opposite-strand tRNA within ≤ 500 bp; chooses all candidates.
write_results()	Outputs integrase_trna.tsv.

extract_trna_region.py

Extracts only the closest tRNA per integrase (smallest distance) and writes a FASTA record:

>INT123:contigX:12345-12409:+
ACGT…

Used later as BLAST query.

annotate_mge_region.py  (key filtering step updated)

Function	Purpose
extract_subject_region()	Pulls a ±100 kb window around the tRNA/integrase pair.
run_blast_on_region()	Makes an in-memory BLAST DB for that window, runs blastn -word_size 4 -dust no. Now filters for hits where qstart == 1 and length < len(tRNA) – ensuring attL contains the tRNA 3′ end exactly as discussed.
main()	Iterates every tRNA query, records best hit, accumulates into mge_blast.tsv.

extract_mge_regions.py

Function	Purpose
parse_blast()	Keeps best bitscore per integrase.
parse_trna()	Loads tRNA proximity table.
filter_closest_trna()	Again ensures single tRNA per integrase (extra safety).
build_region_data()	Calculates final element coordinates: if tRNA +, region = [tRNA_start, BLAST_end]; if tRNA –, reverse.
extract_regions()	Writes DNA subsequences → mge_region.fa.

annotate_and_orient_mge.py

Creates the final GenBank:
	•	Reverse-complements region if tRNA on + to keep attL at 5′, attR at 3′
	•	Adds three features: attL, attR, integrase CDS
	•	attachment_sites.tsv summarises their coordinates.

⸻

Customising the distance / filters
	•	Distance threshold – change --max_distance in Snakemake (trna_proximity rule) or via CLI.
	•	Window size – tweak WINDOW_SIZE (100 000 bp default) atop annotate_mge_region.py.
	•	BLAST filters – in the patch above use qstart <= N or add qend >= len(tRNA) - M if you also want perfect 3′ overlap.

⸻

Re-generating combined HMM

If you add/remove Pfam profiles, simply edit config.yaml and rerun:

snakemake build_combined_hmm

All downstream rules will be auto-invalidated.

