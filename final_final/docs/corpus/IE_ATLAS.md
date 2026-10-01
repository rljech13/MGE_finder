# Corpus: IE Atlas classifier

**Status:** living notes for the modular rule-based IE classifier (not ML).  
**Canonical config:** `final_final/02_atlas/ie_classifier/ie_atlas_labels.yaml`  
**Code package:** `final_final/02_atlas/ie_classifier/`  
The old paths under `results_organized/` are links to these files.  
**Last audited:** 2026-08-08  
**Agreed four-class target (phage / conjugative / phage_like / unknown):** see [`IE_ATLAS_FOUR_CLASS.md`](IE_ATLAS_FOUR_CLASS.md).

> Goal of this file: one place that answers *what metrics exist, how they are computed, and how they become labels* — for the locked 6960 cohort and for the current YAML/code.

---

## 1. What it is

Rule-based modular classifier over a representative locus (GBFF CDS + eggNOG/Pfam annotations):

1. Match curated markers (Pfam ID and/or text regex) on each gene.
2. Aggregate **unique** marker weights into module scores ∈ [0, 1].
3. Walk `priority_order`; first label whose gate passes wins → `primary_label`.
4. Map to `coarse_label`, assign `confidence`, optional `subtype` / novel flags.

**Not used at classify time:** HMMER bit-scores/E-values, copy-number of the same marker (unique by marker *key*), Leiden clusters (diagnostic only).

---

## 2. Pipeline (one element)

```
truth row
  → load rep GBFF CDS + eggNOG hydrate (locus_tag × genome_id)
  → annotate_genes          (marker_match.py)
  → [optional] TXSScan / external detector injection
  → compute_module_scores   (scores.py)
  → assign_primary_label    (rules.py)  # + IE_core meta-boosts
  → coarse_label_from_primary
  → confidence / subtype / novel / evidence blob
```

Entry points:

| Script | Role |
|--------|------|
| `results_organized/scripts/run_ie_atlas_classification.py` | Full NB7-style run |
| `final_final/scripts/build_stage5_atlas.py` | Stage-5 universe join (source of **6960** labels) |
| `final_final/scripts/build_blast_window_40kb_atlas.py` | Classify ≤40 kb Bakta window (tRNA→IE) |
| `final_final/scripts/atlas_stage9_io.py` | Stage-9 join + palettes |

Core functions: `classify.classify_from_genes` / `classify_truth_row`.

---

## 3. Inputs

| Input | Why |
|-------|-----|
| Truth TSV (per IE) | IDs, optional `ie_len_nt`, `tRNA_type_assigned`, integrase locus |
| `rep_index.tsv` + `rep_*.gbff` | CDS products / locus tags / att notes |
| `locus_eggnog_merged.tsv` | Pfam + eggNOG text per gene |
| `ie_atlas_labels.yaml` | markers, weights, thresholds, priority |
| TXSScan / ISEScan / IntegronFinder / DefenseFinder / geNomad / VirSorter (optional) | Extra marker injection |

Matching fields (`marker_match.py`): `pfam_id`, `emapper_PFAMs`, plus text blob of `product`, `eggNOG_Description`, `seed_eggNOG_ortholog`, `KEGG_ko`, `emapper_PFAMs`.

---

## 4. Score metrics — how they work

### 4.1 Aggregation (`score_aggregation`)

```yaml
mode: marker_sum          # default
sum_scale: 1.22
ie_core_boost: 1.12
max_score_cap: 1.0
```

For a module:

\[
\text{score} = \min\!\bigl(\texttt{max\_score\_cap},\; \tfrac{\sum_{\text{unique keys}} w}{\texttt{sum\_scale}}\bigr)
\]

- Each marker **key** counts once per locus, even if present on many genes.
- Optional legacy `gene_damped`: divide by \(\sqrt{1+(n-1)\cdot\text{slope}}\times\text{divisor}\).

### 4.2 Public module scores (`compute_module_scores`)

| Score | Categories / formula |
|-------|----------------------|
| `IE_core_score` | unique `integrase` × `ie_core_boost`, then meta-boosts (below) |
| `phage_score` | packaging + head + tail + connector (**not** lysis/reg) |
| `inovirus_score` | `inovirus` |
| `conjugation_score` | `conjugation` |
| `AICE_DSS_score` | `1.0·AICE + 0.3·DSS + 0.2·repABC` |
| `replication_score` | `replication` |
| `helicase_only_score` | helicase weights if **no** pol markers; else 0 |
| `plasmid_score` | `plasmid` |
| `host_contam_score` | `host` |
| `pici_score` | `max(PICI, phage_reg aggregate)` |
| `mobilome_score` | `mobilome` |
| `is_score` | `IS` |
| `pai_score` / `rm_island_score` / `ta_island_score` | PAI / RM / max(TA_island, ta_system) |
| `island_module_score` | `min(1, pai + 0.85·rm + 0.75·ta)` |
| `competence_score` | `competence` |
| `conj_non_t4ss_score` | `conj_non_t4ss` |
| `satellite_score` / `gta_score` | satellite; `gta = max(satellite, gta_marker)` |
| `mobilization_score` | `0.9 · conjugation_score` |
| `annotation_coverage` | fraction of genes with eggNOG description or Pfam |

Internal helpers: `_phage_structural_score`, conjugation flags (relaxase/t4cp/virb4), phage structural flags (TerL, HK97 MCP, packaging-only).

### 4.3 IE_core meta-boosts (`rules._boost_ie_core` + att)

Applied **after** base score, before gates:

| Condition | Δ |
|-----------|--:|
| `tRNA_type_assigned` present | +0.15 |
| integrase locus found among CDS | +0.10 |
| attL **and** attR notes in GBFF | +0.08 |

All capped at 1.0.

**Example:** single integrase (w=0.54) → base `min(1, 0.54×1.12)=0.605`. With tRNA + locus + att → ≈0.935.

**Example caudoviral:** TerL(0.50)+HK97(0.54)=1.04 → `1.04/1.22≈0.85` phage_score.

### 4.4 Known dead channel

YAML defines **SCC** markers (`ccr`, `cch`, `cch2`, `ccpol`, `mp_partner`) and `scores.py` has synteny helpers, but **no module score aggregates category `SCC`** and **no label rule** uses them. Currently inert for classification.

---

## 5. From scores → labels

### 5.1 Priority (`priority_order`)

First match wins (conj before phage — transfer machinery outranks phage call):

1. `conj.T4SS_ICE` → `conj.IME_mobilizable` → `conj.partial_ICE_or_CIME` → `conj.non_T4SS_machinery` → `conj.competence_island` → `conj.Thermus_DSS_ICETh1_like`
2. `phage.caudoviral.complete` → `phage.inovirus_like` → `phage.PICI_like` → `phage.satellite_or_dependent` → `phage.GTA_like_or_particle_like` → `phage.cryptic_prophage_footprint`
3. `mob.*` (transposon / integron / IS-heavy / group2–casposon)
4. `island.*` (PAI / RM / TA)
5. `rep.AEP_PrimPol_helicase` → `rep.PolA_or_PolB_helicase` → `rep.helicase_only`
6. `other.unclassified_IE`

### 5.2 Key gates (`label_rules` + `anti_phage_overcall`)

| Label | Main gates |
|-------|------------|
| `phage.caudoviral.complete` | HK97 MCP **and** TerL; `phage≥0.20`; `IE_core≥0.15`; `ie_len≥10 kb`; blocked if `conjugation≥0.11` or relaxase/t4cp/virb4 |
| packaging-only TerL (no HK97) | capped → usually `phage.cryptic_prophage_footprint`, confidence ≤ medium |
| `phage.inovirus_like` | `inovirus≥0.20`, require Zot, ≥1 inov marker, `ie_len≤20 kb` |
| `phage.PICI_like` | `pici≥0.18`, integrase, Cpt **or** ≥2 `phage_reg` |
| `phage.GTA_like_or_particle_like` | `phage≥0.14` **or** `gta≥0.12`; not TerL+HK97 |
| `conj.T4SS_ICE` | relaxase **and** t4cp **and** virb4; `conjugation≥0.27`; `IE_core≥0.18` |
| `conj.IME_mobilizable` | relaxase+integrase; no strong VirB4; `conjugation≥0.21`; defer if PrimPol–helicase synteny + `replication≥0.25` |
| `conj.partial_ICE_or_CIME` | `conjugation≥0.11` + partial signature |
| `rep.AEP_PrimPol_helicase` | PrimPol + helicase within `synteny.window_genes: 8`; `replication≥0.27` |
| `rep.helicase_only` | `helicase_only≥0.17`, no pols |
| mob / island | typically scores ~0.20–0.30 + signature / integrase constraints |

### 5.3 Coarse taxonomy (`coarse_label_from_primary`)

| coarse | primaries |
|--------|-----------|
| `phage_caudoviral` | `phage.caudoviral.complete` |
| `phage_inovirus` | `phage.inovirus_like` |
| `phage_gta` | GTA / satellite / cryptic |
| `phage_pici` | `phage.PICI_like` |
| `conj` | any `conj.*` |
| `mob_element` | mob.* |
| `island_module` | island.* |
| `rep` | any `rep.*` |
| `unclassified` | `other.unclassified_IE` (+ unknown phage.*) |

### 5.4 Confidence

```yaml
high:   { min_ie_core: 0.32, min_specialized_module: 0.38 }
medium: { min_ie_core: 0.22, min_specialized_module: 0.18 }
low:    default
```

`top_spec = max(phage, inovirus, conjugation, AICE_DSS, replication, pici, mobilome, island_module, gta)`.  
Caudoviral without TerL+HK97 forced to medium. Packaging-alone cap can force medium.

### 5.5 Subtype / novel

- Subtype: if primary starts with `phage.` and `plasmid_score > 0.35` → `phage_plasmid_like`.
- Novel (`novel.py`): `nonstandard_IE_candidate`, `novel_tier` ∈ {high_interest, medium, watchlist, none} from `novel_outlier_criteria` (strong IE core while unclassified; ≥2 strong specialized modules; weak-module unclassified with adequate annotation coverage).

---

## 6. Markers (58 keys) — by category

Weights from current YAML (2026-08-08). Full Pfam/regex lists: `ie_atlas_labels.yaml` → `module_markers`.

| Category | Keys (weight) |
|----------|---------------|
| IE_core | integrase(0.54) |
| phage_packaging | terl(0.50), ters(0.35) |
| phage_head | hk97_mcp(0.54), portal(0.54), mcp(0.30) |
| phage_tail | tail(0.54) |
| phage_connector | connector(0.44) |
| phage_lysis | lysis(0.34) — **not** in `phage_score` |
| phage_reg | phage_reg(0.24) — feeds `pici_score` |
| inovirus | inov_zot(0.55), inov_rcr(0.40), inov_rep(0.35), inov_ssb(0.30), inov_coat(0.30), inov_assembly(0.25) |
| conjugation | relaxase(0.54), t4cp(0.54), virb4(0.54), mpf(0.30) |
| AICE / DSS | ftsk_traB(0.52); tdtA_like(0.52), nurA(0.416), rm(0.312) |
| replication | primpol_aep(0.54), dna_pol_a/b(0.54), helicase(0.43), repABC(0.48), ssb(0.28) |
| plasmid / host | par_partition, ta_system; ribosomal, metabolic |
| PICI / satellite | pici_cpt(0.48); gta_marker, satellite_p4, genomad_provirus, virsorter_viral |
| mobilome / IS | tn3_resolvase, tn_transposase, inti/attc integron, group2_rt, casposon; is_transposase |
| island | pai_virulence, rm_island, ta_rich |
| competence / non-T4SS | com_ea, com_machinery; conj_tad, conj_t4ap |
| SCC | ccr, cch, cch2, ccpol, mp_partner — **unused in scoring/rules** |

---

## 7. Locked cohort 6960 — which atlas is “live”

| Artifact | Path | Role |
|----------|------|------|
| Cohort definition | `artifacts/attL_complement_blast_window/cohort_6960_DEFINITION.md` | Locked set |
| Atlas join | `.../cohort_6960_atlas_join.tsv` | **Canonical labels for 6960** |
| Source | `stage5_atlas_join.tsv` ← mostly `ie_instance_classification_full25591_strict.tsv` | Older taxonomy |
| Fig-4 collapse | `plot_cohort_6960_four_class_lengths.py` | 4 paper classes |
| Figure | `figures/attL_complement_blast_window/11_four_func_class_attL_ie_lengths.png` | |

**6960 coarse (stage5 join):** unclassified 3738 (53.7%), conj 963, phage_gta 956, rep 951, phage_caudoviral 351, phage_inovirus 1.  
**`mob_element` / `island_module` / `phage_pici` = 0** on this join — those labels exist in current YAML/NB7/stage9, but the strict/stage5 table predates them.

Four-class paper merge:

| Paper class | Coarse sources |
|-------------|----------------|
| Tailed phage | phage_caudoviral + phage_inovirus |
| Conjugative | conj |
| Phage-like | phage_gta + rep |
| Unknown | unclassified |

Stage9 / current NB7 (`stage9_atlas_join.tsv`, n≈7507) still has mob/island/pici populated; ~1735/6663 overlap disagree with 6960 join on coarse — **classifier-version mismatch**, not just filter effect.

---

## 8. File map

### Config / docs
- `.../notebook7_extended/ie_atlas_labels.yaml` — **source of truth for thresholds**
- `.../IE_ATLAS_CLASSIFIER.md` — overview (priority order in text is partly stale vs YAML)
- `.../ie_atlas_markers_reference.md` — marker table (may lag weights/`sum_scale`)
- `final_final/supplementary/diploma_text_IE_atlas.txt`
- `final_final/supplementary/S_Table_IE_atlas_markers.xlsx`
- **This corpus:** `final_final/docs/corpus/IE_ATLAS.md`

### Package `ie_classifier/`
| File | Responsibility |
|------|----------------|
| `classify.py` | orchestration |
| `marker_match.py` | Pfam/regex hits |
| `scores.py` | module scores, flags, synteny |
| `rules.py` | priority gates, coarse, confidence, packaging cap |
| `novel.py` | out-of-bin tiers |
| `schema.py` | public columns / evidence |
| `loaders.py` | GBFF + eggNOG hydrate |
| `txsscan.py` / `external_detectors.py` | optional evidence |
| `atlas_learn.py` | stub only |

---

## 9. Open inconsistencies (do not paper over)

1. **6960 vs current YAML/NB7** — missing mob/island/PICI on locked join; re-run current classifier on 6960 loci to restore.
2. **SCC incomplete** — markers exist, no score/rule path.
3. **Docs drift** — old canvas/reference still mention `phage.caudoviral.defective`, SCC labels, older weights (0.52 vs 0.54) and `sum_scale: 1.2`.
4. **Diploma %** match stage9/NB7 (~7507), not the 6960 join.
5. Optional detectors often off in stage5 rebuild path.

---

## 10. Quick re-run recipes

```bash
# Current classifier on a truth table
python results_organized/scripts/run_ie_atlas_classification.py \
  --truth PATH/truth.tsv \
  --locus-merge results_organized/genome_clustering/project/notebook7_extended/locus_eggnog_merged.tsv \
  --rep-index results_organized/merged_deduplicated/rep_index.tsv \
  --gbff-dir results_organized/merged_deduplicated \
  --config results_organized/genome_clustering/project/notebook7_extended/ie_atlas_labels.yaml \
  --out-tsv PATH/ie_instance_classification.tsv

# Stage5 rebuild
python final_final/scripts/build_stage5_atlas.py

# BLAST-window atlas
python final_final/scripts/build_blast_window_40kb_atlas.py
```

---

## 11. Next corpus slices (project-wide, later)

Suggested siblings under `final_final/docs/corpus/`:

- `FINDER_PIPELINE.md` — Snakemake IE search / attL
- `COHORT_6960.md` — lock definition, filters, numbers
- `ISPOLATOV_MODEL.md` — random-match null
- `TAXONOMY_GTDB.md` — Fig 3 host taxonomy
- `TANDEM_IE.md` — neighbour/array geometry
