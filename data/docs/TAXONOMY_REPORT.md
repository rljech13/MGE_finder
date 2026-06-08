# Deinococcales genome taxonomy report

## a) Download process status

Historical snapshot from a bulk download run:

- Tmux session: `deinococcales_download`
- Target set: **362 genomes** (at time of report)
- Live progress: `logs/download_progress.log`

### Check status:

```bash
cd IE_finder/data
find ncbi/Deinococcales -name "*.fna.gz" | wc -l
tail -f logs/download_progress.log
```

## b) Taxonomic composition

### Order: Deinococcales

Based on analysis of the first 100 genomes from a 362-genome run:

### Main genera and species:

1. **Deinococcus** (main genus)
   - Deinococcus radiodurans (3 genomes)
   - Deinococcus deserti
   - Deinococcus aquaticus (2 genomes)
   - Deinococcus sonorensis (2 genomes)
   - Deinococcus metalli (2 genomes)
   - Deinococcus saxicola (2 genomes)
   - Deinococcus wulumuqiensis (2 genomes)
   - Deinococcus altitudinis
   - Deinococcus antarcticus
   - Deinococcus budaensis
   - Deinococcus caeni
   - Deinococcus cellulosilyticus
   - Deinococcus ficus
   - Deinococcus frigens
   - Deinococcus grandis
   - Deinococcus hohokamensis
   - Deinococcus lacus
   - Deinococcus malanensis
   - Deinococcus marmoris
   - Deinococcus multiflagellatus
   - Deinococcus navajonensis
   - Deinococcus oregonensis
   - Deinococcus petrolearius
   - Deinococcus radiomollis
   - Deinococcus radiophilus
   - Deinococcus radiopugnans
   - Deinococcus rufus
   - Deinococcus soli
   - Deinococcus taklimakanensis
   - Deinococcus yunweiensis
   - Deinococcus sp. (multiple undescribed strains)

2. **Uncultured bacteria**
   - Deinococcaceae bacterium (6 genomes)
   - Deinococcales bacterium (13 genomes)
   - uncultured Deinococcaceae bacterium (7 genomes)
   - uncultured Deinococcales bacterium (4 genomes)
   - uncultured Deinococcus sp. (4 genomes)

3. **Marinithermaceae bacterium** (8 genomes)
   - Representatives of family Marinithermaceae

### Families in order Deinococcales:

1. **Deinococcaceae** — main family
   - Genus: Deinococcus
   - Uncultured representatives

2. **Marinithermaceae** — present in sample
   - Uncultured bacteria

3. **Trueperaceae** — may appear in full set
   - Genus: Truepera (not seen in first 100 genomes)

### Statistics (from 100-genome sample):

- **Unique species**: ~48
- **Main genus**: Deinococcus
- **Uncultured**: ~36 genomes (36%)
- **Described species**: ~64 genomes (64%)

### Notes:

1. Order Deinococcales includes several families:
   - Deinococcaceae (Deinococcus)
   - Trueperaceae (Truepera)
   - Marinithermaceae
   - and others

2. The 100-genome sample is dominated by **Deinococcus**

3. Full analysis of all 362 genomes may reveal additional genera (e.g. Truepera)

4. Many genomes are uncultured or undescribed species

## Commands for verification:

```bash
cd IE_finder/data
tail -f logs/download_progress.log
find ncbi/Deinococcales -name "*.fna.gz" | wc -l
python3 scripts/check_taxonomy.py
```
