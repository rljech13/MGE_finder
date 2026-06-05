# End-to-end stress test panel

Six genomes selected to exercise known edge cases in the IE_finder tool (`finder_pipeline/`):

| Sample | Source | Stress case |
|--------|--------|-------------|
| `e2e_positive_thermus` | GCA_900604845.1_TTHNAR1 | Known Thermus reference with integrative elements |
| `e2e_reference_thermus` | GCA_000008125.1_ASM812v1 | Second Thermus positive control (representative in atlas) |
| `e2e_ambiguous_n` | GCA_000376665.1_ASM37666v1 | Assembly contains ambiguous N (IE N-filter) |
| `e2e_hyperfragmented` | GCA_043732285.1_ASM4373228v1 | 1772 contigs (fragmentation / coordinate handling) |
| `e2e_long_contig` | GCA_965249745.1_glEphLana1 | 5.1 Mbp single contig (Prodigal / long ORFs) |
| `e2e_small_genome` | GCA_019091385.1_ASM1909138v1 | Smallest local Deinococcales draft (~1.2 Mbp) |

Run:

```bash
cd finder_pipeline
CONFIG=e2e_test/ie_finder_config_e2e.yaml ./run.sh --cores=6
```
