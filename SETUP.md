# Setup Guide for MGE Finder

This guide will help you set up MGE Finder from scratch.

## Prerequisites

- Conda or Mamba
- Python 3.10+
- At least 50GB free disk space (for results and databases)

## Step 1: Clone the Repository

```bash
git clone <repository-url>
cd MGE_finder
```

## Step 2: Create Conda Environment

```bash
mamba env create -f envs/MGE_finder.yaml
conda activate MGE_finder
```

## Step 3: Install Additional Dependencies

### For Bakta Pipeline

```bash
conda create -n bakta -c bioconda -c conda-forge bakta
conda activate bakta
# Download Bakta database
bakta_db download --output /path/to/bakta/db
```

### For PADLOC Pipeline

```bash
conda create -n padloc -c bioconda padloc
conda activate padloc
padloc download_db
```

### For vContact Classification

```bash
conda create -n vContact3 -c bioconda vcontact2
# or for vContact3:
conda create -n vContact3 -c bioconda vcontact3
```

## Step 4: Download Pfam HMM Profiles

Download Pfam HMM profiles and place them in the `pfam/` directory:

```bash
# Example: Download PF00589 (Integrase core domain)
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00589/hmm -O pfam/PF00589.hmm

# Example: Download PF22022 (Integrase zinc-binding domain)
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF22022/hmm -O pfam/PF22022.hmm
```

## Step 5: Configure Pipelines

### Finder Pipeline

```bash
cd finder_pipeline
cp finder_config.yaml.example finder_config.yaml
# Edit finder_config.yaml with your paths
nano finder_config.yaml
```

Update the following in `finder_config.yaml`:
- `input_sources`: Paths to your genome directories
- `pfam_profiles`: Paths to your Pfam HMM files

### Bakta Pipeline

```bash
cd bakta_pipeline
cp bakta_config.yaml.example bakta_config.yaml
# Edit bakta_config.yaml
nano bakta_config.yaml
```

Update:
- `db`: Path to your Bakta database

### Other Pipelines

Similarly, copy and configure:
- `padloc_pipeline/padloc_config.yaml.example`
- `classification_investigation/clusterMGE/mge_cluster_config.yaml.example`
- `classification_investigation/vcontact/vcontact_config.yaml.example`
- `protein_clusterization/classification_config.yaml.example`

## Step 6: Prepare Input Genomes

Place your genome FASTA files in the directories specified in `finder_config.yaml` under `input_sources`.

Supported formats:
- NCBI format: `*_genomic.fna`
- Hybracter format: `barcode*.fastq_final.fasta`

## Step 7: Run the Pipeline

See the main [README.md](README.md) for detailed instructions on running each module.

## Troubleshooting

### Conda Environment Issues

If Snakemake complains about conda environments, ensure you're using the correct conda prefix:

```bash
snakemake --conda-prefix /path/to/your/conda/envs
```

### Path Issues

All paths in configuration files should be either:
- Absolute paths (starting with `/`)
- Relative paths from the pipeline directory

### Missing Dependencies

If you encounter missing dependencies, check that:
1. The conda environment is activated
2. All required tools are installed in the environment
3. The environment name in the config matches the actual environment name

## Next Steps

- Read the main [README.md](README.md) for detailed module documentation
- Check individual pipeline README files for module-specific information
- Review example configuration files for guidance

