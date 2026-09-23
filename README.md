# Olabisi

Reproducible bulk RNA-seq analysis of iPSC-derived podocytes for the OCEAN
manuscript. The workflow performs preprocessing, a paired IFNG-response
interaction analysis, empirical-null p-value calibration, CollecTRI/PROGENy
ULM enrichment, and generation of manuscript Figures 5B, 5C, 5E, 5F, and 5G.

## Objective

- Main research question: do podocytes from high-risk cases respond differently
  to IFNG stimulation than podocytes from high-risk controls?

- Analyses include QC, differential expression, TF/pathway activity inference,
  and manuscript figure/source-data generation.

## Data

- Human donors from two groups: APOL1 high-risk genotype FSGS donor (HR cases), APOL1 high-risk genotype healthy controls (HR controls).

- PBMCs were reprogrammed to iPSCs and differentiated into podocyte-like cells.

- Bulk RNA-seq dataset.

- Each donor has IFNG- and IFNG+ samples.

## Software environments

Snakemake 9 is installed in a small orchestration environment. Create and
activate it once with:

```bash
mamba env create --file workflow/envs/snakemake9.yaml
conda activate snakemake9
mamba config set channel_priority strict
snakemake --version
```

The workflow maintains exactly two rule-level environments:

- `workflow/envs/pyenv.yaml` for Python DE analysis and plotting.
- `workflow/envs/renv.yaml` for R p-value calibration.

Snakemake creates those environments automatically.

## Run the paper workflow

From the repository root, preview and then execute the complete paper target:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 paper --software-deployment-method conda -n
snakemake --snakefile workflow/Snakefile --cores 1 paper --software-deployment-method conda
```

To retrieve a missing prior-knowledge snapshot explicitly:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 download_prior_knowledge --software-deployment-method conda
```

The timestamped Ensembl, CollecTRI, and PROGENy files distributed with the
archived data are preferred. Retrieval checks the configured SHA-256 and fails
if an upstream resource has changed, preventing an unnoticed change in the
analysis.