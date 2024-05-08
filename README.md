# MeQTL Analysis of EPIC data with blood samples from UK cohorts

This repository contains code and documentation related to the analysis presented in [Villicaña et al. (2023), *Gen Bio* 24, 176](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03011-x). The analysis focuses on the investigation of DNA methylation quantitative trait loci (meQTLs) using data from 2,358 blood samples across three UK cohorts, profiled with the Illumina Infinium MethylationEPIC array.

## Analysis Overview

The analysis pipeline is structured into several modules. Below is a brief overview of each module:

### MeQTL Analysis by Cohort

- `1_regression`: Conducts DNAm rank-based inverse normal transformation and adjusts for covariates.
- `2_meQTL`: Performs meQTL analysis using the `MatrixEQTL` package. The analysis is executed in chunks by chromosomes. Files from the same cohort are merged and converted to `GWAMA` format using scripts `merge_files.sh` and `qtl2gwama.R` located in the `general` directory.

### MeQTL Meta-Analysis

- `3_meta_analysis`: Conducts `GWAMA` random-effects meta-analysis separately for *cis* and *trans* meQTLs. Then, annotates and filters results by *P*-value estimated with permutations from `5_permutations_ma`. Plug-in FDR computation is performed using `general/graphs.R`.

### Permutations

- `4_permutations`: Performs twenty meQTL analyses with shuffled samples.
- `5_permutations_ma`: Conducts meta-analysis of permutation meQTL results.

### Follow-up Analysis

- `6_annotation`: Annotation of SNPs and CpGs, enrichment analyses, and LD-clumping.
- `7_SMR`: Implements Summary-based Mendelian Randomisation for GWAS traits and *cis*-eQTLs.
- `8_validation`: Validates CpGs with meQTLs based on GoDMC catalogue.
- `9_replication`: Replicates eight CpG meQTLs using MeDIP-seq data.
- `10_450K_split_ma`: Calculates plug-in FDR for sets of 450K legacy and EPIC-only probes separately.
- `11_cell_interaction`: Investigates cell interacting meQTLs with monocytes or CD4T cells.

### Heritability

The `heritability` directory contains the code for heritability estimation using `OpenMX`, utilising twin pairs from TwinsUK.

## Citation

If you use the code or data from this repository, please cite:

Villicaña, S., Castillo-Fernandez, J., Hannon, E. et al. Genetic impacts on DNA methylation help elucidate regulatory genomic processes. *Genome Biol* **24**, 176 (2023), [DOI](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03011-x).

