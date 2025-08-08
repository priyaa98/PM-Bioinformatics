# Whole Exome Sequencing (WES) Pipeline

This folder contains a **complete** whole exome sequencing (WES) pipeline for human variant discovery, annotation, and filtering.  
It is intended to be run **with your own sequencing data** and the appropriate reference and annotation resources.

## Overview

The pipeline follows a standard short-read variant calling and annotation workflow:

1. **Trim Galore** – adapter/quality trimming  
2. **BWA-MEM** – alignment to reference  
3. **Picard MarkDuplicates** – duplicate marking and indexing  
4. **GATK BaseRecalibrator + ApplyBQSR** – base quality score recalibration  
5. **GATK HaplotypeCaller** – per-sample variant calling in GVCF mode  
6. **Variant Annotation with VEP** – functional consequence prediction, population frequency, and pathogenicity scoring  
7. **Variant Filtering & Segregation Analysis with VASE** – inheritance modelling, ClinVar integration, and final case-level reports  

> ⚠️ This pipeline **cannot run without real data**. You must provide your own FASTQs, reference genome, and known-sites resources.  
> The annotation and filtering steps require large datasets such as gnomAD, CADD, and ClinVar.

---

## Requirements

These must be installed and available on your PATH:

- `bash`
- `python` (≥3.8)
- `bwa`
- `samtools`
- `picard`
- `trim_galore`
- `gatk`
- `bcftools`
- `vep` (Ensembl VEP with GRCh38 cache and plugins: gnomAD, CADD, REVEL, AlphaMissense)
- `vase` (with ClinVar data and prepared pedigree files)

---

## Reference Resource Versions

The pipeline was developed and tested between **February 2024** and **September 2024** using the following resource versions:

| Resource        | Version |
|-----------------|---------|
| **Reference genome** | GRCh38 (Homo_sapiens_assembly38) |
| **dbSNP**       | v155 |
| **gnomAD**      | v4.0 |
| **CADD**        | v1.7 |
| **REVEL**       | 2020 release (latest for GRCh38 at the time) |
| **AlphaMissense** | 2023 GRCh38 release |
| **ClinVar**     | July 2024 release |
| **GATK**        | v4.6.0.0 |
| **Picard**      | v2.27.5 |
| **Trim Galore** | v0.6.10 |
| **BWA**         | v0.7.17 |
| **samtools**    | v1.18 |
| **bcftools**    | v1.18 |
| **VEP**         | v110 (with GRCh38 cache from 2024-07) |
| **VASE**        | 2024.05.x (mid-2024 release) |

> ⚠️ If you use newer versions of these tools or resources, some command-line options or output formats may change.  
> Adjust the workflow scripts accordingly.
