# Differential Expression Analysis (RNA-seq)

This folder contains the **Differential Expression (DE) workflow** used for our RNA-seq project. 

---

## 1. Overview

- The goal was to identify **differentially expressed genes (DEGs)** between experimental conditions  
  (e.g. resistant vs sensitive groups).  
- Analysis is performed at the **transcript → gene → functional annotation** level.  
- Tools used follow standard RNA-seq DE pipelines:  
  - **HPC preprocessing**: FastQC, Trim Galore, Salmon  
  - **R analysis**: tximport, DESeq2, visualization, functional annotation  

---

## 2. Methods

### 2.1 HPC preprocessing (arc4)
- Raw FASTQs were processed on the **arc4 HPC cluster** using batch scripts (`qsub` jobs).  
- Steps:
  1. **FastQC** — raw quality check  
  2. **Trim Galore** — adapter/quality trimming (cutadapt backend)  
  3. **Salmon** — transcript-level quantification with selective alignment  
- These outputs (`quant.sf`) were then imported into R for downstream analysis.  
> ⚠️ The actual bash scripts are not included here, but the workflow is fully described.

### 2.2 R-based DE analysis
- Conducted in **R (≥4.1)** with the following packages:
  - `DESeq2`, `tximport`, `ggplot2`, `pheatmap`, `EnhancedVolcano`,  
    `dplyr`, `tidyr`, `clusterProfiler`, `org.Hs.eg.db`
- Workflow:
  1. Import transcript-level data with **tximport**  
  2. Build a **DESeq2 dataset** with sample metadata  
  3. Perform normalization and variance stabilizing transformation  
  4. Run DESeq2 for differential expression testing  
  5. Generate diagnostic and result plots:
     - PCA plot
     - Sample distance heatmap
     - MA plots (before/after shrinkage)
     - Volcano plot
     - Top-gene heatmaps  
  6. Functional annotation with **clusterProfiler** (GO enrichment for BP/MF/CC)
