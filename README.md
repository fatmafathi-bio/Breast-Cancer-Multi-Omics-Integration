# TCGA Breast Cancer Multi-Omics Integration using MOFA+

## Project Overview
This project performs a multi-omics integration of TCGA Breast Cancer (BRCA) data using **Multi-Omics Factor Analysis (MOFA+)**. The goal is to identify latent factors that drive biological variation across transcriptomics (RNA-seq), copy number variations (CNV), and proteomics (RPPA), and link these factors to patient survival.

## Data Sources
The analysis integrates three distinct data modalities from TCGA:
* **mRNA Expression (RNA):** 5,000 highly variable genes (normalized via DESeq2).
* **Copy Number Variation (CNV):** 5,000 variable sites (Gistic scores).
* **Proteomics (RPPA):** 200 proteins (Z-score normalized).

## Workflow
1.  **Preprocessing:** * Data cleaning and barcode standardization.
    * Removal of Y-chromosome genes to prevent gender bias.
    * Normalization (DESeq2 for RNA, Z-score for CNV/RPPA).
2.  **Model Training:** * MOFA+ model trained with 20 factors.
3.  **Downstream Analysis:**
    * Variance decomposition per view.
    * Gene Set Enrichment Analysis (GSEA) using Reactome pathways.
    * Survival Analysis (Cox Proportional Hazards & K-Means clustering).

## Key Findings
Based on the analysis results:
* **Factor 1** captures the largest variance in RNA and identifies Basal/TNBC tumors (Keratinization) with HDAC pathway enrichment, suggesting vulnerability to HDAC inhibitors (e.g., Vorinostat).
* **Pathway Enrichment:** Factor 2,3 is significantly enriched for **Immune System**, **Interferon Signaling**,**TCR signaling**,**Costimulation by the CD28 family** and **PD-1 Signaling** pathways witch can predict immuotherapy efficacy 

**Factor 4** Unique variance in **CNV** layer and it specifically highlights **MET activates PTK2 signaling**,**Signaling by MET**, and **MET promotes cell motility** which can predict that patients enriched in **Factor 4** are potential candidates for **MET inhibitors**. 

**Factor 5** is likely invasive and drug-resistant due to **fibrotic tumor microenvironment**, which typically promotes tumor aggression and protects cancer cells from chemotherapy (drug delivery barrier).

