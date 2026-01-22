# ==============================================================================
# PROJECT: MOFA+ Multi-Omics Integration (TCGA Breast Cancer)
# AUTHOR: [Your Name/GitHub Handle]
# DATE: [Current Date]
# ==============================================================================
# DESCRIPTION:
# This script performs multi-omics integration of TCGA Breast Cancer data using MOFA+.
#
# PIPELINE OVERVIEW:
# 1. Data Cleaning: Formats barcodes, removes version numbers, filters Y-chr.
# 2. Normalization: DESeq2 (RNA), Z-score (CNV/RPPA).
# 3. Model Training: MOFA+ with 20 factors.
# 4. Downstream Analysis: Variance decomposition, Feature weights, Clustering.
# 5. Biological Interpretation: GSEA (Reactome) with custom binary matrix mapping.
# ==============================================================================

# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
# setwd("/path/to/your/project") # <--- UNCOMMENT AND SET WORKING DIRECTORY

set.seed(123) # Ensure reproducibility

# List of required packages
packages <- c("MOFA2", "DESeq2", "matrixStats", "maftools", "dplyr", 
              "survival", "survminer", "ggplot2", "reshape2", "biomaRt", 
              "msigdbr", "org.Hs.eg.db")

# Install missing packages automatically
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

# Load libraries
suppressPackageStartupMessages({
  lapply(packages, library, character.only = TRUE)
})

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# A. Clean TCGA Barcodes (Standardizes to TCGA-XX-XXXX format)
clean_barcodes <- function(x) {
  cols <- if(is.data.frame(x) || is.matrix(x)) colnames(x) else x
  substr(gsub("\\.", "-", cols), 1, 12)
}

# B. Clean Gene IDs (Removes version numbers: ENSG...1 -> ENSG...)
clean_gene_ids <- function(mat) {
  clean_ids <- sub("\\..*", "", rownames(mat))
  # Aggregates counts for duplicate IDs (e.g. PAR regions)
  mat_clean <- rowsum(mat, group = clean_ids)
  return(mat_clean)
}

# C. Robust Normalization (DESeq2 with fallback to robust variance stabilization)
normalize_deseq2 <- function(counts_mat, name) {
  message(paste("Normalizing", name, "with DESeq2..."))
  counts_mat <- round(as.matrix(counts_mat)) # Integers required for DESeq2
  
  colData <- data.frame(row.names = colnames(counts_mat))
  dds <- DESeqDataSetFromMatrix(counts_mat, colData, design = ~1)
  
  # Filter low counts (keep genes with >10 counts total)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Try fast VST, fallback to robust if it fails
  vst_out <- tryCatch({
    vst(dds, blind = TRUE)
  }, error = function(e) {
    message(paste("  (Fast VST failed for", name, "- switching to robust variance stabilization)"))
    varianceStabilizingTransformation(dds, blind = TRUE)
  })
  return(assay(vst_out))
}

# D. CNV Normalization (Winsorization + Z-score)
normalize_cnv <- function(mat) {
  message("Normalizing CNV...")
  mat <- as.matrix(mat)
  mat[mat > 3] <- 3   # Cap extreme gains
  mat[mat < -3] <- -3 # Cap extreme deletions
  return(t(scale(t(mat)))) # Z-score standardization per gene
}

# E. RPPA Normalization (Log2 + Z-score)
normalize_rppa <- function(mat) {
  message("Normalizing RPPA...")
  mat <- as.matrix(mat)
  if(max(mat, na.rm=T) > 50) mat <- log2(mat + 1) # Log transform if raw values
  return(t(scale(t(mat))))
}

# F. Feature Selection (Select top N variable features)
select_top <- function(mat, n) {
  vars <- rowVars(mat)
  mat[order(vars, decreasing=TRUE)[1:min(n, nrow(mat))], ]
}

# ==============================================================================
# 3. DATA LOADING & PREPROCESSING
# ==============================================================================
message("Loading raw data...")
# Ensure these files exist in your directory
rna_raw  <- read.csv("all_mrna.csv", row.names = 1)
cnv_raw  <- read.csv("all_cnv.csv", row.names = 1)
rppa_raw <- read.csv("protein.csv", row.names = 1)

# Standardize Sample IDs
colnames(rna_raw)  <- clean_barcodes(rna_raw)
colnames(cnv_raw)  <- clean_barcodes(cnv_raw)
colnames(rppa_raw) <- clean_barcodes(rppa_raw)

# Clean Feature IDs
rna_raw <- clean_gene_ids(rna_raw)

# Fix Row Names for CNV and RPPA (Handle specific input formats)
# (Assumes specific column structure from user input)
if("Gene-Symbol" %in% colnames(cnv_raw)) {
  rownames(cnv_raw) <- cnv_raw$`Gene-Symbol`
  cnv_raw <- cnv_raw[ , -which(names(cnv_raw) %in% c("Gene-Symbol"))]
}
if("sample" %in% colnames(rppa_raw)) {
  rownames(rppa_raw) <- rppa_raw$sample
  rppa_raw <- rppa_raw[ , -which(names(rppa_raw) %in% c("sample"))]
}

# ------------------------------------------------------------------------------
# Y-CHROMOSOME FILTERING (Quality Control)
# ------------------------------------------------------------------------------
message("Filtering Y-chromosome features to remove gender bias...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
y_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
                 filters = "chromosome_name", values = "Y", mart = mart)

# Filter RNA (Ensembl IDs)
rna_raw <- rna_raw[!rownames(rna_raw) %in% y_genes$ensembl_gene_id, ]

# Filter CNV (Symbol or Ensembl)
if(startsWith(rownames(cnv_raw)[1], "ENSG")) {
  cnv_raw <- cnv_raw[!rownames(cnv_raw) %in% y_genes$ensembl_gene_id, ]
} else {
  cnv_raw <- cnv_raw[!rownames(cnv_raw) %in% y_genes$hgnc_symbol, ]
}

# ------------------------------------------------------------------------------
# INTERSECTION & NORMALIZATION
# ------------------------------------------------------------------------------
# Find common samples across all views
common_samples <- Reduce(intersect, list(colnames(rna_raw), colnames(cnv_raw), colnames(rppa_raw)))
message(paste("Common samples retained:", length(common_samples)))

# Subset data
rna_mat  <- rna_raw[, common_samples]
cnv_mat  <- cnv_raw[, common_samples]
rppa_mat <- rppa_raw[, common_samples]

# Normalize
rna_norm  <- normalize_deseq2(rna_mat, "RNA")
cnv_norm  <- normalize_cnv(cnv_mat)
rppa_norm <- normalize_rppa(rppa_mat)

# Select Highly Variable Features (HVGs)
rna_final  <- select_top(rna_norm, 5000)
cnv_final  <- select_top(cnv_norm, 5000)
rppa_final <- select_top(rppa_norm, 200)

# ==============================================================================
# 4. MOFA+ MODEL TRAINING
# ==============================================================================
# Create MOFA Object
MOFA_obj <- create_mofa(list("RNA" = rna_final, 
                             "CNV" = cnv_final, 
                             "RPPA" = rppa_final))

# Define Options
data_opts <- get_default_data_options(MOFA_obj)
model_opts <- get_default_model_options(MOFA_obj)
model_opts$num_factors <- 20
model_opts$likelihoods <- c("gaussian", "gaussian", "gaussian")
train_opts <- get_default_training_options(MOFA_obj)
train_opts$convergence_mode <- "medium" 

# Prepare and Run
MOFA_model <- prepare_mofa(MOFA_obj, data_opts, model_opts, train_opts)
MOFA_model <- run_mofa(MOFA_model, outfile = "project_MOFA_Final.hdf5")

# Save Reference
saveRDS(MOFA_model, "project_MOFA_Final.rds")

# ==============================================================================
# 5. METADATA INTEGRATION
# ==============================================================================
# Load clinical data
clinical <- read.csv("annotations.csv") 
clinical$patient <- clean_barcodes(clinical$patient)

# Extract sample names from the model (to ensure order)
model_samples <- samples_names(MOFA_model)$group1

# Align metadata to model samples
match_indices <- match(model_samples, clinical$patient)
clinic_filtered <- clinical[match_indices, ]

# Assign correct rownames and sample column
rownames(clinic_filtered) <- model_samples
clinic_filtered$sample <- model_samples # Mandatory column for MOFA

# Attach to model
samples_metadata(MOFA_model) <- clinic_filtered







message("Metadata successfully attached.")

# ==============================================================================
# 6. DOWNSTREAM ANALYSIS & VISUALIZATION
# ==============================================================================

# A. Variance Decomposition
# ------------------------------------------------------------------------------
# Check how much variance is explained by each factor per view
r2_plot <- plot_variance_explained(MOFA_model, x="view", y="factor", plot_total = TRUE)
print(r2_plot[[2]]) # Print total variance plot

# B. Dimension Reduction (t-SNE)
# ------------------------------------------------------------------------------
# Embed samples based on factors
MOFA_model <- run_tsne(MOFA_model)
# We filter out real NAs and literal "NA" strings just in case
samples_to_keep <- clinic_filtered$sample[ !is.na(clinic_filtered$BRCA_Subtype_PAM50) & clinic_filtered$BRCA_Subtype_PAM50 != "NA" ]

# 3. Create a clean model by removing the NA samples
MOFA_model_clean <- subset_samples(MOFA_model, samples_to_keep)

plot_dimred(MOFA_model_clean,
            method = "TSNE",
            color_by = "BRCA_Subtype_PAM50")

# C. Factor Visualization
# ------------------------------------------------------------------------------
# Explore specific factors of interest (e.g., Factors 4-6)
plot_factor(MOFA_model_clean, 
            factor = 1:6,
            color_by = "BRCA_Subtype_PAM50",
            shape_by = "vital_status")

# D. Feature Weights (Identify Drivers)
# ------------------------------------------------------------------------------
# Plot top weights for RNA Factor 1
plot_weights(MOFA_model,
             view = "RNA",
             factor = 1,
             nfeatures = 10,
             scale = TRUE)
# Plot top weights for CNV Factor 1
plot_weights(MOFA_model,
             view = "CNV",
             factor = 4,
             nfeatures = 10,
             scale = TRUE)

# Plot top weights for RPPA Factor 1
plot_weights(MOFA_model,
             view = "RPPA",
             factor = 1,
             nfeatures = 10,
             scale = TRUE)
# ==============================================================================
# 7. GENE SET ENRICHMENT ANALYSIS (GSEA)
# ==============================================================================

# A. Convert Ensembl IDs to Gene Symbols
# ------------------------------------------------------------------------------
# We convert the model features to Symbols for readable plots
message("Converting Ensembl IDs to Gene Symbols...")
current_ids <- features_names(MOFA_model)[["RNA"]]

symbols <- mapIds(org.Hs.eg.db,
                  keys = current_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Handle NAs and Duplicates
symbols[is.na(symbols)] <- current_ids[is.na(symbols)]
symbols <- make.unique(symbols)

# Update model
features_names(MOFA_model)[["RNA"]] <- symbols

# B. Run Enrichment (Robust Binary Matrix Method)
# ------------------------------------------------------------------------------
# Note: Since we just changed the model to Symbols, we must use Symbol-based gene sets.
# If using 'reactomeGS' (which is Ensembl), we cannot use the Symbol-updated model directly.
# Below uses 'msigdbr' to fetch SYMBOL-based pathways to match our new model names.

message("Running GSEA...")
# perform Enrichment Analysis on mRNA data using pre-build Reactome gene sets
data("reactomeGS", package = "MOFAdata")
fsea_rna.results <- run_enrichment(MOFA_model, view="RNA", feature.sets=reactomeGS)


# Visualize Results
plot_enrichment_heatmap(fsea_rna.results, n_features = 10, log_p_value = TRUE)



plot_enrichment(fsea_rna.results, 
                factor = 2, 
                max.pathways = 15
)



plot_enrichment_detailed(fsea_rna.results, 
                         factor = 1, 
                         max.genes = 8, 
                         max.pathways = 5
)



message("Analysis Complete. Script finished successfully.")


