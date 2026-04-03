# EasyDE Vignette: Differential Expression for One Contrast x Stratum

This vignette walks through a single EasyDE analysis from start to finish:
**Diabetic_vs_ND** contrast in **Beta** cells. Each step is explained with
the exact R code the pipeline runs internally. Copy-paste into an R console
or Jupyter notebook with an R kernel.

---

## Overview

EasyDE runs pseudobulk differential expression using DESeq2, with optional
RUVseq batch correction and fGSEA pathway enrichment. For one contrast in
one cell type, the steps are:

1. **Validate** inputs (config, metadata, count matrix)
2. **Prepare coldata** — filter samples, deduplicate donors, preflight checks
3. **Base DESeq2** — fit negative binomial GLM, extract DE results
4. **RUVseq** — estimate latent unwanted variation, select safe factors
5. **Final DESeq2** — re-fit with W factors from RUVseq
6. **fGSEA** — gene set enrichment on DE results

---

## Setup

```r
library(yaml)
library(data.table)
library(dplyr)
library(DESeq2)
library(RUVSeq)
library(EDASeq)
library(fgsea)
library(ggplot2)
library(ggrepel)

# Working directory: the EasyDE root
# setwd("/path/to/EasyDE")
```

---

## Configuration

```r
# Load pipeline config
cfg <- yaml.load_file("config/config_traits.yaml")

# Load contrasts
contrasts <- fread("contrasts/contrasts.csv")

# Our example: Diabetic_vs_ND in Beta cells
CONTRAST_ID  <- "Diabetic_vs_ND"
STRATUM      <- "Beta"
ctr          <- contrasts[contrast_id == CONTRAST_ID]

contrast_var <- ctr$contrast_var       # "Derived diabetes status"
control_grp  <- ctr$control_grp        # "Normal"
exp_grp      <- ctr$experimental_grp   # "Diabetes"
biosample_col <- ctr$biosample_id_col  # "sample_accession"

# Parse pipe-separated fields
covariates   <- strsplit(ctr$covariates, "\\|")[[1]]
strata       <- strsplit(ctr$strata, "\\|")[[1]]
cat("Covariates:", paste(covariates, collapse = ", "), "\n")
cat("Strata:", paste(strata, collapse = ", "), "\n")
```

---

## Step 1: Load and Inspect Inputs

```r
# Load sample metadata
meta <- fread(cfg$inputs$sample_metadata)
cat("Samples:", nrow(meta), "\n")
cat("Columns:", paste(head(colnames(meta), 15), collapse = ", "), "...\n")

# Load count matrix for Beta cells
counts_file <- file.path(cfg$inputs$counts_dir, paste0(STRATUM, ".counts.csv"))
counts_raw <- fread(counts_file, header = TRUE)
genes <- counts_raw[[1]]
counts_mat <- as.matrix(counts_raw[, -1])
rownames(counts_mat) <- genes
cat("Genes:", nrow(counts_mat), " Samples:", ncol(counts_mat), "\n")
```

---

## Step 2: Prepare Coldata (Sample Filtering)

This is the most complex step — it decides which samples survive into DE.

```r
# ── 2a. Filter metadata to this stratum ──────────────────────────
# Match stratum column in metadata (may be called "stratum", "cell_type", etc.)
stratum_col <- grep("stratum|cell.type|biosample_type", colnames(meta),
                     ignore.case = TRUE, value = TRUE)[1]
meta_stratum <- meta[meta[[stratum_col]] == STRATUM]
cat("Samples in", STRATUM, ":", nrow(meta_stratum), "\n")
```

```r
# ── 2b. Apply contrast-specific filters ──────────────────────────
# (e.g., filter to "no_treatment", "Control Without Diabetes", etc.)
for (i in 1:3) {
  fc <- ctr[[paste0("filter_col_", i)]]
  fv <- ctr[[paste0("filter_val_", i)]]
  if (!is.na(fc) && fc != "" && !is.na(fv) && fv != "") {
    allowed <- strsplit(fv, "\\|")[[1]]
    before <- nrow(meta_stratum)
    meta_stratum <- meta_stratum[meta_stratum[[fc]] %in% allowed]
    cat("  Filter", fc, "in [", paste(allowed, collapse=", "),
        "]: ", before, "->", nrow(meta_stratum), "\n")
  }
}
```

```r
# ── 2c. Subset to control + experimental groups ─────────────────
meta_stratum <- meta_stratum[
  meta_stratum[[contrast_var]] %in% c(control_grp, exp_grp)
]
cat("After group subsetting:", nrow(meta_stratum), "\n")
```

```r
# ── 2d. Deduplicate donors (keep first sample per donor) ────────
donor_col <- cfg$filtering$donor_col
meta_stratum <- meta_stratum[!duplicated(meta_stratum[[donor_col]]), ]
cat("After deduplication:", nrow(meta_stratum), "unique donors\n")
```

```r
# ── 2e. Align counts to surviving samples ────────────────────────
sample_ids <- meta_stratum[[biosample_col]]
common_samples <- intersect(sample_ids, colnames(counts_mat))
counts_filt <- counts_mat[, common_samples]
coldata <- as.data.frame(meta_stratum[match(common_samples, sample_ids), ])
rownames(coldata) <- common_samples
cat("Final: ", ncol(counts_filt), "samples x", nrow(counts_filt), "genes\n")
```

```r
# ── 2f. Preflight checks ────────────────────────────────────────
n_ctrl <- sum(coldata[[contrast_var]] == control_grp)
n_exp  <- sum(coldata[[contrast_var]] == exp_grp)
cat("Control (", control_grp, "):", n_ctrl, "samples\n")
cat("Experimental (", exp_grp, "):", n_exp, "samples\n")

if (n_ctrl < 3 || n_exp < 3) {
  stop("Too few samples in one group (need >= 3 per group)")
}
```

> **If preflight fails**: Common causes are aggressive filtering (too few
> donors pass all filters) or rare conditions. Options:
>
>   - Relax `filter_col` restrictions (e.g., include treated samples)
>   - Reduce `min_cells_per_sample` threshold
>   - Merge related conditions (e.g., "T1D" + "T2D" into "Diabetes")

```r
# ── 2g. Sanitize variable names for R formula ───────────────────
# DESeq2 formulas can't handle spaces or special characters
sanitize <- function(x) gsub("[^[:alnum:]]", "_", x)

orig_names <- colnames(coldata)
colnames(coldata) <- sanitize(colnames(coldata))
contrast_var_safe <- sanitize(contrast_var)
covariates_safe   <- sanitize(covariates)

# Drop covariates with < 2 unique values
for (cv in covariates_safe) {
  if (cv %in% colnames(coldata)) {
    n_unique <- length(unique(na.omit(coldata[[cv]])))
    if (n_unique < 2) {
      cat("  Dropping covariate", cv, "(only", n_unique, "unique value)\n")
      covariates_safe <- setdiff(covariates_safe, cv)
    }
  }
}

# Coerce character columns to factors
for (col in colnames(coldata)) {
  if (is.character(coldata[[col]])) {
    coldata[[col]] <- as.factor(coldata[[col]])
  }
}
```

> **If design matrix is rank-deficient**: This means some covariates are
> perfectly collinear. Options:
>
>   - Remove the redundant covariate from the `covariates` field in contrasts.csv
>   - Set `scale_numeric_covariates: true` in config to z-score numeric covariates
>   - Check if a factor has only one level after filtering (auto-dropped above)

---

## Step 3: Base DESeq2

```r
# ── 3a. Build design formula ────────────────────────────────────
# Covariates first, contrast variable last (DESeq2 convention)
active_covs <- intersect(covariates_safe, colnames(coldata))
formula_str <- paste("~", paste(c(active_covs, contrast_var_safe), collapse = " + "))
design_formula <- as.formula(formula_str)
cat("Design formula:", formula_str, "\n")
```

```r
# ── 3b. Create DESeqDataSet ─────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = design_formula
)

# Gene-level filtering: keep genes with >= 5 counts in >= 2 samples per group
keep <- rowSums(counts(dds) >= cfg$filtering$min_gene_counts) >= 2
dds <- dds[keep, ]
cat("Genes after filtering:", sum(keep), "\n")
```

```r
# ── 3c. Run DESeq2 ──────────────────────────────────────────────
dds <- DESeq(dds)
```

> **If DESeq2 fails with "every gene contains at least one zero"**: This
> means the count matrix is heavily zero-inflated (common in single-cell
> pseudobulk). Fix:
>
>   - Use `sfType = "poscounts"` in `DESeq()`:
>     `dds <- DESeq(dds, sfType = "poscounts")`
>   - This uses a positive-count estimator for size factors instead of the
>     median-of-ratios method

```r
# ── 3d. Extract results ─────────────────────────────────────────
res_base <- results(dds, contrast = c(contrast_var_safe, exp_grp, control_grp))
res_base <- as.data.frame(res_base)
res_base$gene <- rownames(res_base)

fdr_thr <- cfg$deseq2$fdr_threshold     # 0.05
lfc_thr <- cfg$deseq2$l2fc_threshold    # 0

n_degs <- sum(res_base$padj < fdr_thr &
              abs(res_base$log2FoldChange) > lfc_thr, na.rm = TRUE)
cat("DEGs (base model):", n_degs, "\n")
```

```r
# ── 3e. Optional: LFC shrinkage ─────────────────────────────────
# Shrinks noisy log-fold-changes toward zero (better for ranking/visualization)
if (cfg$deseq2$lfc_shrinkage) {
  library(apeglm)
  coef_name <- resultsNames(dds)[length(resultsNames(dds))]
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  res_base$log2FoldChange_shrunk <- res_shrunk$log2FoldChange
}
```

```r
# ── 3f. Volcano plot ────────────────────────────────────────────
ggplot(res_base, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > lfc_thr),
             size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("grey70", "red3")) +
  geom_hline(yintercept = -log10(0.05), lty = 2, alpha = 0.3) +
  labs(title = paste(CONTRAST_ID, "-", STRATUM, "(base)"),
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")
```

---

## Step 4: RUVseq (Removal of Unwanted Variation)

RUVseq estimates latent technical factors (W) that are not captured by known
covariates. It uses empirical negative control genes (genes that are NOT
differentially expressed) to learn these factors.

```r
# ── 4a. Identify empirical negative control genes ────────────────
# Genes with high p-values in the base model are likely not DE
pval_threshold <- cfg$ruvseq$empirical_pval_threshold  # 0.5
neg_control_genes <- rownames(res_base)[
  !is.na(res_base$pvalue) & res_base$pvalue > pval_threshold
]
cat("Negative control genes:", length(neg_control_genes), "\n")
```

```r
# ── 4b. Run RUVg for k = 1 to max_k ────────────────────────────
max_k <- min(cfg$ruvseq$max_latent_factors,
             ncol(counts_filt) - length(active_covs) - 2)
cat("Testing k = 1 to", max_k, "\n")

# Prepare normalized count set for RUVSeq
set <- newSeqExpressionSet(
  counts = as.matrix(counts_filt[rownames(dds), ]),
  phenoData = AnnotatedDataFrame(coldata)
)
set <- betweenLaneNormalization(set, which = "upper")

ruv_results <- list()
for (k in 1:max_k) {
  ruv_results[[k]] <- RUVg(set, neg_control_genes, k = k)
}
```

```r
# ── 4c. Select best k using RLE-ANOVA elbow ─────────────────────
# For each k, compute RLE and ANOVA F-stat (lower = better normalization)
f_stats <- sapply(1:max_k, function(k) {
  norm_counts <- normCounts(ruv_results[[k]])
  # RLE: row median subtracted from each row, then per-column distribution
  rle <- sweep(log2(norm_counts + 1), 1,
               rowMedians(log2(norm_counts + 1)))
  rle_long <- data.frame(
    value  = as.vector(rle),
    sample = rep(colnames(rle), each = nrow(rle))
  )
  anova_fit <- aov(value ~ sample, data = rle_long)
  summary(anova_fit)[[1]][["F value"]][1]
})

# Elbow: point of maximum second derivative
if (max_k > 2) {
  d2 <- diff(diff(f_stats))
  best_k <- which.max(abs(d2)) + 1
} else {
  best_k <- max_k
}
cat("Best k:", best_k, "\n")
```

```r
# ── 4d. Check for disease-associated W factors ──────────────────
# Any W correlated with the contrast variable could absorb real biology
W_mat <- pData(ruv_results[[best_k]])[, paste0("W_", 1:best_k), drop = FALSE]

safe_Ws <- c()
disease_Ws <- c()
for (w in colnames(W_mat)) {
  pval <- anova(
    lm(W_mat[[w]] ~ 1),
    lm(W_mat[[w]] ~ coldata[[contrast_var_safe]])
  )$`Pr(>F)`[2]
  if (pval < 0.05) {
    disease_Ws <- c(disease_Ws, w)
    cat("  ", w, "correlates with", contrast_var, "(p =", round(pval, 4), ") — EXCLUDED\n")
  } else {
    safe_Ws <- c(safe_Ws, w)
    cat("  ", w, "is safe (p =", round(pval, 4), ")\n")
  }
}
```

> **If ALL W factors are disease-associated**: This means the dominant
> technical variation is confounded with the biology. The pipeline falls
> back to the base DESeq2 results (step 3). This is not uncommon when
> the contrast is a strong effect (e.g., disease status).
>
>   - The base model results are still valid — just without RUV correction
>   - Consider adding more covariates to the model to capture batch effects

---

## Step 5: Final DESeq2 with Safe W Factors

```r
if (length(safe_Ws) > 0) {
  # Add safe W factors to coldata
  coldata_ruv <- cbind(coldata, W_mat[, safe_Ws, drop = FALSE])

  # Extended formula: covariates + safe Ws + contrast
  formula_ruv <- paste("~", paste(c(active_covs, safe_Ws, contrast_var_safe),
                                   collapse = " + "))
  cat("RUV formula:", formula_ruv, "\n")

  dds_ruv <- DESeqDataSetFromMatrix(
    countData = counts_filt[rownames(dds), ],
    colData   = coldata_ruv,
    design    = as.formula(formula_ruv)
  )
  dds_ruv <- DESeq(dds_ruv)

  res_ruv <- results(dds_ruv, contrast = c(contrast_var_safe, exp_grp, control_grp))
  res_ruv <- as.data.frame(res_ruv)
  res_ruv$gene <- rownames(res_ruv)

  n_degs_ruv <- sum(res_ruv$padj < fdr_thr &
                    abs(res_ruv$log2FoldChange) > lfc_thr, na.rm = TRUE)
  cat("DEGs (RUV model):", n_degs_ruv, "(vs", n_degs, "base)\n")
} else {
  cat("No safe W factors — using base results as final\n")
  res_ruv <- res_base
}
```

---

## Step 6: fGSEA (Gene Set Enrichment)

```r
# ── 6a. Load gene sets ──────────────────────────────────────────
gmt_file <- cfg$fgsea$gene_sets_file
pathways <- gmtPathways(gmt_file)
cat("Loaded", length(pathways), "pathways\n")

# Filter by size
sizes <- sapply(pathways, length)
pathways <- pathways[sizes >= cfg$fgsea$min_gene_set_size &
                     sizes <= cfg$fgsea$max_gene_set_size]
cat("After size filter:", length(pathways), "pathways\n")
```

```r
# ── 6b. Build ranked gene list ──────────────────────────────────
# Use the best available results (RUV if available, else base)
res_final <- if (exists("res_ruv")) res_ruv else res_base

# Exclude ribosomal + mitochondrial genes
exclude_patterns <- cfg$fgsea$exclude_gene_patterns
keep_genes <- !grepl(paste(exclude_patterns, collapse = "|"), res_final$gene)
res_for_gsea <- res_final[keep_genes, ]

# Rank by Wald statistic (default)
ranks <- setNames(res_for_gsea$stat, res_for_gsea$gene)
ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
cat("Genes ranked:", length(ranks), "\n")
```

```r
# ── 6c. Run fGSEA ───────────────────────────────────────────────
fgsea_res <- fgsea(pathways = pathways, stats = ranks)
fgsea_res <- fgsea_res[order(padj)]

n_sig <- sum(fgsea_res$padj < cfg$fgsea$fdr_threshold, na.rm = TRUE)
cat("Significant pathways:", n_sig, "\n")
head(fgsea_res[, .(pathway, NES, padj, size)], 10)
```

---

## Summary

```r
cat("\n", strrep("=", 60), "\n")
cat(" Summary:", CONTRAST_ID, "—", STRATUM, "\n")
cat(strrep("=", 60), "\n")
cat("  Samples:          ", ncol(counts_filt), "\n")
cat("  Genes (filtered): ", nrow(dds), "\n")
cat("  DEGs (base):      ", n_degs, "\n")
if (exists("n_degs_ruv")) {
  cat("  DEGs (RUV):       ", n_degs_ruv, "\n")
  cat("  RUVseq best k:   ", best_k, "\n")
  cat("  Safe W factors:   ", paste(safe_Ws, collapse = ", "), "\n")
}
cat("  Pathways (sig):   ", n_sig, "\n")
```

---

## What the Pipeline Does Differently

When run via Snakemake (`snakemake --configfile config/config_traits.yaml -j8`),
the pipeline:

- Loops over ALL contrasts x ALL strata (e.g., 15 contrasts x 11 cell types)
- Parallelizes independent jobs across cores
- Writes structured logs, skip sentinels, and preflight reports
- Aggregates results per contrast (step 08) and across the full pipeline (step 09)
- Produces status heatmaps showing which contrast x stratum combinations succeeded

The vignette above covers the exact same statistical steps for a single run.
