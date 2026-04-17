# =============================================================================
# stats_utils.R
# Small statistical helper functions
#
# These are pure utility functions with no side effects.
# =============================================================================


#' Convert a p-value to an asterisk significance label
#'
#' Used for annotating correlation plots with significance stars.
#'
#' @param x Numeric p-value (or vector of p-values)
#' @return Character string: "", "*", "**", or "***"
pvalue_to_stars <- function(x) {
    dplyr::case_when(
        x >= 0.05             ~ "",
        x < 0.05 & x >= 0.01 ~ "*",
        x < 0.01 & x >= 0.001 ~ "**",
        x < 0.001             ~ "***"
    )
}


#' Correlate RUVseq latent variables against known metadata variables
#'
#' For each W x covariate pair, tests whether the covariate explains
#' significant variation in the latent variable. Uses the omnibus F-test
#' (ANOVA: full model vs intercept-only) so that multi-level factors are
#' tested as a whole rather than picking one arbitrary contrast.
#'
#' For numeric variables and 2-level factors, the F-test is mathematically
#' equivalent to the t-test on the single coefficient (F = t^2, same p-value).
#' For 3+-level factors, it tests all levels jointly — avoiding the old bug
#' where only one arbitrary contrast (row 2 of the coefficient table) was tested.
#'
#' Returns t-statistics (for heatmap display) and omnibus F-test p-values
#' (for exclusion decisions) as matrices (covariates x latent vars).
#'
#' @param meta         data.frame of sample metadata (must contain latent var columns W_1, W_2, ...)
#' @param latent_names Character vector of latent variable column names (e.g. c("W_1", "W_2"))
#' @param correlate_vars Character vector of metadata column names to correlate against
#' @return Named list with two matrices: $t_mat and $p_mat (covariates x latent vars).
#'         $t_mat contains signed t-statistics for heatmap display (for multi-level
#'         factors, uses the coefficient with the largest absolute t-value).
#'         $p_mat contains omnibus F-test p-values for exclusion decisions.
correlate_latent_vars <- function(meta, latent_names, correlate_vars) {

    n_lat  <- length(latent_names)
    n_vars <- length(correlate_vars)

    t_mat <- matrix(nrow = n_vars, ncol = n_lat,
                    dimnames = list(correlate_vars, latent_names))
    p_mat <- matrix(nrow = n_vars, ncol = n_lat,
                    dimnames = list(correlate_vars, latent_names))

    for (v in correlate_vars) {
        for (x in latent_names) {
            fit_full <- lm(as.formula(paste0(x, " ~ ", v)), data = meta, na.action = na.omit)
            fit_null <- lm(as.formula(paste0(x, " ~ 1")),   data = meta, na.action = na.omit)
            coefs    <- summary(fit_full)$coefficients
            # Guard: if the variable was dropped (e.g. only 1 level), record NA
            if (nrow(coefs) < 2) {
                t_mat[v, x] <- NA
                p_mat[v, x] <- NA
            } else {
                # Omnibus F-test p-value for the exclusion decision:
                # tests whether the covariate explains ANY variance in W
                a <- anova(fit_null, fit_full)
                p_mat[v, x] <- a$`Pr(>F)`[2]

                # For heatmap display: use the coefficient with the largest
                # absolute t-statistic (preserves sign/direction for visual)
                non_intercept <- coefs[-1, , drop = FALSE]
                best_row <- which.max(abs(non_intercept[, 3]))
                t_mat[v, x] <- non_intercept[best_row, 3]
            }
        }
    }

    return(list(t_mat = t_mat, p_mat = p_mat))
}


#' Find the best number of RUV latent factors (k) from ANOVA results
#'
#' Uses the "elbow" method: maximizes the second derivative of the
#' F-statistic curve across k values. This finds the point of
#' diminishing returns in variance reduction.
#'
#' @param anova_df data.frame with columns k_num (integer) and f_statistic (numeric)
#' @return Integer: the best k value
find_best_k <- function(anova_df) {

    df     <- anova_df[order(anova_df$k_num), ]
    f_stat <- df$f_statistic

    d1 <- c(NA, diff(f_stat))
    d2 <- c(diff(d1), NA)

    best_k <- df$k_num[which.max(d2)]

    return(best_k)
}


#' Count differentially expressed genes using FDR and LFC thresholds
#'
#' Centralized DEG counting used by 03, 05, and 07. Avoids duplicating
#' the same threshold logic in multiple scripts.
#'
#' @param result_df      data.frame with columns `padj` and `log2FoldChange`
#' @param fdr_threshold  Numeric: maximum adjusted p-value
#' @param l2fc_threshold Numeric: minimum absolute log2 fold change
#' @return Integer count of significant genes
count_degs <- function(result_df, fdr_threshold, l2fc_threshold) {
    sum(!is.na(result_df$padj) &
        result_df$padj < fdr_threshold &
        abs(result_df$log2FoldChange) > l2fc_threshold,
        na.rm = TRUE)
}


#' Sanitize column/variable names for use in R formulas
#'
#' Replaces non-alphanumeric characters with underscores.
#' We use underscores, NOT dots — dots have implicit method-dispatch meaning
#' in R (e.g. print.data.frame) and can cause subtle bugs in model formulas.
#'
#' @param x Character vector of names to sanitize
#' @return Character vector with non-alphanumeric characters replaced by "_"
sanitize_names <- function(x) {
    gsub("[^[:alnum:]]", "_", x)
}


#' Normalize a name for fuzzy matching
#'
#' Strips everything except letters/digits, then lowercases.
#' Used to match stratum names across naming conventions:
#'   "Active Stellate"  ->  "activestellate"
#'   "ActiveStellate"   ->  "activestellate"
#'   "Gamma + Epsilon"  ->  "gammaepsilon"
#'   "Gamma_Epsilon"    ->  "gammaepsilon"
#'
#' Unlike sanitize_names() (which produces R-formula-safe names with underscores),
#' this produces a pure-alnum lowercase key meant ONLY for equality/prefix matching.
#'
#' @param x Character vector of names
#' @return Character vector: lowercased, non-alnum stripped
normalize_name <- function(x) {
    tolower(gsub("[^[:alnum:]]", "", x))
}


#' Canonicalize sample IDs for join safety
#'
#' Mirrors R's make.names() convention: replaces any character not in
#' [A-Za-z0-9_.] with ".", collapses consecutive ".", trims leading/trailing ".".
#' Apply to BOTH count-matrix column names and metadata biosample_id_col values
#' before any intersect/match — the Pseudobulk tool applies make.names() when
#' writing column headers, so this keeps both sides in sync regardless of what
#' special characters ("+", "-", spaces, "/", etc.) appear in treatment IDs.
#'
#' Idempotent: sanitize_sample_id(sanitize_sample_id(x)) == sanitize_sample_id(x)
#'
#' @param x Character vector of sample IDs
#' @return Character vector with special characters replaced by "."
sanitize_sample_id <- function(x) {
    x <- gsub("[^A-Za-z0-9_.]", ".", x)
    x <- gsub("\\.{2,}", ".", x)
    gsub("^\\.|\\.$", "", x)
}


#' Find the metadata cell-count column that corresponds to a stratum name
#'
#' Resolves naming mismatches between how strata appear in contrasts.csv
#' (compact: "ActiveStellate", "Gamma_Epsilon", "Immune") versus how the
#' corresponding cell-count columns appear in sample_metadata.csv
#' (spaced/symbolic: "Active Stellate", "Gamma + Epsilon", "Immune (Macrophages)").
#'
#' Strategy (most-specific first):
#'   1. Exact normalized match  — handles spaces/underscores/symbols
#'   2. Normalized-prefix match — handles parenthetical suffixes like
#'      "Immune" -> "Immune (Macrophages)"
#'      Only triggers if (a) prefix is >= 4 chars and (b) exactly one column
#'      matches (or picks shortest on ties).
#'
#' @param stratum     Stratum name from contrasts.csv (e.g. "Immune")
#' @param col_names   Character vector: all column names of sample_metadata
#' @return Named list: $col_name (matched column or NULL), $method ("exact"|"prefix"|NA),
#'         $ambiguous (logical), $candidates (character vector of all prefix matches, if any)
match_stratum_column <- function(stratum, col_names) {

    clean_target <- normalize_name(stratum)
    clean_cols   <- normalize_name(col_names)

    # 1. Exact normalized match
    exact_idx <- which(clean_cols == clean_target)
    if (length(exact_idx) == 1) {
        return(list(col_name = col_names[exact_idx], method = "exact",
                    ambiguous = FALSE, candidates = col_names[exact_idx]))
    }
    if (length(exact_idx) > 1) {
        return(list(col_name = NULL, method = NA_character_,
                    ambiguous = TRUE, candidates = col_names[exact_idx]))
    }

    # 2. Prefix match (only if target is long enough to be meaningful)
    if (nchar(clean_target) >= 4) {
        prefix_idx <- which(startsWith(clean_cols, clean_target))
        if (length(prefix_idx) == 1) {
            return(list(col_name = col_names[prefix_idx], method = "prefix",
                        ambiguous = FALSE, candidates = col_names[prefix_idx]))
        }
        if (length(prefix_idx) > 1) {
            # Tiebreak: pick the shortest cleaned name (closest match)
            best <- prefix_idx[which.min(nchar(clean_cols[prefix_idx]))]
            return(list(col_name = col_names[best], method = "prefix",
                        ambiguous = TRUE, candidates = col_names[prefix_idx]))
        }
    }

    # No match
    list(col_name = NULL, method = NA_character_,
         ambiguous = FALSE, candidates = character(0))
}


#' Re-coerce factors lost during CSV round-trip and drop constant covariates
#'
#' CSV serialization strips factor metadata. This restores factors for columns
#' that DESeq2 needs as factors, and drops any covariates with <2 observed
#' levels (which cause rank-deficient design matrices).
#'
#' Used by scripts 03, 04, and 05 after loading coldata from CSV.
#'
#' @param coldata       data.frame of sample metadata (modified in place via reference semantics)
#' @param contrast_var  Sanitized contrast variable name
#' @param covariates    Character vector of sanitized covariate names
#' @param logger        Logger instance (optional, for warning messages)
#' @return Named list: $coldata (updated), $covariates (surviving covariates),
#'         $dropped (character vector of dropped covariate names)
coerce_and_drop_covariates <- function(coldata, contrast_var, covariates, logger = NULL) {

    # Re-coerce character/logical columns to factors
    for (v in c(contrast_var, covariates)) {
        if (v %in% colnames(coldata) && (is.character(coldata[[v]]) || is.logical(coldata[[v]]))) {
            coldata[[v]] <- as.factor(coldata[[v]])
        }
    }

    # Drop constant-value covariates
    dropped <- character(0)
    for (v in covariates) {
        if (v %in% colnames(coldata)) {
            n_unique <- if (is.factor(coldata[[v]])) {
                length(levels(droplevels(coldata[[v]])))
            } else {
                length(unique(na.omit(coldata[[v]])))
            }
            if (n_unique < 2) {
                dropped <- c(dropped, v)
            }
        }
    }
    if (length(dropped) > 0) {
        covariates <- setdiff(covariates, dropped)
        if (!is.null(logger)) {
            logger$warn("drop_covariate",
                sprintf("covariates dropped (constant after filtering): %s",
                        paste(dropped, collapse = ", ")))
        }
    }

    list(coldata = coldata, covariates = covariates, dropped = dropped)
}


#' Extract DESeq2 results for a contrast or continuous variable
#'
#' Handles both categorical (contrast vector) and numeric (name-based) contrasts.
#' If lfc_shrinkage is enabled, runs lfcShrink (apeglm) and merges the shrunk
#' estimates as an additional column `log2FoldChange_shrunk` in the same output
#' data.frame, producing ONE file with both LFC estimates side by side.
#' pvalue, padj, and stat are always from the regular (unshrunk) results.
#'
#' Used by scripts 03 and 05.
#'
#' @param dds              DESeqDataSet after running DESeq()
#' @param contrast_var     Sanitized contrast variable name
#' @param control_grp      Control group value (or "" for numeric)
#' @param experimental_grp Experimental group value (or "" for numeric)
#' @param lfc_shrinkage    Logical: attempt lfcShrink and add shrunk column?
#' @param fdr_threshold    Numeric: FDR cutoff for counting DEGs
#' @param l2fc_threshold   Numeric: minimum absolute LFC for counting DEGs
#' @param logger           Logger instance (optional, for shrinkage warnings)
#' @return Named list: $results (data.frame with log2FoldChange and
#'         log2FoldChange_shrunk), $n_degs (integer), $formula (string),
#'         $shrinkage_ran (logical)
extract_deseq_results <- function(dds,
                                  contrast_var,
                                  control_grp,
                                  experimental_grp,
                                  lfc_shrinkage  = TRUE,
                                  fdr_threshold  = 0.05,
                                  l2fc_threshold = 0,
                                  logger         = NULL) {

    is_categorical <- !is.na(control_grp) && control_grp != ""

    if (is_categorical) {
        result <- DESeq2::results(dds,
                                  contrast = c(contrast_var, experimental_grp, control_grp))
    } else {
        result <- DESeq2::results(dds, name = contrast_var)
    }

    result_df        <- as.data.frame(result)
    result_df$gene   <- rownames(result_df)
    result_df        <- result_df[order(result_df$pvalue, na.last = TRUE), ]

    n_degs <- count_degs(result_df, fdr_threshold, l2fc_threshold)

    # LFC shrinkage: merged as a second column, not a separate file.
    shrinkage_ran <- FALSE
    result_df$log2FoldChange_shrunk <- NA_real_

    if (lfc_shrinkage && is_categorical) {
        coef_name <- paste0(
            make.names(contrast_var), "_",
            make.names(experimental_grp), "_vs_",
            make.names(control_grp)
        )
        available_coefs <- DESeq2::resultsNames(dds)
        if (!coef_name %in% available_coefs) {
            match_coef <- available_coefs[grepl(make.names(experimental_grp), available_coefs, fixed = TRUE)]
            if (length(match_coef) == 1) coef_name <- match_coef
        }
        result_shrunk <- tryCatch(
            DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
            error = function(e) {
                if (!is.null(logger)) {
                    logger$warn("lfc_shrinkage",
                        sprintf("lfcShrink failed (coef=%s): %s", coef_name, e$message))
                }
                NULL
            }
        )
        if (!is.null(result_shrunk)) {
            shrunk_df <- as.data.frame(result_shrunk)
            idx <- match(result_df$gene, rownames(shrunk_df))
            result_df$log2FoldChange_shrunk <- shrunk_df$log2FoldChange[idx]
            shrinkage_ran <- TRUE
        }
    }

    list(
        results       = result_df,
        n_degs        = n_degs,
        formula       = deparse(DESeq2::design(dds)),
        shrinkage_ran = shrinkage_ran
    )
}
