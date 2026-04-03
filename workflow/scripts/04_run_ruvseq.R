# =============================================================================
# 04_run_ruvseq.R
# EasyDE — RUVseq Normalization and k Selection
#
# PURPOSE:
#   Uses the base DESeq2 results to identify empirical negative control genes
#   (genes with high p-values — likely not DE), then runs RUVg across k=1..K
#   latent factors. Selects the best k using the elbow of the RLE F-statistic
#   curve (second derivative method). Identifies which W factors correlate with
#   the contrast variable and excludes them from the final formula to avoid
#   regressing out disease signal.
#
# INPUTS:
#   - intermediates/dds_base.rds
#   - finals/results_base.tsv
#   - intermediates/coldata.csv
#   - intermediates/counts_filtered.csv
#   - intermediates/name_map.csv
#
# OUTPUTS:
#   - intermediates/ruvseq_best_coldata.csv   coldata + best W factors
#   - intermediates/ruvseq_summary.csv        k selection + W pruning decisions
#   - plots/ruvseq_diagnostics.pdf            correlation heatmap + RLE elbow plot
#
# USAGE:
#   Rscript workflow/scripts/04_run_ruvseq.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND \
#       --stratum Beta
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(RUVSeq)
    library(EDASeq)
    library(optparse)
})

# Resolve utils directory relative to this script
.script_dir <- tryCatch({
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("--file=", "", file_arg[1])))
    } else {
        dirname(normalizePath(sys.frame(1)$ofile))
    }
}, error = function(e) getwd())

source(file.path(.script_dir, "utils", "io_utils.R"))
source(file.path(.script_dir, "utils", "logging_utils.R"))
source(file.path(.script_dir, "utils", "stats_utils.R"))
source(file.path(.script_dir, "utils", "plot_utils.R"))


# =============================================================================
# Script-specific helpers
# =============================================================================

#' Run RUVg across k = 1..max_k and return all W factor sets
#'
#' Identifies empirical negative control genes (high p-value in base DESeq2),
#' normalizes counts between lanes, then runs RUVg for each k.
#'
#' @param counts_mat          Integer count matrix (genes x samples)
#' @param base_results        data.frame of base DESeq2 results (must have $pvalue)
#' @param max_k               Maximum number of latent factors to test
#' @param empirical_pval_thresh P-value threshold above which genes are considered
#'                             negative controls (not DE)
#' @return Named list indexed by k: each element has $W (factor matrix) and
#'         $norm_counts (normalized count matrix)
run_ruvseq_all_k <- function(counts_mat,
                              base_results,
                              max_k,
                              empirical_pval_thresh = 0.5) {

    # Empirical negative control genes: those with high p-values
    not_sig    <- rownames(base_results)[
        !is.na(base_results$pvalue) & base_results$pvalue > empirical_pval_thresh
    ]

    seq_set    <- EDASeq::newSeqExpressionSet(as.matrix(counts_mat))
    seq_set    <- EDASeq::betweenLaneNormalization(seq_set, which = "upper")
    empirical  <- rownames(seq_set)[rownames(seq_set) %in% not_sig]

    all_results <- list()
    for (k in seq_len(max_k)) {
        ruv_set          <- RUVSeq::RUVg(seq_set, empirical, k = k)
        all_results[[k]] <- list(
            W           = pData(ruv_set),
            norm_counts = ruv_set@assayData$normalizedCounts
        )
    }

    return(all_results)
}


#' Compute RLE-based ANOVA F-statistic for each k
#'
#' For each k, calculates the Relative Log Expression (RLE) matrix,
#' then runs ANOVA (RLE value ~ sample) to get an F-statistic measuring
#' between-sample variance. Lower F = better normalization.
#'
#' @param ruvseq_results  Output of run_ruvseq_all_k()
#' @param biosample_id_col Column name for sample IDs
#' @param contrast_id     For labeling the output data.frame
#' @param stratum         For labeling the output data.frame
#' @return data.frame with columns: k, k_num, f_statistic, residual_variance
compute_rle_anova <- function(ruvseq_results,
                               biosample_id_col,
                               contrast_id,
                               stratum) {

    anova_rows <- lapply(seq_along(ruvseq_results), function(k) {

        norm_counts <- ruvseq_results[[k]]$norm_counts
        W_pheno     <- ruvseq_results[[k]]$W %>%
                           tibble::rownames_to_column(biosample_id_col)

        # RLE: subtract per-gene median from log counts
        log_counts  <- log(norm_counts + 1)
        gene_meds   <- apply(log_counts, 1, median)
        rle_mat     <- apply(log_counts, 2, function(x) x - gene_meds)

        rle_long               <- t(rle_mat) %>% as.data.frame()
        gene_names             <- colnames(rle_long)
        rle_long[[biosample_id_col]] <- rownames(rle_long)
        rle_long               <- rle_long %>%
                                      merge(W_pheno, by = biosample_id_col) %>%
                                      pivot_longer(cols = all_of(gene_names))

        fit            <- aov(as.formula(paste("value ~", biosample_id_col)),
                              data = rle_long)
        anova_summary  <- summary(fit)[[1]]
        f_stat         <- anova_summary[biosample_id_col, "F value"]
        res_var        <- anova_summary["Residuals", "Mean Sq"]

        data.frame(
            contrast         = contrast_id,
            stratum          = stratum,
            k                = paste0("k_", k),
            k_num            = k,
            f_statistic      = f_stat,
            residual_variance = res_var,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, anova_rows)
}


#' Identify W factors associated with any known covariate
#'
#' W factors that significantly correlate with ANY variable we are testing
#' for correlation should be excluded from the final model formula:
#'   - Ws correlated with the contrast variable carry disease/trait signal —
#'     including them would suppress real biology.
#'   - Ws correlated with known covariates (Age, BMI, etc.) are redundant
#'     with those covariates and risk introducing multicollinearity.
#'
#' Two uses of correlation results:
#'   HEATMAP DISPLAY: Bonferroni-adjusted p-values across all variable × W
#'     combinations, returned as $p_mat_adj for the diagnostic plot.
#'   EXCLUSION DECISION: Two-tier strategy.
#'     (a) CONTRAST VARIABLE — raw p < 0.05. We use the most liberal threshold
#'         here because Ws that carry even a hint of disease/trait signal MUST
#'         be removed; a false-positive exclusion only costs statistical power,
#'         while a false-negative lets the W absorb real biology.
#'     (b) OTHER COVARIATES — BH-adjusted (per-W) p < 0.05. For each W, its
#'         p-values against the non-contrast variables are BH-corrected. The W
#'         is excluded if ANY adjusted p < 0.05. This controls per-W FDR
#'         without over-correcting across independent W keep/exclude decisions.
#'
#' @param coldata_with_W   data.frame of coldata merged with W factors
#' @param latent_names     Character vector of W column names (e.g. "W_1", "W_2")
#' @param latent_corr_vars Character vector of variables to correlate against
#' @param contrast_var     Sanitized contrast variable name (used for raw-p exclusion)
#' @param biosample_id_col Column name for sample IDs
#' @return Named list: $t_mat, $p_mat_display (BH-adjusted, for heatmap),
#'         $disease_associated_Ws (excluded from formula), $safe_Ws
identify_disease_ws <- function(coldata_with_W,
                                 latent_names,
                                 latent_corr_vars,
                                 contrast_var,
                                 biosample_id_col) {

    # Run correlations — returns raw p-values (omnibus F-test)
    corr_result <- correlate_latent_vars(
        meta           = coldata_with_W,
        latent_names   = latent_names,
        correlate_vars = latent_corr_vars
    )

    # --- Build the display p-value matrix for the diagnostic heatmap ---
    # BH-adjusted per-W across ALL variables (informational only).
    # Stars on the heatmap help users see which Ws correlate with what,
    # but only the contrast variable drives exclusion decisions.
    p_mat_raw <- corr_result$p_mat
    p_mat_display <- p_mat_raw
    for (w in colnames(p_mat_raw)) {
        pvals_raw <- p_mat_raw[, w, drop = TRUE]
        valid     <- !is.na(pvals_raw)
        pvals_adj <- rep(NA_real_, length(pvals_raw))
        if (sum(valid) > 0) {
            pvals_adj[valid] <- p.adjust(pvals_raw[valid], method = "BH")
        }
        p_mat_display[, w] <- pvals_adj
    }

    # --- W exclusion: contrast variable only ---
    # Exclude any W with raw p < 0.05 for the contrast variable.
    # A false-positive exclusion only costs power; a false-negative
    # lets the W absorb real biological signal.
    # Covariate correlations are shown on the heatmap for information
    # but do NOT drive exclusion — covariates are already in the model.
    associated_ws <- character(0)
    if (contrast_var %in% rownames(p_mat_raw)) {
        for (w in colnames(p_mat_raw)) {
            pval <- p_mat_raw[contrast_var, w]
            if (!is.na(pval) && pval < 0.05) {
                associated_ws <- c(associated_ws, w)
            }
        }
    }

    safe_ws <- latent_names[!latent_names %in% associated_ws]

    list(
        t_mat                = corr_result$t_mat,
        p_mat_display        = p_mat_display,
        disease_associated_Ws = associated_ws,
        safe_Ws              = safe_ws
    )
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id, stratum) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df
    contrast_row <- contrasts_df[contrasts_df$contrast_id == contrast_id, ]

    stratum_safe <- sanitize_names(stratum)
    log_file     <- file.path(cfg$outputs$logs_dir, contrast_id,
                               paste0(stratum_safe, ".log"))
    logger       <- make_logger(contrast_id = contrast_id,
                                biosample   = stratum,
                                log_file    = log_file,
                                level       = cfg$logging$level)

    logger$info("start", "script=04_run_ruvseq")

    # Define inter_dir early and register skip sentinel
    inter_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    register_skip_sentinel(inter_dir, c("ruvseq_best_coldata.csv", "ruvseq_summary.csv"))

    # --- Define inter_dir early so sentinel fires on any skip ---
    on.exit({
        for (.f in c("ruvseq_best_coldata.csv", "ruvseq_summary.csv")) {
            .p <- file.path(inter_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
    }, add = TRUE)

    # --- Paired mode: skip RUVSeq entirely ---
    # In paired analysis the donor is already a blocking factor in the DESeq2
    # formula, making RUVSeq latent factors redundant and potentially harmful
    # (they can absorb genuine biological signal, reducing DEG counts).
    if (isTRUE(cfg$steps$paired_analysis)) {
        logger$skip("paired_skip",
            "paired_analysis enabled — RUVSeq skipped (donor already a blocking factor)")
        return(invisible(NULL))
    }

    # --- Load inputs ---
    # Bail out gracefully if upstream step was skipped
    if (is_skip_sentinel(file.path(inter_dir, "coldata.csv"))) {
        logger$skip("upstream_skip", "upstream step was skipped — propagating skip")
        return(invisible(NULL))
    }

    # Check preflight status from 02
    preflight_status <- read_preflight_status(inter_dir)
    if (preflight_status == "fail") {
        logger$skip("preflight_fail", "preflight validation failed in step 02 — propagating skip")
        return(invisible(NULL))
    }
    if (preflight_status == "warn") {
        logger$warn("preflight_status", "preflight reported warnings — proceeding with caution")
    }

    coldata <- read.csv(file.path(inter_dir, "coldata.csv"), check.names = FALSE, row.names = NULL)
    rownames(coldata) <- coldata[[1]]

    counts_df          <- fread(file.path(inter_dir, "counts_filtered.csv")) %>% as.data.frame()
    gene_names         <- counts_df[[1]]
    counts_mat         <- as.matrix(counts_df[, -1])
    storage.mode(counts_mat) <- "integer"
    rownames(counts_mat)     <- gene_names

    finals_dir   <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "finals")
    base_results <- fread(file.path(finals_dir, "results_base.tsv")) %>% as.data.frame()
    rownames(base_results) <- base_results$gene

    name_map_df  <- fread(file.path(inter_dir, "name_map.csv")) %>% as.data.frame()
    name_map     <- setNames(name_map_df$sanitized, name_map_df$original)

    # --- Resolve names ---
    contrast_var       <- trimws(contrast_row$contrast_var)
    biosample_id_col   <- trimws(contrast_row$biosample_id_col)
    covariates_raw     <- parse_pipe_field(contrast_row$covariates)
    latent_corr_raw    <- parse_pipe_field(contrast_row$latent_corr_vars)

    safe_contrast      <- name_map[[contrast_var]]
    safe_covs          <- sanitize_names(covariates_raw)
    safe_covs          <- safe_covs[safe_covs %in% colnames(coldata)]

    # Re-coerce factors (CSV strips metadata) and drop constant covariates
    cad <- coerce_and_drop_covariates(coldata, safe_contrast, safe_covs, logger)
    coldata   <- cad$coldata
    safe_covs <- cad$covariates

    # Sanitize and drop constant latent_corr_vars (reuse same function)
    safe_corr_vars <- sanitize_names(latent_corr_raw)
    safe_corr_vars <- safe_corr_vars[safe_corr_vars %in% colnames(coldata)]
    cad_corr <- coerce_and_drop_covariates(coldata, safe_contrast, safe_corr_vars, logger = NULL)
    coldata        <- cad_corr$coldata
    safe_corr_vars <- cad_corr$covariates
    if (length(cad_corr$dropped) > 0) {
        logger$warn("drop_corr_var",
            sprintf("latent_corr_vars dropped (constant after filtering): %s",
                    paste(cad_corr$dropped, collapse = ", ")))
    }

    if (length(safe_corr_vars) == 0) {
        logger$warn("ruvseq_setup",
            "no latent_corr_vars found in coldata — W pruning and heatmap will be skipped")
    }

    # --- Determine max k ---
    n_samples  <- ncol(counts_mat)
    n_covs     <- length(safe_covs)
    max_k      <- min(cfg$ruvseq$max_latent_factors, n_samples - n_covs - 1)

    if (max_k <= 0) {
        logger$skip("ruvseq_max_k",
            sprintf("max_k = %d — not enough samples for RUVseq given %d covariates",
                    max_k, n_covs))
        return(invisible(NULL))
    }

    logger$info("ruvseq_setup", sprintf("max_k=%d | n_samples=%d | n_covariates=%d",
                                         max_k, n_samples, n_covs))

    # --- Run RUVseq ---
    logger$info("ruvseq_run", sprintf("running RUVg for k=1..%d", max_k))
    ruvseq_results <- try_logged(
        run_ruvseq_all_k(
            counts_mat            = counts_mat,
            base_results          = base_results,
            max_k                 = max_k,
            empirical_pval_thresh = cfg$ruvseq$empirical_pval_threshold
        ),
        logger  = logger,
        step    = "ruvseq_run",
        success = sprintf("RUVg completed for k=1..%d", max_k)
    )
    if (is.null(ruvseq_results)) return(invisible(NULL))

    # --- Compute RLE ANOVA to find best k ---
    logger$info("ruvseq_k_select", "computing RLE ANOVA across k values")
    anova_df <- compute_rle_anova(ruvseq_results, biosample_id_col,
                                   contrast_id, stratum)

    best_k <- if (max_k <= 2) max_k else find_best_k(anova_df)
    # find_best_k uses a second-derivative elbow method that requires >= 3 k values.
    # With max_k <= 2, just use max_k directly.
    if (length(best_k) == 0 || is.na(best_k)) best_k <- max_k
    logger$info("ruvseq_k_select", sprintf("best_k = %d", best_k))

    # --- Build best coldata with W factors ---
    best_W       <- ruvseq_results[[best_k]]$W
    latent_names <- colnames(best_W)

    coldata_with_W <- best_W %>%
        tibble::rownames_to_column(biosample_id_col) %>%
        merge(coldata, by = biosample_id_col)
    rownames(coldata_with_W) <- coldata_with_W[[biosample_id_col]]

    # --- Identify disease-associated W factors ---
    ws_result <- if (length(safe_corr_vars) > 0) {
        identify_disease_ws(
            coldata_with_W  = coldata_with_W,
            latent_names    = latent_names,
            latent_corr_vars = safe_corr_vars,
            contrast_var    = safe_contrast,
            biosample_id_col = biosample_id_col
        )
    } else {
        list(t_mat = NULL, p_mat_adj = NULL,
             disease_associated_Ws = character(0),
             safe_Ws = latent_names)
    }

    logger$info("ruvseq_w_pruning",
        sprintf("contrast-associated Ws (excluded from formula): %s",
                if (length(ws_result$disease_associated_Ws) == 0) "none"
                else paste(ws_result$disease_associated_Ws, collapse = ", ")))

    logger$info("ruvseq_w_pruning",
        sprintf("safe Ws (kept in formula): %s",
                if (length(ws_result$safe_Ws) == 0) "none"
                else paste(ws_result$safe_Ws, collapse = ", ")))

    # --- Write outputs ---
    write_result(coldata_with_W,
                 file.path(inter_dir, "ruvseq_best_coldata.csv"), sep = ",")

    ruv_summary <- data.frame(
        contrast_id           = contrast_id,
        stratum               = stratum,
        best_k                = best_k,
        disease_associated_Ws = paste(ws_result$disease_associated_Ws, collapse = "|"),
        safe_Ws               = paste(ws_result$safe_Ws, collapse = "|"),
        stringsAsFactors      = FALSE
    )
    write_result(ruv_summary,
                 file.path(inter_dir, "ruvseq_summary.csv"), sep = ",")

    write_result(anova_df,
                 file.path(inter_dir, "ruvseq_anova.csv"), sep = ",")

    logger$info("write_outputs", "ruvseq_best_coldata.csv, ruvseq_summary.csv written")

    # --- Diagnostic plots ---
    if (isTRUE(cfg$write_outputs$ruvseq_diagnostic_plots)) {

        plot_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "plots")
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
        pdf(file.path(plot_dir, "ruvseq_diagnostics.pdf"), width = 10, height = 7)

        # RLE elbow plot
        k_levels   <- paste0("k_", seq_len(max_k))
        anova_df$k <- factor(anova_df$k, levels = k_levels)
        print(plot_rle_elbow(anova_df, best_k,
                             title = paste(stratum, contrast_id)))

        # Correlation heatmap (only if corr vars available)
        if (!is.null(ws_result$t_mat)) {
            print(plot_latent_correlation(
                t_mat         = ws_result$t_mat,
                p_mat_display = ws_result$p_mat_display,
                disease_ws    = ws_result$disease_associated_Ws,
                contrast_var  = contrast_var,
                title         = paste(stratum, contrast_id, "— latent factor correlations")
            ))
        }

        dev.off()
        logger$info("write_outputs", "ruvseq_diagnostics.pdf written")
    }

    logger$info("complete",
        sprintf("best_k=%d | safe_Ws=%s | disease_Ws=%s",
                best_k,
                paste(ws_result$safe_Ws, collapse = ","),
                paste(ws_result$disease_associated_Ws, collapse = ",")))

    invisible(list(
        ruvseq_results = ruvseq_results,
        best_k         = best_k,
        coldata_with_W = coldata_with_W,
        ws_result      = ws_result,
        anova_df       = anova_df
    ))
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {

    opt_list <- list(
        make_option("--config",   type = "character", default = "config/config.yaml"),
        make_option("--contrast", type = "character"),
        make_option("--stratum",  type = "character")
    )
    opts <- parse_args(OptionParser(option_list = opt_list))
    if (is.null(opts$contrast)) stop("--contrast is required")
    if (is.null(opts$stratum))  stop("--stratum is required")

    main(opts$config, opts$contrast, opts$stratum)
}