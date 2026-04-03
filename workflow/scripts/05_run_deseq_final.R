# =============================================================================
# 05_run_deseq_final.R
# EasyDE - Final DESeq2 with RUVseq Correction
#
# PURPOSE:
#   Runs DESeq2 again using the RUVseq-corrected coldata (with W factors).
#   The formula is extended to include the safe W factors identified in
#   04_run_ruvseq.R (those that do NOT correlate with the contrast variable).
#   Produces the final corrected DE results.
#
# INPUTS:  (all from intermediates/)
#   - coldata.csv                original coldata from 02
#   - counts_filtered.csv        filtered counts from 02
#   - name_map.csv               variable name mapping from 02
#   - ruvseq_best_coldata.csv    coldata + W factors from 04
#   - ruvseq_summary.csv         best_k + safe_Ws from 04
#
# OUTPUTS: (all to intermediates/ and plots/)
#   - results_ruv.tsv            Full DESeq2 results after RUV correction
#   - results_ruv_shrunk.tsv     LFC-shrunk results (if enabled)
#   - plots/volcano_ruv.pdf      Volcano plot
#
# USAGE:
#   Rscript workflow/scripts/05_run_deseq_final.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND \
#       --stratum Beta
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(DESeq2)
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

#' Build the RUV-extended DESeq2 formula
#'
#' Appends safe W factors to the original covariate list.
#' Formula structure: ~ covariate_1 + ... + W_1 + ... + contrast_var
#' The contrast variable is always last (standard DESeq2 convention).
#'
#' @param covariates   Character vector of sanitized covariate names (active)
#' @param safe_ws      Character vector of safe W factor names (e.g. "W_1", "W_2")
#' @param contrast_var Sanitized contrast variable name
#' @return Formula string
build_ruv_formula <- function(covariates, safe_ws, contrast_var) {
    terms <- c(covariates, safe_ws, contrast_var)
    paste0("~ ", paste(terms, collapse = " + "))
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

    logger$info("start", "script=05_run_deseq_final")

    # Define inter_dir and finals_dir early and register skip sentinels
    inter_dir  <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    finals_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "finals")
    register_skip_sentinel(finals_dir, c("results_ruv.tsv"))
    register_skip_sentinel(inter_dir,  c("model_info_ruv.csv"))

    # --- Sentinel fires on any skip ---
    on.exit({
        for (.f in c("results_ruv.tsv")) {
            .p <- file.path(finals_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
        for (.f in c("model_info_ruv.csv")) {
            .p <- file.path(inter_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
    }, add = TRUE)



    # --- Load inputs ---
    # Bail out gracefully if upstream step was skipped
    if (is_skip_sentinel(file.path(inter_dir, "ruvseq_best_coldata.csv"))) {
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

    coldata_ruv <- try_logged(
        read.csv(file.path(inter_dir, "ruvseq_best_coldata.csv"), check.names = FALSE, row.names = NULL),
        logger = logger, step = "load_coldata", success = "ruvseq_best_coldata loaded"
    )
    if (is.null(coldata_ruv)) return(invisible(NULL))
    rownames(coldata_ruv) <- coldata_ruv[[1]]

    counts_df <- try_logged(
        fread(file.path(inter_dir, "counts_filtered.csv")) %>% as.data.frame(),
        logger = logger, step = "load_counts", success = "counts loaded"
    )
    if (is.null(counts_df)) return(invisible(NULL))
    gene_names         <- counts_df[[1]]
    counts_mat         <- as.matrix(counts_df[, -1])
    storage.mode(counts_mat) <- "integer"
    rownames(counts_mat)     <- gene_names

    name_map_df <- fread(file.path(inter_dir, "name_map.csv")) %>% as.data.frame()
    name_map    <- setNames(name_map_df$sanitized, name_map_df$original)

    ruv_summary <- fread(file.path(inter_dir, "ruvseq_summary.csv")) %>% as.data.frame()

    # --- Resolve names ---
    contrast_var     <- trimws(contrast_row$contrast_var)
    control_grp      <- trimws(contrast_row$control_grp)
    experimental_grp <- trimws(contrast_row$experimental_grp)
    covariates_raw   <- parse_pipe_field(contrast_row$covariates)

    safe_contrast <- name_map[[contrast_var]]
    safe_covs     <- sanitize_names(covariates_raw)
    safe_covs     <- safe_covs[safe_covs %in% colnames(coldata_ruv)]

    # Re-coerce factors (CSV strips metadata) and drop constant covariates
    cad <- coerce_and_drop_covariates(coldata_ruv, safe_contrast, safe_covs, logger)
    coldata_ruv <- cad$coldata
    safe_covs   <- cad$covariates

    # Parse safe Ws from the ruvseq summary
    safe_ws_str <- ruv_summary$safe_Ws[1]
    safe_ws     <- if (!is.na(safe_ws_str) && trimws(safe_ws_str) != "") {
        trimws(strsplit(safe_ws_str, "\\|")[[1]])
    } else {
        character(0)
    }

    if (length(safe_ws) == 0) {
        logger$warn("deseq_final",
            "no safe W factors - running final DESeq2 without RUV correction (same as base)")
    }

    # --- Build formula and run DESeq2 ---
    formula_str <- build_ruv_formula(safe_covs, safe_ws, safe_contrast)
    logger$info("deseq_final", sprintf("formula: %s", formula_str))

    counts_mat <- counts_mat[, rownames(coldata_ruv), drop = FALSE]

    # --- Check preflight for zero-inflation flag ---
    preflight <- read_preflight_report(inter_dir)
    use_poscounts <- !is.null(preflight) &&
        "zero_inflated" %in% colnames(preflight) &&
        isTRUE(preflight$zero_inflated[1])

    dds_ruv <- try_logged({
        dds <- suppressWarnings(suppressMessages(DESeq2::DESeqDataSetFromMatrix(
            countData = counts_mat,
            colData   = coldata_ruv,
            design    = as.formula(formula_str)
        )))
        if (use_poscounts) {
            logger$info("deseq_final", "zero-inflated counts — using sfType='poscounts'")
            suppressMessages(DESeq2::DESeq(dds, sfType = "poscounts"))
        } else {
            suppressMessages(DESeq2::DESeq(dds))
        }
    },
    logger  = logger,
    step    = "deseq_final",
    success = "DESeq2 (RUV) completed"
    )
    if (is.null(dds_ruv)) return(invisible(NULL))

    # --- Extract results (shared function from stats_utils.R) ---
    deseq_out <- try_logged(
        extract_deseq_results(
            dds              = dds_ruv,
            contrast_var     = safe_contrast,
            control_grp      = control_grp,
            experimental_grp = experimental_grp,
            lfc_shrinkage    = isTRUE(cfg$deseq2$lfc_shrinkage),
            fdr_threshold    = cfg$deseq2$fdr_threshold,
            l2fc_threshold   = cfg$deseq2$l2fc_threshold,
            logger           = logger
        ),
        logger  = logger,
        step    = "extract_results",
        success = "results extracted"
    )
    if (is.null(deseq_out)) return(invisible(NULL))

    result_df     <- deseq_out$results
    n_degs        <- deseq_out$n_degs
    shrinkage_ran <- deseq_out$shrinkage_ran

    logger$info("deseq_final",
        sprintf("DEGs (FDR<%.2f, |LFC|>%.1f): %d",
                cfg$deseq2$fdr_threshold, cfg$deseq2$l2fc_threshold, n_degs))

    # --- Write outputs ---
    dir.create(finals_dir, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(cfg$write_outputs$ruvseq_results)) {
        write_result(result_df, file.path(finals_dir, "results_ruv.tsv"))
        logger$info("write_outputs",
            sprintf("results_ruv.tsv written (shrunk_col=%s)",
                    ifelse(shrinkage_ran, "yes", "no")))
    }

    saveRDS(dds_ruv, file.path(inter_dir, "dds_ruv.rds"))
    logger$info("write_outputs", "dds_ruv.rds saved")

    model_info_ruv <- data.frame(
        formula = formula_str,
        n_degs  = n_degs,
        stringsAsFactors = FALSE
    )
    write_result(model_info_ruv, file.path(inter_dir, "model_info_ruv.csv"), sep = ",")
    logger$info("write_outputs", "model_info_ruv.csv written")

    if (isTRUE(cfg$write_outputs$volcano_plots)) {
        plot_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "plots")
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

        plot_df <- result_df
        if (shrinkage_ran) {
            plot_df$log2FoldChange <- plot_df$log2FoldChange_shrunk
        }

        p <- plot_volcano(
            deseq_df = plot_df,
            fdr      = cfg$deseq2$fdr_threshold,
            title    = formula_str,
            subtitle = sprintf("DEGs: %d | %s - %s vs %s (RUV corrected)",
                               n_degs, stratum, experimental_grp, control_grp)
        )
        pdf(file.path(plot_dir, "volcano_ruv.pdf"), width = 8, height = 6)
        print(p)
        dev.off()
        logger$info("write_outputs", "volcano_ruv.pdf written")
    }

    logger$info("complete",
        sprintf("n_degs_ruv=%d | formula=%s", n_degs, formula_str))

    invisible(list(
        dds_ruv       = dds_ruv,
        results       = result_df,
        shrinkage_ran = shrinkage_ran,
        n_degs        = n_degs,
        formula       = formula_str
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