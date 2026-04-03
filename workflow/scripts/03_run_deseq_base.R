# =============================================================================
# 03_run_deseq_base.R
# EasyDE - Base DESeq2 Analysis
#
# PURPOSE:
#   Runs the initial DESeq2 model on the filtered coldata and count matrix
#   produced by 02_prepare_coldata.R. This is the "baseline" model - no
#   RUVseq correction yet. Results are used directly if ruvseq is disabled,
#   and as the empirical gene set input for RUVseq if it is enabled.
#
# INPUTS:
#   - results/{contrast_id}/{stratum}/intermediates/coldata.csv
#   - results/{contrast_id}/{stratum}/intermediates/counts_filtered.csv
#   - results/{contrast_id}/{stratum}/intermediates/name_map.csv
#
# OUTPUTS:
#   - intermediates/dds_base.rds             DESeqDataSet object
#   - intermediates/results_base.tsv         Full DESeq2 results table
#   - intermediates/results_base_shrunk.tsv  LFC-shrunk results (if enabled)
#   - plots/volcano_base.pdf                 Volcano plot
#
# USAGE:
#   Rscript workflow/scripts/03_run_deseq_base.R \
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

#' Build the DESeq2 model formula string
#'
#' Constructs the formula from the active covariates and contrast variable.
#' Covariates that were dropped (single-value) are excluded.
#' Format: ~ covariate_1 + covariate_2 + contrast_var
#'
#' @param covariates   Character vector of active sanitized covariate names
#' @param contrast_var Sanitized contrast variable name
#' @return Formula string (e.g. "~ Sex + BMI + disease_status")
build_deseq_formula <- function(covariates, contrast_var) {
    terms <- c(covariates, contrast_var)
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

    logger$info("start", "script=03_run_deseq_base")

    # Define inter_dir and finals_dir early and register skip sentinels
    inter_dir  <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    finals_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "finals")
    register_skip_sentinel(finals_dir, c("results_base.tsv"))
    register_skip_sentinel(inter_dir,  c("model_info_base.csv"))

    # --- Sentinel fires on any skip ---
    on.exit({
        for (.f in c("results_base.tsv")) {
            .p <- file.path(finals_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
        for (.f in c("model_info_base.csv")) {
            .p <- file.path(inter_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
    }, add = TRUE)


    # --- Load intermediates from 02 ---

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

    coldata  <- try_logged(
        read.csv(file.path(inter_dir, "coldata.csv"), check.names = FALSE, row.names = NULL),
        logger = logger, step = "load_coldata", success = "coldata loaded"
    )
    if (is.null(coldata)) return(invisible(NULL))
    rownames(coldata) <- coldata[[1]]

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

    # --- Resolve sanitized names ---
    contrast_var     <- trimws(contrast_row$contrast_var)
    control_grp      <- trimws(contrast_row$control_grp)
    experimental_grp <- trimws(contrast_row$experimental_grp)
    covariates_raw   <- parse_pipe_field(contrast_row$covariates)

    safe_contrast <- name_map[[contrast_var]]
    safe_covs     <- sanitize_names(covariates_raw)

    # Only keep covariates that survived into coldata (not dropped for single value)
    safe_covs <- safe_covs[safe_covs %in% colnames(coldata)]

    # Re-coerce factors (CSV strips metadata) and drop constant covariates
    cad <- coerce_and_drop_covariates(coldata, safe_contrast, safe_covs, logger)
    coldata   <- cad$coldata
    safe_covs <- cad$covariates

    # --- Build formula and run DESeq2 ---
    formula_str <- build_deseq_formula(safe_covs, safe_contrast)
    logger$info("deseq_base", sprintf("formula: %s", formula_str))

    # --- Check preflight for zero-inflation flag ---
    preflight <- read_preflight_report(inter_dir)
    use_poscounts <- !is.null(preflight) &&
        "zero_inflated" %in% colnames(preflight) &&
        isTRUE(preflight$zero_inflated[1])

    # Align column order
    counts_mat <- counts_mat[, rownames(coldata), drop = FALSE]

    dds <- try_logged({
        dds <- suppressWarnings(suppressMessages(DESeq2::DESeqDataSetFromMatrix(
            countData = counts_mat,
            colData   = coldata,
            design    = as.formula(formula_str)
        )))
        if (use_poscounts) {
            logger$info("deseq_base", "zero-inflated counts — using sfType='poscounts'")
            suppressMessages(DESeq2::DESeq(dds, sfType = "poscounts"))
        } else {
            suppressMessages(DESeq2::DESeq(dds))
        }
    },
    logger  = logger,
    step    = "deseq_base",
    success = "DESeq2 completed"
    )
    if (is.null(dds)) return(invisible(NULL))

    logger$info("deseq_base", "extracting results")

    # --- Extract results ---
    deseq_out <- try_logged(
        extract_deseq_results(
            dds              = dds,
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

    logger$info("deseq_base",
        sprintf("DEGs (FDR<%.2f, |LFC|>%.1f): %d",
                cfg$deseq2$fdr_threshold, cfg$deseq2$l2fc_threshold,
                deseq_out$n_degs))
    if (deseq_out$shrinkage_ran) {
        logger$info("deseq_base", "log2FoldChange_shrunk column added to results")
    }

    # --- Write outputs ---
    out_dir   <- inter_dir
    plot_dir  <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "plots")
    dir.create(finals_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(cfg$write_outputs$deseq_base_results)) {
        write_result(deseq_out$results,
                     file.path(finals_dir, "results_base.tsv"))
        logger$info("write_outputs",
            sprintf("results_base.tsv written (shrunk_col=%s)",
                    ifelse(deseq_out$shrinkage_ran, "yes", "no")))
    }

    # Save dds object for downstream RUVseq step
    saveRDS(dds, file.path(out_dir, "dds_base.rds"))
    logger$info("write_outputs", "dds_base.rds saved")

    # Write lightweight model info for downstream aggregation (avoids loading .rds)
    model_info_base <- data.frame(
        formula = formula_str,
        n_degs  = deseq_out$n_degs,
        stringsAsFactors = FALSE
    )
    write_result(model_info_base, file.path(out_dir, "model_info_base.csv"), sep = ",")
    logger$info("write_outputs", "model_info_base.csv written")

    if (isTRUE(cfg$write_outputs$volcano_plots)) {
        # For volcano: prefer shrunk LFC for x-axis if available, else regular LFC
        plot_df <- deseq_out$results
        if (deseq_out$shrinkage_ran) {
            plot_df$log2FoldChange <- plot_df$log2FoldChange_shrunk
        }

        p <- plot_volcano(
            deseq_df = plot_df,
            fdr      = cfg$deseq2$fdr_threshold,
            title    = formula_str,
            subtitle = sprintf("DEGs: %d | %s - %s vs %s",
                               deseq_out$n_degs, stratum, experimental_grp, control_grp)
        )
        pdf(file.path(plot_dir, "volcano_base.pdf"), width = 8, height = 6)
        print(p)
        dev.off()
        logger$info("write_outputs", "volcano_base.pdf written")
    }

    logger$info("complete",
        sprintf("n_degs_base=%d | formula=%s", deseq_out$n_degs, formula_str))

    invisible(list(dds = dds, results = deseq_out))
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