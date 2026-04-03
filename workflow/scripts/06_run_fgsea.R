# =============================================================================
# 06_run_fgsea.R
# EasyDE - Gene Set Enrichment Analysis (fGSEA)
#
# PURPOSE:
#   Runs fGSEA on DE results. Respects the fgsea.run_on config setting:
#     "base"   - runs only on base DESeq2 results
#     "ruvseq" - runs only on RUV-corrected results
#     "both"   - runs on both, writes separate output files for each
#
#   If steps.ruvseq is false, always falls back to base results regardless
#   of run_on setting.
#
#   Before ranking genes, removes ribosomal, mitochondrial, and other
#   user-specified gene lists to avoid pathway inflation from housekeeping genes.
#
# INPUTS:  (from intermediates/)
#   - results_base.tsv / results_base_shrunk.tsv
#   - results_ruv.tsv / results_ruv_shrunk.tsv (if ruvseq enabled)
#   - resources/gsea_files/*.gmt.txt
#
# OUTPUTS: (to intermediates/ and plots/)
#   - fgsea_base_all.tsv / fgsea_ruv_all.tsv
#   - fgsea_base_sig.tsv / fgsea_ruv_sig.tsv
#   - plots/fgsea_base.pdf / plots/fgsea_ruv.pdf
#
# USAGE:
#   Rscript workflow/scripts/06_run_fgsea.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND \
#       --stratum Beta
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(fgsea)
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

#' Load and merge gene exclusion lists into one character vector
#'
#' Reads each CSV file in exclude_gene_lists and extracts gene symbols from
#' the "Approved symbol" column. Also applies regex pattern exclusions.
#' Returns a combined character vector of all genes to exclude.
#'
#' @param exclude_gene_lists  Character vector of CSV file paths
#' @param exclude_patterns    Character vector of regex patterns to exclude
#' @param all_genes           Character vector of all gene names in the results
#'                            (used to apply pattern matching)
#' @return Character vector of gene symbols to exclude
build_exclusion_list <- function(exclude_gene_lists,
                                  exclude_patterns,
                                  all_genes) {

    excluded <- character(0)

    # From gene list CSV files
    for (f in exclude_gene_lists) {
        if (!file.exists(f)) next
        gene_df   <- tryCatch(
            fread(f) %>% as.data.frame(),
            error = function(e) NULL
        )
        if (is.null(gene_df)) next
        if ("Approved symbol" %in% colnames(gene_df)) {
            excluded <- c(excluded, na.omit(gene_df[["Approved symbol"]]))
        }
    }

    # From regex patterns applied to gene names
    for (pattern in exclude_patterns) {
        matched  <- all_genes[grepl(pattern, all_genes)]
        excluded <- c(excluded, matched)
    }

    unique(excluded)
}


#' Build a ranked gene list for fGSEA from DESeq2 results
#'
#' Supports three ranking statistics, configured via fgsea.ranking_stat:
#'
#'   "stat"        - DESeq2 Wald statistic. Signed and standardized; already
#'                   captures direction and significance. Default and original
#'                   pipeline behavior (res$stat).
#'   "lfc_pvalue"  - sign(log2FoldChange) * -log10(pvalue). Weights heavily
#'                   by fold change magnitude.
#'   "lfc_shrunk"  - sign(log2FoldChange_shrunk) * -log10(pvalue). Uses the
#'                   apeglm-shrunk LFC for direction. Falls back to lfc_pvalue
#'                   if log2FoldChange_shrunk column is absent or all NA.
#'
#' Genes in the exclusion list and genes with NA stat values are dropped.
#' Ties broken alphabetically by gene name for reproducibility.
#'
#' @param deseq_results   data.frame with columns: gene, stat, log2FoldChange,
#'                        pvalue, and optionally log2FoldChange_shrunk
#' @param exclude_genes   Character vector of gene symbols to remove
#' @param ranking_stat    One of "stat", "lfc_pvalue", "lfc_shrunk"
#' @return Named numeric vector sorted descending (name = gene symbol)
build_gene_ranks <- function(deseq_results,
                              exclude_genes  = character(0),
                              ranking_stat   = "stat") {

    df <- deseq_results[!deseq_results$gene %in% exclude_genes, ]

    rank_values <- switch(ranking_stat,

        "stat" = {
            df <- df[!is.na(df$stat), ]
            df$stat
        },

        "lfc_pvalue" = {
            df <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue), ]
            sign(df$log2FoldChange) * -log10(df$pvalue + 1e-300)
        },

        "lfc_shrunk" = {
            # Use shrunk LFC column if present and non-missing; otherwise fall back
            has_shrunk <- "log2FoldChange_shrunk" %in% colnames(df) &&
                          !all(is.na(df$log2FoldChange_shrunk))
            if (has_shrunk) {
                df <- df[!is.na(df$log2FoldChange_shrunk) & !is.na(df$pvalue), ]
                sign(df$log2FoldChange_shrunk) * -log10(df$pvalue + 1e-300)
            } else {
                df <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue), ]
                sign(df$log2FoldChange) * -log10(df$pvalue + 1e-300)
            }
        },

        # Unknown value: fall back to stat with a warning
        {
            warning(sprintf("Unknown fgsea.ranking_stat='%s', falling back to 'stat'",
                            ranking_stat))
            df <- df[!is.na(df$stat), ]
            df$stat
        }
    )

    df$rank_stat <- rank_values
    df <- df[order(-df$rank_stat, df$gene), ]
    setNames(df$rank_stat, df$gene)
}


#' Load DE results for fGSEA from the single merged output file
#'
#' Results files now contain both log2FoldChange (MLE) and log2FoldChange_shrunk
#' (NA if shrinkage was not run) in one file. No separate _shrunk.tsv exists.
#'
#' @param inter_dir   Path to the intermediates directory
#' @param suffix      One of "base" or "ruv"
#' @return data.frame of DE results, or NULL if file not found
load_de_results_for_fgsea <- function(finals_dir, suffix) {
    path <- file.path(finals_dir, sprintf("results_%s.tsv", suffix))
    if (!file.exists(path)) return(NULL)
    if (is_skip_sentinel(path)) return(NULL)
    return(fread(path) %>% as.data.frame())
}


#' Run one fGSEA pass (base or ruv) and write outputs
#'
#' Run fGSEA and return all + significant results
#'
#' @param ranks         Named numeric vector of gene ranks
#' @param pathways      Named list of gene sets (from fgsea::gmtPathways)
#' @param min_size      Minimum gene set size
#' @param max_size      Maximum gene set size
#' @param fdr_threshold Adjusted p-value cutoff for significant results
#' @param logger        Logger instance (optional)
#' @param label         Step label for logging
#' @return Named list: $all (data.frame), $sig (data.frame), or NULL on failure
run_fgsea_analysis <- function(ranks,
                                pathways,
                                min_size      = 10,
                                max_size      = 500,
                                fdr_threshold = 0.10,
                                logger        = NULL,
                                label         = "fgsea") {

    result <- tryCatch(
        suppressWarnings(
            fgsea::fgseaMultilevel(
                pathways = pathways,
                stats    = ranks,
                eps      = 0.0,
                minSize  = min_size,
                maxSize  = max_size
            )
        ),
        error = function(e) {
            if (!is.null(logger))
                logger$warn(label, sprintf("fgseaMultilevel failed: %s", e$message))
            NULL
        }
    )

    if (is.null(result)) return(NULL)

    all_df <- as.data.frame(result)[order(as.data.frame(result)$pval), ]
    sig_df <- all_df[!is.na(all_df$padj) & all_df$padj < fdr_threshold, ]

    if (!is.null(logger))
        logger$info(label, sprintf("total=%d | significant=%d (FDR<%.2f)",
                                   nrow(all_df), nrow(sig_df), fdr_threshold))

    # leadingEdge is a list column - collapse to pipe-separated string for writing
    flatten_leading_edge <- function(df) {
        if ("leadingEdge" %in% colnames(df))
            df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = "|")
        df
    }
    list(all = flatten_leading_edge(all_df),
         sig = flatten_leading_edge(sig_df),
         n_sig = nrow(sig_df))
}


#' Promoted from an inner closure to a top-level function so all dependencies
#' are explicit parameters rather than captured from the enclosing environment.
#'
#' @param suffix           "base" or "ruv"
#' @param inter_dir        Path to intermediates directory
#' @param plot_dir         Path to plots directory
#' @param cfg              Full config list
#' @param pathways         Named list of gene sets from gmtPathways()
#' @param static_excluded  Gene symbols excluded from file-based lists (built once)
#' @param stratum          Stratum name (for plot titles)
#' @param contrast_id      Contrast ID (for plot titles)
#' @param logger           Logger instance
#' @return Named list: $n_sig (integer or NA if skipped)
run_one_fgsea_pass <- function(suffix,
                                inter_dir,
                                finals_dir,
                                plot_dir,
                                cfg,
                                pathways,
                                static_excluded,
                                stratum,
                                contrast_id,
                                logger) {

    de_results <- load_de_results_for_fgsea(finals_dir, suffix)

    if (is.null(de_results)) {
        logger$warn(sprintf("fgsea_%s", suffix),
            sprintf("no DE results file found for suffix '%s' - skipping", suffix))
        return(list(n_sig = NA))
    }

    # Pattern-based exclusions re-applied to actual gene names in this result set
    pattern_excluded <- build_exclusion_list(
        exclude_gene_lists = character(0),
        exclude_patterns   = cfg$fgsea$exclude_gene_patterns %||% character(0),
        all_genes          = de_results$gene
    )
    all_excluded <- unique(c(static_excluded, pattern_excluded))

    logger$info(sprintf("fgsea_%s", suffix),
        sprintf("genes excluded: %d | genes remaining: %d",
                length(all_excluded),
                sum(!de_results$gene %in% all_excluded)))

    ranking_stat <- cfg$fgsea$ranking_stat %||% "stat"
    ranks <- build_gene_ranks(
        deseq_results = de_results,
        exclude_genes = all_excluded,
        ranking_stat  = ranking_stat
    )
    logger$info(sprintf("fgsea_%s", suffix),
        sprintf("ranking_stat=%s | genes ranked: %d", ranking_stat, length(ranks)))

    # Guard: all genes excluded — no ranks to analyze
    if (length(ranks) == 0) {
        logger$warn(sprintf("fgsea_%s", suffix),
            "all genes were excluded — no ranked genes available for fGSEA")
        return(list(n_sig = NA))
    }

    result <- run_fgsea_analysis(
        ranks         = ranks,
        pathways      = pathways,
        min_size      = cfg$fgsea$min_gene_set_size,
        max_size      = cfg$fgsea$max_gene_set_size,
        fdr_threshold = cfg$fgsea$fdr_threshold,
        logger        = logger,
        label         = sprintf("fgsea_%s", suffix)
    )
    if (is.null(result)) return(list(n_sig = NA))

    if (isTRUE(cfg$write_outputs$fgsea_results_all)) {
        write_result(result$all,
                     file.path(finals_dir, sprintf("fgsea_%s.tsv", suffix)))
        logger$info("write_outputs",
            sprintf("fgsea_%s.tsv written (%d pathways)", suffix, nrow(result$all)))
    }

    p <- plot_fgsea_barplot(
        fgsea_sig_df = result$sig,
        top_n        = cfg$fgsea$n_top_pathways_plot,
        title        = sprintf("Top pathways - %s %s (%s)", stratum, contrast_id, suffix)
    )
    if (!is.null(p)) {
        pdf(file.path(plot_dir, sprintf("fgsea_%s.pdf", suffix)), width = 12, height = 8)
        print(p)
        dev.off()
        logger$info("write_outputs", sprintf("fgsea_%s.pdf written", suffix))
    }

    list(n_sig = result$n_sig)
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

    logger$info("start", "script=06_run_fgsea")

    # Define inter_dir early and register skip sentinel
    inter_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    register_skip_sentinel(inter_dir, c("fgsea.done"))

    # Check preflight status from 02
    preflight_status <- read_preflight_status(inter_dir)
    if (preflight_status == "fail") {
        logger$skip("preflight_fail", "preflight validation failed in step 02 — propagating skip")
        return(invisible(NULL))
    }
    if (preflight_status == "warn") {
        logger$warn("preflight_status", "preflight reported warnings — proceeding with caution")
    }

    finals_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "finals")
    plot_dir   <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "plots")
    dir.create(finals_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Determine which results to run fGSEA on ---
    run_on        <- cfg$fgsea$run_on
    ruvseq_active <- isTRUE(cfg$steps$ruvseq) && !isTRUE(cfg$steps$paired_analysis)

    # If ruvseq not active (disabled or paired mode), always fall back to base
    if (!ruvseq_active && run_on %in% c("ruvseq", "both")) {
        logger$warn("fgsea_routing",
            sprintf("fgsea.run_on='%s' but RUVSeq inactive (disabled or paired mode) - running on base results only",
                    run_on))
        run_on <- "base"
    }

    run_on_base  <- run_on %in% c("base", "both")
    run_on_ruv   <- run_on %in% c("ruvseq", "both")

    # --- Load gene sets ---
    logger$info("load_gene_sets",
        sprintf("loading gene sets from: %s", cfg$fgsea$gene_sets_file))
    pathways <- try_logged(
        suppressWarnings(fgsea::gmtPathways(cfg$fgsea$gene_sets_file)),
        logger  = logger,
        step    = "load_gene_sets"
    )
    if (is.null(pathways)) return(invisible(NULL))

    logger$info("load_gene_sets", sprintf("%d pathways loaded", length(pathways)))

    # --- Build exclusion list from files (done once; patterns applied per pass) ---
    static_excluded <- build_exclusion_list(
        exclude_gene_lists = cfg$fgsea$exclude_gene_lists %||% character(0),
        exclude_patterns   = character(0),   # patterns applied per-pass against actual genes
        all_genes          = character(0)
    )

    # --- Execute passes ---
    n_sig_base <- if (run_on_base) {
        run_one_fgsea_pass(
            suffix          = "base",
            inter_dir       = inter_dir,
            finals_dir      = finals_dir,
            plot_dir        = plot_dir,
            cfg             = cfg,
            pathways        = pathways,
            static_excluded = static_excluded,
            stratum         = stratum,
            contrast_id     = contrast_id,
            logger          = logger
        )$n_sig
    } else NA

    n_sig_ruv <- if (run_on_ruv) {
        run_one_fgsea_pass(
            suffix          = "ruv",
            inter_dir       = inter_dir,
            finals_dir      = finals_dir,
            plot_dir        = plot_dir,
            cfg             = cfg,
            pathways        = pathways,
            static_excluded = static_excluded,
            stratum         = stratum,
            contrast_id     = contrast_id,
            logger          = logger
        )$n_sig
    } else NA

    logger$info("complete",
        sprintf("n_sig_base=%s | n_sig_ruv=%s",
                ifelse(is.na(n_sig_base), "skipped", n_sig_base),
                ifelse(is.na(n_sig_ruv),  "skipped", n_sig_ruv)))

    invisible(list(n_sig_base = n_sig_base, n_sig_ruv = n_sig_ruv))
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