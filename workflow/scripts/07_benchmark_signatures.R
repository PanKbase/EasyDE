# =============================================================================
# 07_benchmark_signatures.R
# EasyDE - Benchmark DE results against curated gene signatures
#
# PURPOSE:
#   Tests whether DEGs from base and RUV models recover known cell-type
#   gene signatures. Runs two complementary tests:
#
#     1. FORA (overrepresentation) — are UP or DOWN DEGs enriched for
#        signature genes? Each gene set (base, ruv, ruv_only) is split
#        into upregulated (LFC > 0) and downregulated (LFC < 0) subsets,
#        and FORA is run on each direction separately.
#
#     2. fGSEA (ranking-based) — are signature genes enriched at the
#        top/bottom of the full ranked gene list? Uses the Wald statistic.
#        Run on base and RUV rankings separately.
#
#   Designed for benchmarking with cell-type identity markers, disease
#   signatures, pathway gene sets, and negative controls.
#
# INPUTS:
#   - results/{contrast}/{stratum}/finals/results_base.tsv
#   - results/{contrast}/{stratum}/finals/results_ruv.tsv
#   - resources/benchmarking/master_gene_signatures.tsv (single master file)
#
# OUTPUTS:
#   - results/{contrast}/{stratum}/benchmarking/benchmark_signatures.tsv
#
# USAGE:
#   Rscript workflow/scripts/07_benchmark_signatures.R \
#       --config config/config.yaml \
#       --contrast ND_vs_T1D \
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


# =============================================================================
# Script-specific helpers
# =============================================================================

#' Load gene signatures from the master TSV file and filter to a given cell type
#'
#' Reads the master_gene_signatures.tsv file, filters by tissue (matching the
#' config tissue or pan_tissue entries), then matches signatures to the stratum.
#' Cell type matching normalizes both sides: stratum "Beta" matches cell_type
#' "beta_cell", and "all_cells*" patterns match any stratum.
#'
#' Master file columns: signature_id, tissue, cell_type, control_type, category,
#' class, subcategory, description, genes, num_genes, confidence, source_file
#'
#' @param sig_file   Path to the master gene signatures TSV file
#' @param stratum    Cell type to match (e.g., "Beta", "Alpha")
#' @param tissue     Tissue filter from config (e.g., "Pancreas")
#' @return Named list: $pathways (named list of gene vectors for fGSEA),
#'         $meta (data.frame of signature metadata), or NULL if no matches
load_signatures <- function(sig_file, stratum, tissue) {
    if (!file.exists(sig_file)) return(NULL)

    master <- fread(sig_file, sep = "\t") %>% as.data.frame()

    # Filter by tissue: keep pan_tissue + matching tissue
    master <- master[
        grepl("^pan_tissue$", master$tissue, ignore.case = TRUE) |
        tolower(master$tissue) == tolower(tissue),
    ]

    if (nrow(master) == 0) return(NULL)

    stratum_norm <- tolower(stratum)
    ct_base <- tolower(gsub("_cell$", "", master$cell_type))
    is_all  <- grepl("^all_cells|^all$", master$cell_type, ignore.case = TRUE)
    matched <- master[is_all | ct_base == stratum_norm, ]

    if (nrow(matched) == 0) return(NULL)

    # Deduplicate by description
    matched <- matched[!duplicated(matched$description), ]

    # Parse pipe-separated gene lists
    pathways <- setNames(
        lapply(matched$genes, function(g) {
            genes <- trimws(unlist(strsplit(g, "\\|")))
            genes[genes != ""]
        }),
        matched$description
    )

    meta <- matched[, c("description", "control_type", "category", "class",
                         "subcategory", "confidence", "num_genes", "source_file")]

    list(pathways = pathways, meta = meta)
}


#' Define directional DEG gene sets for FORA testing
#'
#' For each of the 3 gene sets (base, ruv, ruv_only), splits DEGs into
#' upregulated (LFC > 0) and downregulated (LFC < 0) subsets.
#'
#' @param base_res  Base DESeq2 results data.frame (or NULL)
#' @param ruv_res   RUV DESeq2 results data.frame (or NULL)
#' @param fdr       FDR threshold for calling DEGs
#' @return Named list with entries: base_all, base_up, base_down, ruv_all,
#'         ruv_up, ruv_down, ruv_only_all, ruv_only_up, ruv_only_down.
#'         Each is a character vector or NULL.
define_directional_deg_sets <- function(base_res, ruv_res, fdr = 0.05) {

    get_degs_directional <- function(res) {
        if (is.null(res)) return(list(up = NULL, down = NULL, all = NULL))
        sig <- res[!is.na(res$padj) & res$padj < fdr, ]
        list(
            up   = sig$gene[sig$log2FoldChange > 0],
            down = sig$gene[sig$log2FoldChange < 0],
            all  = sig$gene
        )
    }

    base <- get_degs_directional(base_res)
    ruv  <- get_degs_directional(ruv_res)

    # RUV-only: significant in RUV but NOT in base
    ruv_only_all  <- NULL
    ruv_only_up   <- NULL
    ruv_only_down <- NULL
    if (!is.null(ruv$all) && !is.null(base$all)) {
        ruv_only_genes <- setdiff(ruv$all, base$all)
        if (length(ruv_only_genes) > 0 && !is.null(ruv_res)) {
            ruv_only_all  <- ruv_only_genes
            ruv_only_rows <- ruv_res[ruv_res$gene %in% ruv_only_genes, ]
            ruv_only_up   <- ruv_only_rows$gene[ruv_only_rows$log2FoldChange > 0]
            ruv_only_down <- ruv_only_rows$gene[ruv_only_rows$log2FoldChange < 0]
        }
    }

    list(
        base_all      = base$all,
        base_up       = base$up,
        base_down     = base$down,
        ruv_all       = ruv$all,
        ruv_up        = ruv$up,
        ruv_down      = ruv$down,
        ruv_only_all  = ruv_only_all,
        ruv_only_up   = ruv_only_up,
        ruv_only_down = ruv_only_down
    )
}


#' Run FORA (overrepresentation analysis) for a set of genes against signatures
#'
#' Uses fgsea::fora() which implements a hypergeometric test.
#'
#' @param genes     Character vector of query genes (e.g., upregulated DEGs)
#' @param pathways  Named list of signature gene sets
#' @param universe  Character vector of all tested genes
#' @return data.frame with columns: pathway, pval, padj, overlap, size,
#'         overlapGenes. Or NULL if no results.
run_fora <- function(genes, pathways, universe) {
    if (length(genes) == 0 || length(universe) == 0) return(NULL)

    # Intersect genes with universe
    genes_in_universe <- intersect(genes, universe)
    if (length(genes_in_universe) == 0) return(NULL)

    result <- tryCatch(
        fgsea::fora(
            pathways = pathways,
            genes    = genes_in_universe,
            universe = universe,
            minSize  = 1,
            maxSize  = 500
        ),
        error = function(e) NULL
    )
    if (is.null(result) || nrow(result) == 0) return(NULL)

    df <- as.data.frame(result)
    # Flatten overlapGenes list column
    if ("overlapGenes" %in% colnames(df)) {
        df$overlapGenes <- sapply(df$overlapGenes, paste, collapse = "|")
    }
    df
}


#' Build ranked gene list for fGSEA from DESeq2 results
#'
#' @param deseq_results  data.frame with gene and stat columns
#' @return Named numeric vector sorted descending
build_ranks <- function(deseq_results) {
    df <- deseq_results[!is.na(deseq_results$stat), ]
    df <- df[order(-df$stat, df$gene), ]
    setNames(df$stat, df$gene)
}


#' Run fGSEA with custom signatures (small gene sets allowed)
#'
#' @param ranks     Named numeric vector of gene ranks
#' @param pathways  Named list of gene sets
#' @return data.frame of fGSEA results, or NULL on failure
run_fgsea_sigs <- function(ranks, pathways) {
    result <- tryCatch(
        suppressWarnings(
            fgsea::fgseaMultilevel(
                pathways = pathways,
                stats    = ranks,
                eps      = 0.0,
                minSize  = 2,
                maxSize  = 500
            )
        ),
        error = function(e) NULL
    )
    if (is.null(result) || nrow(result) == 0) return(NULL)

    df <- as.data.frame(result)[order(result$pval), ]

    # Flatten leadingEdge list column
    if ("leadingEdge" %in% colnames(df)) {
        df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = "|")
    }
    df
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id, stratum) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df

    stratum_safe <- sanitize_names(stratum)
    log_file     <- file.path(cfg$outputs$logs_dir, contrast_id,
                               paste0(stratum_safe, ".log"))
    logger       <- make_logger(contrast_id = contrast_id,
                                biosample   = stratum,
                                log_file    = log_file,
                                level       = cfg$logging$level)

    logger$info("start", "script=07_benchmark_signatures")

    # --- Paths ---
    finals_dir <- file.path(cfg$outputs$results_dir, contrast_id,
                             stratum_safe, "finals")
    bench_dir  <- file.path(cfg$outputs$results_dir, contrast_id,
                             stratum_safe, "benchmarking")
    dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)

    # Register skip sentinel so Snakemake output contract is satisfied
    register_skip_sentinel(bench_dir, c("benchmark.done"))

    # --- Load DE results ---
    base_res <- NULL
    ruv_res  <- NULL

    base_path <- file.path(finals_dir, "results_base.tsv")
    ruv_path  <- file.path(finals_dir, "results_ruv.tsv")

    if (file.exists(base_path) && !is_skip_sentinel(base_path)) {
        base_res <- fread(base_path) %>% as.data.frame()
        logger$info("load_results", sprintf("base results: %d genes", nrow(base_res)))
    }
    if (file.exists(ruv_path) && !is_skip_sentinel(ruv_path)) {
        ruv_res <- fread(ruv_path) %>% as.data.frame()
        logger$info("load_results", sprintf("RUV results: %d genes", nrow(ruv_res)))
    }

    if (is.null(base_res) && is.null(ruv_res)) {
        logger$skip("no_results", "no base or RUV results found — skipping benchmarking")
        return(invisible(NULL))
    }

    # --- Load signatures from master file ---
    sig_file <- cfg$benchmarking$signature_file
    tissue   <- cfg$benchmarking$tissue %||% "Pancreas"

    if (is.null(sig_file) || sig_file == "") {
        logger$info("no_signatures", "no signature_file configured — skipping")
        return(invisible(NULL))
    }

    sig_data <- load_signatures(sig_file, stratum, tissue)
    if (is.null(sig_data)) {
        logger$info("no_signatures",
            sprintf("no matching signatures for stratum '%s' (tissue=%s) — skipping",
                    stratum, tissue))
        return(invisible(NULL))
    }
    logger$info("load_signatures",
        sprintf("%d signatures loaded for %s", length(sig_data$pathways), stratum))

    # --- Define directional DEG sets ---
    fdr <- cfg$deseq2$fdr_threshold %||% 0.05
    deg_sets <- define_directional_deg_sets(base_res, ruv_res, fdr)

    for (nm in names(deg_sets)) {
        n <- if (is.null(deg_sets[[nm]])) "NA" else length(deg_sets[[nm]])
        logger$info("deg_sets", sprintf("%s: %s genes", nm, n))
    }

    # --- FORA (directional overrepresentation) ---
    fora_rows <- list()

    # Map each directional set to its parent gene_set label and source results
    # "all" direction = combined up+down DEGs (for overall overlap fraction)
    fora_configs <- list(
        list(name = "base_all",      gene_set = "base",     direction = "all",
             genes = deg_sets$base_all,      res = base_res),
        list(name = "base_up",       gene_set = "base",     direction = "up",
             genes = deg_sets$base_up,       res = base_res),
        list(name = "base_down",     gene_set = "base",     direction = "down",
             genes = deg_sets$base_down,     res = base_res),
        list(name = "ruv_all",       gene_set = "ruv",      direction = "all",
             genes = deg_sets$ruv_all,       res = ruv_res),
        list(name = "ruv_up",        gene_set = "ruv",      direction = "up",
             genes = deg_sets$ruv_up,        res = ruv_res),
        list(name = "ruv_down",      gene_set = "ruv",      direction = "down",
             genes = deg_sets$ruv_down,      res = ruv_res),
        list(name = "ruv_only_all",  gene_set = "ruv_only", direction = "all",
             genes = deg_sets$ruv_only_all,  res = ruv_res),
        list(name = "ruv_only_up",   gene_set = "ruv_only", direction = "up",
             genes = deg_sets$ruv_only_up,   res = ruv_res),
        list(name = "ruv_only_down", gene_set = "ruv_only", direction = "down",
             genes = deg_sets$ruv_only_down, res = ruv_res)
    )

    for (fc in fora_configs) {
        if (is.null(fc$genes) || is.null(fc$res)) next

        universe <- fc$res$gene[!is.na(fc$res$padj)]
        fora_out <- run_fora(fc$genes, sig_data$pathways, universe)

        if (!is.null(fora_out) && nrow(fora_out) > 0) {
            fora_out$gene_set  <- fc$gene_set
            fora_out$direction <- fc$direction
            fora_out$n_degs    <- length(fc$genes)
            fora_rows[[length(fora_rows) + 1]] <- fora_out
        }
    }

    fora_df <- if (length(fora_rows) > 0) do.call(rbind, fora_rows) else NULL

    # --- fGSEA (ranking-based enrichment) ---
    fgsea_rows <- list()

    fgsea_configs <- list(
        base = base_res,
        ruv  = ruv_res
    )

    for (gs_name in names(fgsea_configs)) {
        res <- fgsea_configs[[gs_name]]
        if (is.null(res)) next

        ranks     <- build_ranks(res)
        fgsea_out <- run_fgsea_sigs(ranks, sig_data$pathways)

        if (!is.null(fgsea_out) && nrow(fgsea_out) > 0) {
            fgsea_out$gene_set <- gs_name
            fgsea_rows[[length(fgsea_rows) + 1]] <- fgsea_out
        }
    }

    fgsea_df <- if (length(fgsea_rows) > 0) do.call(rbind, fgsea_rows) else NULL

    # --- Assemble output ---

    # Standardize FORA columns
    if (!is.null(fora_df)) {
        fora_out <- fora_df %>%
            left_join(sig_data$meta, by = c("pathway" = "description")) %>%
            transmute(
                contrast       = contrast_id,
                stratum        = stratum,
                signature      = pathway,
                control_type   = control_type,
                category       = category,
                class          = class,
                subcategory    = subcategory,
                confidence     = confidence,
                n_genes_sig    = num_genes,
                gene_set       = gene_set,
                direction      = direction,
                test           = "fora",
                statistic      = overlap / size,  # proportion of signature found
                statistic_name = "overlap_fraction",
                NES            = NA_real_,
                pvalue         = pval,
                padj           = padj,
                n_overlap      = overlap,
                n_sig_in_universe = size,
                n_degs         = n_degs,
                overlap_genes  = overlapGenes
            )
    } else {
        fora_out <- NULL
    }

    # Standardize fGSEA columns
    if (!is.null(fgsea_df)) {
        fgsea_out_df <- fgsea_df %>%
            left_join(sig_data$meta, by = c("pathway" = "description")) %>%
            transmute(
                contrast       = contrast_id,
                stratum        = stratum,
                signature      = pathway,
                control_type   = control_type,
                category       = category,
                class          = class,
                subcategory    = subcategory,
                confidence     = confidence,
                n_genes_sig    = num_genes,
                gene_set       = gene_set,
                direction      = NA_character_,
                test           = "fgsea",
                statistic      = NES,
                statistic_name = "NES",
                NES            = NES,
                pvalue         = pval,
                padj           = padj,
                n_overlap      = size,
                n_sig_in_universe = NA_integer_,
                n_degs         = NA_integer_,
                overlap_genes  = leadingEdge
            )
    } else {
        fgsea_out_df <- NULL
    }

    combined <- rbind(fora_out, fgsea_out_df)

    if (is.null(combined) || nrow(combined) == 0) {
        logger$info("complete", "no benchmark results produced")
        return(invisible(NULL))
    }

    # --- Write output ---
    out_file <- file.path(bench_dir, "benchmark_signatures.tsv")
    fwrite(combined, out_file, sep = "\t")
    logger$info("write_output", sprintf("benchmark_signatures.tsv: %d rows", nrow(combined)))

    # --- Log summary ---
    if (!is.null(fora_out)) {
        n_fora_sig <- sum(!is.na(fora_out$padj) & fora_out$padj < 0.05, na.rm = TRUE)
        logger$info("summary_fora",
            sprintf("FORA tests: %d total, %d significant (padj<0.05)",
                    nrow(fora_out), n_fora_sig))
    }
    if (!is.null(fgsea_out_df)) {
        n_fgsea_sig <- sum(!is.na(fgsea_out_df$padj) & fgsea_out_df$padj < 0.05, na.rm = TRUE)
        logger$info("summary_fgsea",
            sprintf("fGSEA tests: %d total, %d significant (padj<0.05)",
                    nrow(fgsea_out_df), n_fgsea_sig))
    }

    logger$info("complete", sprintf("benchmarking done -- %d total rows written", nrow(combined)))

    # --- Per-stratum benchmark plots ---
    tryCatch({
        source(file.path(.script_dir, "utils", "plot_utils.R"))

        bench_fora <- combined[combined$test == "fora", ]
        has_ct <- "control_type" %in% colnames(bench_fora)

        if (nrow(bench_fora) > 0) {
            bench_collapsed <- collapse_benchmark_by_class(bench_fora)

            # 1. Positive controls
            save_heatmap(plot_fora_heatmap(
                bench_fora,
                control_type_val = if (has_ct) "positive" else NULL,
                categories_filter = if (!has_ct) "cell_identity" else NULL,
                title = sprintf("%s / %s -- Positive Controls", contrast_id, stratum),
                fill_color = "#2E8B57", fill_limit = 0.6
            ), file.path(bench_dir, "benchmark_positive.pdf"), logger)

            # 2. Negative controls
            save_heatmap(plot_fora_heatmap(
                bench_fora,
                control_type_val = if (has_ct) "negative" else NULL,
                categories_filter = if (!has_ct) c("artifact", "subject_confounder") else NULL,
                title = sprintf("%s / %s -- Negative Controls", contrast_id, stratum),
                fill_color = "#D94040", fill_limit = 0.5
            ), file.path(bench_dir, "benchmark_negative.pdf"), logger)

            # 3. Collapsed (class-level)
            if (!is.null(bench_collapsed) && nrow(bench_collapsed) > 0) {
                save_heatmap(plot_fora_heatmap(
                    bench_collapsed,
                    title = sprintf("%s / %s -- Class-Level Collapsed", contrast_id, stratum),
                    fill_color = "#7B3294", fill_limit = 0.3, collapsed = TRUE
                ), file.path(bench_dir, "benchmark_collapsed.pdf"), logger)
            }
        }
    }, error = function(e) {
        logger$warn("benchmark_plots", sprintf("could not generate plots: %s", e$message))
    })

    invisible(combined)
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
    opts <- parse_args(OptionParser(option_list = opt_list),
                       args = commandArgs(trailingOnly = TRUE))
    if (is.null(opts$contrast)) stop("--contrast is required")
    if (is.null(opts$stratum))  stop("--stratum is required")

    main(opts$config, opts$contrast, opts$stratum)
}
