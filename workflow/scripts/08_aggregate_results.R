# =============================================================================
# 08_aggregate_results.R
# EasyDE -- Results Aggregation and Run Summary
#
# PURPOSE:
#   Runs ONCE per contrast (not per stratum). Walks all stratum subdirectories
#   under results/{contrast_id}/, collects per-stratum outputs, and assembles:
#
#     1. contrast_summary.csv                -- one row per stratum, audit trail
#     2. preflight_summary.csv               -- concatenated preflight reports
#     3. {contrast}_deseq_base.tsv           -- all strata DESeq2 base results
#     4. {contrast}_deseq_ruv.tsv            -- all strata DESeq2 RUV results
#     5. {contrast}_fgsea_base.tsv           -- all strata fGSEA base results
#     6. {contrast}_fgsea_ruv.tsv            -- all strata fGSEA RUV results
#     7. {contrast}_benchmark_signatures.tsv -- all strata benchmark results
#     8. errors.log / warnings.log           -- split consolidated logs
#     9. log_summary.tsv                     -- aggregated flagged log lines
#
# INPUTS:
#   - results/{contrast_id}/{stratum}/finals/*.tsv
#   - results/{contrast_id}/{stratum}/benchmarking/benchmark_signatures.tsv
#   - results/{contrast_id}/{stratum}/intermediates/*.csv
#   - logs/{contrast_id}/*.log
#
# USAGE:
#   Rscript workflow/scripts/08_aggregate_results.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
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
# Script-level helpers
# =============================================================================

#' Safe file reader -- returns NULL if file missing or is a skip sentinel
safe_read <- function(path) {
    if (!file.exists(path)) return(NULL)
    if (is_skip_sentinel(path)) return(NULL)
    tryCatch(fread(path) %>% as.data.frame(), error = function(e) NULL)
}


#' Collect summary statistics for one stratum
#'
#' Reads whatever output files exist -- missing files just become NA values.
#' This means a partially-completed stratum still contributes a row to the
#' summary, showing exactly how far it got before failing.
#'
#' @param inter_dir     Path to {contrast_id}/{stratum}/intermediates/
#' @param finals_dir    Path to {contrast_id}/{stratum}/finals/
#' @param contrast_id   Contrast identifier string
#' @param stratum       Stratum name
#' @param fdr_threshold FDR threshold for counting DEGs and pathways
#' @param l2fc_threshold Minimum |LFC| for counting DEGs
#' @param fgsea_fdr     FDR threshold for counting significant pathways
#' @param log_dir       Path to logs/{contrast_id}/ directory (for error messages)
#' @return Single-row data.frame (a summary row)
collect_stratum_summary <- function(inter_dir,
                                     finals_dir,
                                     contrast_id,
                                     stratum,
                                     fdr_threshold,
                                     l2fc_threshold,
                                     fgsea_fdr,
                                     log_dir = NULL) {

    # --- coldata: sample counts (from intermediates) ---
    coldata      <- safe_read(file.path(inter_dir, "coldata.csv"))
    n_total      <- if (!is.null(coldata)) nrow(coldata) else NA

    # --- base DESeq2 results (from finals) ---
    base_results <- safe_read(file.path(finals_dir, "results_base.tsv"))

    n_degs_base <- if (!is.null(base_results)) {
        sum(!is.na(base_results$padj) &
            base_results$padj < fdr_threshold &
            abs(base_results$log2FoldChange) > l2fc_threshold, na.rm = TRUE)
    } else NA

    # --- model info: formulas (from intermediates) ---
    model_base   <- safe_read(file.path(inter_dir, "model_info_base.csv"))
    base_formula <- if (!is.null(model_base) && "formula" %in% colnames(model_base)) {
        model_base$formula[1]
    } else NA

    model_ruv    <- safe_read(file.path(inter_dir, "model_info_ruv.csv"))
    ruv_formula  <- if (!is.null(model_ruv) && "formula" %in% colnames(model_ruv)) {
        model_ruv$formula[1]
    } else NA

    # --- RUVseq summary (from intermediates) ---
    ruv_summary  <- safe_read(file.path(inter_dir, "ruvseq_summary.csv"))
    best_k       <- if (!is.null(ruv_summary)) ruv_summary$best_k[1]    else NA
    safe_ws      <- if (!is.null(ruv_summary)) ruv_summary$safe_Ws[1]   else NA
    disease_ws   <- if (!is.null(ruv_summary)) ruv_summary$disease_associated_Ws[1] else NA

    # --- RUV DESeq2 results (from finals) ---
    ruv_results  <- safe_read(file.path(finals_dir, "results_ruv.tsv"))

    n_degs_ruv <- if (!is.null(ruv_results)) {
        sum(!is.na(ruv_results$padj) &
            ruv_results$padj < fdr_threshold &
            abs(ruv_results$log2FoldChange) > l2fc_threshold, na.rm = TRUE)
    } else NA

    # --- fGSEA results (from finals -- count significant from full table) ---
    fgsea_base <- safe_read(file.path(finals_dir, "fgsea_base.tsv"))
    fgsea_ruv  <- safe_read(file.path(finals_dir, "fgsea_ruv.tsv"))

    n_pathways_base <- if (!is.null(fgsea_base)) {
        sum(!is.na(fgsea_base$padj) & fgsea_base$padj < fgsea_fdr, na.rm = TRUE)
    } else NA
    n_pathways_ruv <- if (!is.null(fgsea_ruv)) {
        sum(!is.na(fgsea_ruv$padj) & fgsea_ruv$padj < fgsea_fdr, na.rm = TRUE)
    } else NA

    # --- Read preflight report (from intermediates) ---
    preflight      <- read_preflight_report(inter_dir)
    pf_status      <- if (!is.null(preflight) && "status" %in% colnames(preflight)) preflight$status[1] else NA
    pf_dropped     <- if (!is.null(preflight) && "dropped_covariates" %in% colnames(preflight)) preflight$dropped_covariates[1] else NA
    pf_rank_ok     <- if (!is.null(preflight) && "rank_ok" %in% colnames(preflight)) preflight$rank_ok[1] else NA
    pf_near_sing   <- if (!is.null(preflight) && "near_singular_vars" %in% colnames(preflight)) preflight$near_singular_vars[1] else NA
    pf_warnings    <- if (!is.null(preflight) && "warnings" %in% colnames(preflight)) preflight$warnings[1] else NA
    pf_errors      <- if (!is.null(preflight) && "errors" %in% colnames(preflight)) preflight$errors[1] else NA
    pf_n_genes     <- if (!is.null(preflight) && "n_genes" %in% colnames(preflight)) preflight$n_genes[1] else NA
    pf_n_control   <- if (!is.null(preflight) && "n_control" %in% colnames(preflight)) preflight$n_control[1] else NA
    pf_n_experimental <- if (!is.null(preflight) && "n_experimental" %in% colnames(preflight)) preflight$n_experimental[1] else NA

    # --- Helper: extract message from a log line ---
    # Log format: [timestamp] [LEVEL] contrast=X | biosample=Y | step=Z | msg=...
    # or:         [timestamp] [SKIP ] ... | reason=R | message text
    extract_log_message <- function(line) {
        # Try msg= first (ERROR lines)
        m <- regmatches(line, regexpr("msg=.+$", line))
        if (length(m) > 0) return(sub("^msg=", "", m))
        # Try reason= (SKIP lines): extract "reason | rest"
        m <- regmatches(line, regexpr("reason=.+$", line))
        if (length(m) > 0) return(sub("^reason=", "", m))
        # Fallback: everything after step=
        m <- regmatches(line, regexpr("step=.+$", line))
        if (length(m) > 0) return(m)
        return(line)
    }

    # --- Read stratum error log (used for all non-success statuses) ---
    err_lines_all  <- character(0)
    skip_lines_all <- character(0)
    if (!is.null(log_dir)) {
        stratum_safe_log <- sanitize_names(stratum)
        err_log_path <- file.path(log_dir, paste0(stratum_safe_log, ".errors.log"))
        if (file.exists(err_log_path)) {
            all_lines      <- readLines(err_log_path, warn = FALSE)
            err_lines_all  <- all_lines[grepl("\\[ERROR\\]", all_lines)]
            skip_lines_all <- all_lines[grepl("\\[SKIP \\]", all_lines)]
        }
    }

    # --- Determine status (fine-grained) ---
    #
    # Statuses (in priority order):
    #   success            — full pipeline: base DESeq2 + RUV DESeq2 both produced results
    #   success_base_only  — base DESeq2 produced results; RUV was intentionally skipped
    #                        (e.g. paired analysis uses donor as blocking factor) or RUV
    #                        failed. error_message explains why RUV was skipped.
    #   error              — DESeq2 crashed unexpectedly. error_message captures the crash
    #                        reason so it can be converted to a preflight check in future
    #                        updates (goal: every crash becomes a predictable skip).
    #   skipped_no_samples — 0 samples after stratum/min_cells filtering, or <3 total
    #   skipped_min_group  — samples exist but <3 per group, or n <= n_coefficients
    #   skipped_no_pairs   — paired analysis: no donors found in both groups
    #   skipped_preflight  — preflight validation failed: design matrix is rank-deficient,
    #                        formula can't be built, covariates are collinear, or count
    #                        matrix has structural issues (see error_message for details)
    #   skipped            — skipped for a reason not in the above categories;
    #                        error_message contains the skip reason from the log
    #   not_run            — stratum directory exists but no outputs and no log entries;
    #                        likely the pipeline was interrupted before this stratum ran
    #
    status    <- "unknown"
    error_msg <- NA

    # Collect ALL error messages for this stratum (for comprehensive reporting)
    all_error_msgs <- if (length(err_lines_all) > 0) {
        sapply(err_lines_all, extract_log_message, USE.NAMES = FALSE)
    } else character(0)

    if (!is.na(pf_status) && pf_status == "fail") {
        # ---- Preflight validation failed ----
        # Preflight checks the design matrix (rank, singularity), count matrix
        # structure, group membership, and sample counts BEFORE running DESeq2.
        status <- "skipped_preflight"
        if (!is.na(pf_errors) && nchar(pf_errors) > 0) {
            error_msg <- pf_errors
        } else if (length(all_error_msgs) > 0) {
            error_msg <- paste(unique(all_error_msgs), collapse = "; ")
        }

    } else if (!is.null(base_results)) {
        # ---- Base DESeq2 completed ----
        if (!is.null(ruv_results)) {
            status <- "success"
        } else {
            status <- "success_base_only"
            # Explain why RUV is missing
            ruv_skip <- skip_lines_all[grepl("paired_skip|ruvseq_max_k", skip_lines_all)]
            if (length(ruv_skip) > 0) {
                error_msg <- extract_log_message(ruv_skip[1])
            } else if (length(skip_lines_all) > 0) {
                ruv_skip2 <- skip_lines_all[grepl("04_run_ruvseq|05_run_deseq_final", skip_lines_all)]
                if (length(ruv_skip2) > 0) {
                    error_msg <- extract_log_message(ruv_skip2[1])
                }
            }
        }

    } else if (!is.null(coldata)) {
        # ---- Coldata produced but DESeq2 crashed ----
        # Report the crash so we can add it as a preflight check in the future.
        status <- "error"
        if (length(all_error_msgs) > 0) {
            # Exclude cascading "a dimension is zero" from RUV — that's always
            # secondary to the real crash. Report the original error.
            primary <- all_error_msgs[!grepl("a dimension is zero", all_error_msgs)]
            if (length(primary) > 0) {
                error_msg <- paste(unique(primary), collapse = "; ")
            } else {
                error_msg <- paste(unique(all_error_msgs), collapse = "; ")
            }
        } else {
            error_msg <- "DESeq2 crashed — see full log for details"
        }

    } else {
        # ---- No coldata: samples were eliminated before DESeq2 ----
        # Parse the log to determine WHY samples were eliminated.

        # First, gather all meaningful skip reasons (ignore upstream_skip propagations)
        real_skips <- skip_lines_all[!grepl("reason=upstream_skip", skip_lines_all)]
        first_skip <- if (length(real_skips) > 0) real_skips[1] else {
            if (length(skip_lines_all) > 0) skip_lines_all[1] else NULL
        }

        if (!is.null(first_skip)) {
            if (grepl("reason=filter_stratum", first_skip) ||
                grepl("only \\d+ samples? remain", first_skip)) {
                # Cell type has 0-2 samples total after min_cells / stratum filter
                status    <- "skipped_no_samples"
                error_msg <- extract_log_message(first_skip)

            } else if (grepl("reason=insufficient_samples", first_skip)) {
                # Samples exist but not enough for the model
                if (grepl("only \\d+ samples? remain", first_skip)) {
                    status <- "skipped_no_samples"
                } else {
                    # <3 per group, n <= n_coefficients, rank-deficient
                    status <- "skipped_min_group"
                }
                error_msg <- extract_log_message(first_skip)
                # Append the specific preflight errors if available
                if (length(all_error_msgs) > 0) {
                    error_msg <- paste(c(error_msg, unique(all_error_msgs)), collapse = "; ")
                }

            } else if (grepl("reason=paired_filter|no donors found in both groups", first_skip) ||
                       length(grep("paired_filter", err_lines_all)) > 0) {
                status    <- "skipped_no_pairs"
                error_msg <- if (length(all_error_msgs) > 0) {
                    paste(unique(all_error_msgs), collapse = "; ")
                } else {
                    extract_log_message(first_skip)
                }

            } else if (grepl("reason=preflight_counts|reason=column_check|reason=preflight_design|reason=preflight_fail",
                             first_skip)) {
                status    <- "skipped_preflight"
                error_msg <- if (length(all_error_msgs) > 0) {
                    paste(unique(all_error_msgs), collapse = "; ")
                } else {
                    extract_log_message(first_skip)
                }

            } else {
                # A skip reason we didn't anticipate — log it transparently
                status    <- "skipped"
                error_msg <- extract_log_message(first_skip)
            }
        } else if (length(all_error_msgs) > 0) {
            # No SKIP lines but there are ERROR lines — something crashed early
            status    <- "error"
            error_msg <- paste(unique(all_error_msgs), collapse = "; ")
        } else {
            # No outputs, no logs — pipeline never reached this stratum
            status    <- "not_run"
        }
    }

    # Build summary row
    row <- make_summary_row(
        biosample            = stratum,
        contrast_id          = contrast_id,
        status               = status,
        error_message        = error_msg,
        n_control_samps      = pf_n_control,
        n_experimental_samps = pf_n_experimental,
        n_samples_total      = n_total,
        base_formula         = base_formula,
        n_degs_base          = n_degs_base,
        ruv_formula          = ruv_formula,
        n_degs_ruv           = n_degs_ruv,
        n_pathways_base      = n_pathways_base,
        n_pathways_ruv       = n_pathways_ruv
    )

    # Append preflight columns
    row$preflight_status     <- pf_status
    row$n_genes_preflight    <- pf_n_genes
    row$dropped_covariates   <- pf_dropped
    row$rank_ok              <- pf_rank_ok
    row$near_singular_vars   <- pf_near_sing
    row$preflight_warnings   <- pf_warnings

    row
}


#' Collect error messages from a stratum's error log file (if it exists)
get_error_log_path <- function(log_dir, stratum) {
    stratum_safe <- sanitize_names(stratum)
    path <- file.path(log_dir, paste0(stratum_safe, ".errors.log"))
    if (file.exists(path) && file.size(path) > 0) path else NULL
}


#' Collect and concatenate a result file across all strata, adding metadata columns
collect_results <- function(strata, results_dir, finals_file, inter_file) {
    rows <- lapply(strata, function(s) {
        finals_dir <- file.path(results_dir, s, "finals")
        inter_dir  <- file.path(results_dir, s, "intermediates")
        df <- safe_read(file.path(finals_dir, finals_file))
        if (is.null(df)) return(NULL)
        mi <- safe_read(file.path(inter_dir, inter_file))
        df$stratum <- s
        df$formula <- if (!is.null(mi) && "formula" %in% colnames(mi)) mi$formula[1] else NA
        df
    })
    do.call(rbind, Filter(Negate(is.null), rows))
}


#' Collect and concatenate fGSEA results across all strata
collect_fgsea <- function(strata, results_dir, file_name) {
    rows <- lapply(strata, function(s) {
        df <- safe_read(file.path(results_dir, s, "finals", file_name))
        if (is.null(df)) return(NULL)
        df$stratum <- s
        df
    })
    do.call(rbind, Filter(Negate(is.null), rows))
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df

    cat(sprintf(
        "\n[%s] [INFO ] contrast=%s | step=08_aggregate_results | starting\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), contrast_id
    ))

    results_dir  <- file.path(cfg$outputs$results_dir, contrast_id)
    log_dir      <- file.path(cfg$outputs$logs_dir, contrast_id)
    summary_dir  <- file.path(results_dir, "summary")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Discover all strata that were run ---
    all_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
    strata   <- all_dirs[all_dirs != "summary"]

    # --- Discover expected strata from contrasts.csv ---
    contrast_row <- contrasts_df[contrasts_df$contrast_id == contrast_id, ]
    if (nrow(contrast_row) > 0) {
        expected_strata_safe <- sanitize_names(parse_pipe_field(contrast_row$strata[1]))
    } else {
        expected_strata_safe <- character(0)
    }
    missing_strata <- setdiff(expected_strata_safe, strata)

    if (length(strata) == 0 && length(missing_strata) == 0) {
        cat(sprintf("[WARN] No stratum directories found under: %s\n", results_dir))
        return(invisible(NULL))
    }

    cat(sprintf("[INFO] Found %d strata with results: %s\n",
                length(strata), paste(strata, collapse = ", ")))
    if (length(missing_strata) > 0) {
        cat(sprintf("[INFO] %d expected strata not run: %s\n",
                    length(missing_strata), paste(missing_strata, collapse = ", ")))
    }

    # --- Collect per-stratum summary rows ---
    fgsea_fdr <- cfg$fgsea$fdr_threshold %||% 0.10

    summary_rows <- lapply(strata, function(stratum) {
        inter_dir  <- file.path(results_dir, stratum, "intermediates")
        finals_dir <- file.path(results_dir, stratum, "finals")

        collect_stratum_summary(
            inter_dir      = inter_dir,
            finals_dir     = finals_dir,
            contrast_id    = contrast_id,
            stratum        = stratum,
            fdr_threshold  = cfg$deseq2$fdr_threshold,
            l2fc_threshold = cfg$deseq2$l2fc_threshold,
            fgsea_fdr      = fgsea_fdr,
            log_dir        = log_dir
        )
    })

    if (length(summary_rows) > 0) {
        contrast_summary <- do.call(rbind, summary_rows)
    } else {
        contrast_summary <- NULL
    }

    # --- Add rows for strata that were expected but never ran ---
    if (length(missing_strata) > 0) {
        missing_rows <- lapply(missing_strata, function(s) {
            row <- make_summary_row(
                biosample   = s,
                contrast_id = contrast_id,
                status      = "not_run"
            )
            row$preflight_status   <- NA
            row$n_genes_preflight  <- NA
            row$dropped_covariates <- NA
            row$rank_ok            <- NA
            row$near_singular_vars <- NA
            row$preflight_warnings <- NA
            row
        })
        missing_df <- do.call(rbind, missing_rows)
        contrast_summary <- if (!is.null(contrast_summary)) {
            rbind(contrast_summary, missing_df)
        } else {
            missing_df
        }
    }

    # =====================================================================
    # Write summary outputs
    # =====================================================================

    # --- 1. contrast_summary.csv ---
    summary_path <- file.path(summary_dir, "contrast_summary.csv")
    write_result(contrast_summary, summary_path, sep = ",")

    # --- 2. preflight_summary.csv ---
    preflight_rows <- lapply(strata, function(s) {
        read_preflight_report(file.path(results_dir, s, "intermediates"))
    })
    preflight_all <- do.call(rbind, Filter(Negate(is.null), preflight_rows))
    if (!is.null(preflight_all) && nrow(preflight_all) > 0) {
        write_result(preflight_all, file.path(summary_dir, "preflight_summary.csv"), sep = ",")
    }

    # --- 3. Consolidated DESeq2 tables ---
    deseq_base <- collect_results(strata, results_dir, "results_base.tsv", "model_info_base.csv")
    if (!is.null(deseq_base) && nrow(deseq_base) > 0) {
        deseq_base$model <- "Base"
        write_result(deseq_base, file.path(summary_dir, paste0(contrast_id, "_deseq_base.tsv")))
    }

    deseq_ruv <- collect_results(strata, results_dir, "results_ruv.tsv", "model_info_ruv.csv")
    if (!is.null(deseq_ruv) && nrow(deseq_ruv) > 0) {
        deseq_ruv$model <- "RUVseq"
        write_result(deseq_ruv, file.path(summary_dir, paste0(contrast_id, "_deseq_ruv.tsv")))
    }

    # --- 3b. Final best-model table: RUVseq where available, Base fallback ---
    # For each stratum, pick the best available model:
    #   - RUVseq if it ran (status = "success")
    #   - Base if only base ran (status = "success_base_only")
    # This gives one definitive results file per contrast.
    ruv_strata  <- if (!is.null(deseq_ruv))  unique(deseq_ruv$stratum)  else character(0)
    base_strata <- if (!is.null(deseq_base)) unique(deseq_base$stratum) else character(0)
    base_only   <- setdiff(base_strata, ruv_strata)

    deseq_final_parts <- list()
    if (!is.null(deseq_ruv) && nrow(deseq_ruv) > 0) {
        deseq_final_parts[["ruv"]] <- deseq_ruv
    }
    if (length(base_only) > 0 && !is.null(deseq_base)) {
        base_fallback <- deseq_base[deseq_base$stratum %in% base_only, , drop = FALSE]
        if (nrow(base_fallback) > 0) {
            deseq_final_parts[["base_fb"]] <- base_fallback
        }
    }
    if (length(deseq_final_parts) > 0) {
        deseq_final <- do.call(rbind, deseq_final_parts)
        rownames(deseq_final) <- NULL
        write_result(deseq_final,
                     file.path(summary_dir, paste0(contrast_id, "_deseq_final.tsv")))
        cat(sprintf("  Final results: %d strata (%d RUVseq, %d Base fallback)\n",
                    length(unique(deseq_final$stratum)),
                    length(ruv_strata),
                    length(base_only)))
    }

    # --- 4. Consolidated fGSEA tables ---
    fgsea_base <- collect_fgsea(strata, results_dir, "fgsea_base.tsv")
    if (!is.null(fgsea_base) && nrow(fgsea_base) > 0) {
        write_result(fgsea_base, file.path(summary_dir, paste0(contrast_id, "_fgsea_base.tsv")))
    }

    fgsea_ruv <- collect_fgsea(strata, results_dir, "fgsea_ruv.tsv")
    if (!is.null(fgsea_ruv) && nrow(fgsea_ruv) > 0) {
        write_result(fgsea_ruv, file.path(summary_dir, paste0(contrast_id, "_fgsea_ruv.tsv")))
    }

    # --- 5. Consolidated benchmark signatures ---
    bench_rows <- lapply(strata, function(s) {
        bench_path <- file.path(results_dir, s, "benchmarking", "benchmark_signatures.tsv")
        safe_read(bench_path)
    })
    bench_all <- do.call(rbind, Filter(Negate(is.null), bench_rows))
    if (!is.null(bench_all) && nrow(bench_all) > 0) {
        write_result(bench_all,
                     file.path(summary_dir, paste0(contrast_id, "_benchmark_signatures.tsv")))
        cat(sprintf("  Benchmark signatures: %d rows across %d strata\n",
                    nrow(bench_all), length(unique(bench_all$stratum))))
    }

    # =====================================================================
    # Console summary
    # =====================================================================

    cat("\n-- Run Summary -------------------------------------------------------\n")
    cat(sprintf("  Contrast: %s\n", contrast_id))
    cat(sprintf("  Strata:   %d total\n\n", nrow(contrast_summary)))

    status_counts <- table(contrast_summary$status)
    for (s in names(status_counts)) {
        cat(sprintf("    %-15s %d\n", s, status_counts[[s]]))
    }

    if (length(missing_strata) > 0) {
        cat(sprintf("\n  Not run:  %d strata (%s)\n",
                    length(missing_strata), paste(missing_strata, collapse = ", ")))
    }

    # --- Preflight summary ---
    if ("preflight_status" %in% colnames(contrast_summary)) {
        cat("\n  Preflight:\n")
        pf_vals <- contrast_summary$preflight_status
        pf_vals[is.na(pf_vals)] <- "not available"
        pf_counts <- table(pf_vals)
        for (s in names(pf_counts)) {
            cat(sprintf("    %-15s %d\n", s, pf_counts[[s]]))
        }
    }

    cat("\n  Per-stratum DEG counts:\n")
    for (i in seq_len(nrow(contrast_summary))) {
        row <- contrast_summary[i, ]
        cat(sprintf("    %-20s  status=%-15s  DEGs_base=%-6s  DEGs_ruv=%-6s  pathways_base=%-6s  pathways_ruv=%s\n",
                    row$biosample,
                    row$status,
                    ifelse(is.na(row$n_degs_base), "--", row$n_degs_base),
                    ifelse(is.na(row$n_degs_ruv),  "--", row$n_degs_ruv),
                    ifelse(is.na(row$n_pathways_base), "--", row$n_pathways_base),
                    ifelse(is.na(row$n_pathways_ruv),  "--", row$n_pathways_ruv)))
    }

    # =====================================================================
    # Log aggregation
    # =====================================================================

    cat("\n-- Log Summary -------------------------------------------------------\n")
    log_summary_path <- file.path(summary_dir, "log_summary.tsv")

    strata_with_errors <- Filter(Negate(is.null),
        lapply(strata, get_error_log_path, log_dir = log_dir))

    if (length(strata_with_errors) > 0) {
        cat(sprintf("  %d stratum/strata had warnings or errors.\n", length(strata_with_errors)))

        summarize_logs(log_dir, log_summary_path)
        cat(sprintf("  Log summary: %s\n", log_summary_path))

        # Split into errors-only and warnings-only
        error_log_files <- list.files(log_dir, pattern = "\\.errors\\.log$",
                                      full.names = TRUE, recursive = FALSE)
        if (length(error_log_files) > 0) {
            all_lines   <- unlist(lapply(error_log_files, readLines, warn = FALSE))
            error_lines <- all_lines[grepl("\\[(ERROR|SKIP )\\]", all_lines)]
            warn_lines  <- all_lines[grepl("\\[WARN \\]", all_lines)]

            if (length(error_lines) > 0) {
                writeLines(error_lines, file.path(summary_dir, "errors.log"))
                cat(sprintf("  Errors:   %s\n", file.path(summary_dir, "errors.log")))
            }
            if (length(warn_lines) > 0) {
                writeLines(warn_lines, file.path(summary_dir, "warnings.log"))
                cat(sprintf("  Warnings: %s\n", file.path(summary_dir, "warnings.log")))
            }
        }
    } else {
        cat("  No warnings or errors found across all strata.\n")
    }

    # =====================================================================
    # Status overview heatmap
    # =====================================================================

    if (!is.null(contrast_summary) && nrow(contrast_summary) > 0) {
        tryCatch({
            source(file.path(.script_dir, "utils", "plot_utils.R"))
            status_plot <- plot_status_heatmap(contrast_summary, contrast_id)
            if (!is.null(status_plot)) {
                # Adaptive width: ~0.6 inches per cell type, minimum 8
                n_types <- length(unique(contrast_summary$biosample))
                pdf_w   <- max(8, n_types * 0.6 + 2)
                pdf(file.path(summary_dir, "status_overview.pdf"),
                    width = pdf_w, height = 4)
                print(status_plot)
                dev.off()
                cat(sprintf("  Status heatmap: %s\n",
                            file.path(summary_dir, "status_overview.pdf")))
            }
        }, error = function(e) {
            cat(sprintf("  [WARN] Could not generate status heatmap: %s\n", e$message))
        })
    }

    # =====================================================================
    # Benchmark signature plots
    # =====================================================================

    if (!is.null(bench_all) && nrow(bench_all) > 0) {
        tryCatch({
            # plot_utils.R already sourced above for status heatmap

            # --- FORA-based specificity heatmaps (summary: strata on x-axis) ---
            bench_fora <- bench_all[bench_all$test == "fora", ]
            n_strata <- length(unique(bench_all$stratum))
            has_control_type <- "control_type" %in% colnames(bench_fora)

            if (nrow(bench_fora) > 0) {

                # Collapsed (class-level) data
                bench_collapsed <- collapse_benchmark_by_class(bench_fora)

                # Summary plots: 3 per model (base, ruv) x (positive, negative, collapsed)
                for (mc in list(list(gs = "base", label = "Base"),
                                list(gs = "ruv",  label = "RUV"))) {

                    save_heatmap(plot_fora_heatmap(
                        bench_fora,
                        control_type_val = if (has_control_type) "positive" else NULL,
                        categories_filter = if (!has_control_type) "cell_identity" else NULL,
                        title = sprintf("%s -- %s -- Positive Controls", contrast_id, mc$label),
                        fill_color = "#2E8B57", fill_limit = 0.6,
                        x_var = "stratum", gene_set_filter = mc$gs
                    ), file.path(summary_dir, sprintf("benchmark_%s_positive.pdf", mc$gs)))

                    save_heatmap(plot_fora_heatmap(
                        bench_fora,
                        control_type_val = if (has_control_type) "negative" else NULL,
                        categories_filter = if (!has_control_type) c("artifact", "subject_confounder") else NULL,
                        title = sprintf("%s -- %s -- Negative Controls", contrast_id, mc$label),
                        fill_color = "#D94040", fill_limit = 0.5,
                        x_var = "stratum", gene_set_filter = mc$gs
                    ), file.path(summary_dir, sprintf("benchmark_%s_negative.pdf", mc$gs)))

                    if (!is.null(bench_collapsed) && nrow(bench_collapsed) > 0) {
                        save_heatmap(plot_fora_heatmap(
                            bench_collapsed,
                            title = sprintf("%s -- %s -- Class-Level Collapsed", contrast_id, mc$label),
                            fill_color = "#7B3294", fill_limit = 0.3,
                            collapsed = TRUE,
                            x_var = "stratum", gene_set_filter = mc$gs
                        ), file.path(summary_dir, sprintf("benchmark_%s_collapsed.pdf", mc$gs)))
                    }
                }
            }

        }, error = function(e) {
            cat(sprintf("  [WARN] Could not generate benchmark plots: %s\n", e$message))
        })
    }

    cat("\n-- Complete ----------------------------------------------------------\n")
    cat(sprintf("  Results: %s\n\n", summary_dir))

    invisible(contrast_summary)
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {

    opt_list <- list(
        make_option("--config",   type = "character", default = "config/config.yaml"),
        make_option("--contrast", type = "character",
                    help = "contrast_id to aggregate (e.g. T1D_vs_ND)")
    )
    opts <- parse_args(OptionParser(option_list = opt_list))
    if (is.null(opts$contrast)) stop("--contrast is required")

    main(opts$config, opts$contrast)
}
