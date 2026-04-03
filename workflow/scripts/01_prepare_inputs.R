# =============================================================================
# 01_prepare_inputs.R
# EasyDE — Input Validation
#
# PURPOSE:
#   Validates all inputs before any compute runs. Reads config.yaml and
#   contrasts.csv, checks that every required file exists, every required
#   column is present, every threshold is sane, and every cross-reference
#   between files is consistent.
#
#   Collects ALL issues before reporting — you see every problem at once.
#   Structural errors (broken CSV, duplicates) hard-stop the pipeline.
#   Per-contrast errors (missing columns, bad groups) are flagged and
#   logged — the affected contrasts will be skipped in step 02.
#
# USAGE:
#   Rscript workflow/scripts/01_prepare_inputs.R --config config/config.yaml
#
# OUTPUT:
#   Prints a structured validation report to stdout.
#   Writes logs/validation.log
#   On success: exits 0 and writes results/ + logs/ directory stubs.
#   On failure: exits 1 with a consolidated error list.
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
source(file.path(.script_dir, "utils", "validation_utils.R"))


# =============================================================================
# Validation helpers
# =============================================================================

#' Check that a required key exists in a nested config list
#' Returns an error string or NULL
check_config_key <- function(cfg, keys, label) {
    val <- tryCatch(
        Reduce(`[[`, keys, cfg),
        error = function(e) NULL
    )
    if (is.null(val)) {
        return(sprintf("config.yaml missing required key: %s", label))
    }
    NULL
}


# =============================================================================
# Section validators
# =============================================================================

#' Validate config.yaml structure, paths, and value ranges
#'
#' @param cfg   Parsed config list from yaml::read_yaml()
#' @return Named list: $errors (character), $warnings (character)
validate_config <- function(cfg) {

    errors   <- character(0)
    warnings <- character(0)

    cat("\n── Validating config.yaml ──────────────────────────────────────\n")

    # --- Required top-level sections ---
    required_sections <- c("inputs", "outputs", "steps", "write_outputs",
                           "filtering", "deseq2", "ruvseq", "fgsea", "logging")
    for (section in required_sections) {
        if (is.null(cfg[[section]])) {
            errors <- add_error(errors, sprintf("Missing required section: '%s'", section))
        }
    }

    # --- Input file paths ---
    if (!is.null(cfg$inputs)) {

        if (!dir.exists(cfg$inputs$counts_dir)) {
            errors <- add_error(errors,
                sprintf("counts_dir not found: '%s'", cfg$inputs$counts_dir))
        } else {
            cat(sprintf("  [OK] counts_dir: %s\n", cfg$inputs$counts_dir))
        }

        if (!file.exists(cfg$inputs$sample_metadata)) {
            errors <- add_error(errors,
                sprintf("sample_metadata not found: '%s'", cfg$inputs$sample_metadata))
        } else {
            cat(sprintf("  [OK] sample_metadata: %s\n", cfg$inputs$sample_metadata))
        }

        if (!file.exists(cfg$inputs$contrasts_file)) {
            errors <- add_error(errors,
                sprintf("contrasts_file not found: '%s'", cfg$inputs$contrasts_file))
        } else {
            cat(sprintf("  [OK] contrasts_file: %s\n", cfg$inputs$contrasts_file))
        }
    }

    # --- steps flags must be logical ---
    if (!is.null(cfg$steps)) {
        for (flag in names(cfg$steps)) {
            if (!is.logical(cfg$steps[[flag]])) {
                errors <- add_error(errors,
                    sprintf("steps.%s must be true or false, got: '%s'", flag, cfg$steps[[flag]]))
            }
        }

        # Logical dependency: fgsea.run_on = "ruvseq" requires steps.ruvseq = true
        if (isFALSE(cfg$steps$ruvseq) && !is.null(cfg$fgsea$run_on) &&
            cfg$fgsea$run_on == "ruvseq") {
            errors <- add_error(errors,
                "fgsea.run_on = 'ruvseq' but steps.ruvseq = false. ",
                "Set fgsea.run_on to 'base' or enable steps.ruvseq.")
        }

        # Warn if ruvseq off but run_on = "both"
        if (isFALSE(cfg$steps$ruvseq) && !is.null(cfg$fgsea$run_on) &&
            cfg$fgsea$run_on == "both") {
            warnings <- add_warning(warnings,
                "fgsea.run_on = 'both' but steps.ruvseq = false. fGSEA will run on base results only.")
        }
    }

    # --- DESeq2 thresholds ---
    if (!is.null(cfg$deseq2)) {
        fdr <- cfg$deseq2$fdr_threshold
        if (!is.null(fdr) && (!is.numeric(fdr) || fdr <= 0 || fdr >= 1)) {
            errors <- add_error(errors,
                sprintf("deseq2.fdr_threshold must be between 0 and 1, got: %s", fdr))
        }

        l2fc <- cfg$deseq2$l2fc_threshold
        if (!is.null(l2fc) && (!is.numeric(l2fc) || l2fc < 0)) {
            errors <- add_error(errors,
                sprintf("deseq2.l2fc_threshold must be >= 0, got: %s", l2fc))
        }
    }

    # --- RUVseq thresholds ---
    if (!is.null(cfg$ruvseq)) {
        k <- cfg$ruvseq$max_latent_factors
        if (!is.null(k) && (!is.numeric(k) || k < 1)) {
            errors <- add_error(errors,
                sprintf("ruvseq.max_latent_factors must be >= 1, got: %s", k))
        }

        pval <- cfg$ruvseq$empirical_pval_threshold
        if (!is.null(pval) && (!is.numeric(pval) || pval <= 0 || pval >= 1)) {
            errors <- add_error(errors,
                sprintf("ruvseq.empirical_pval_threshold must be between 0 and 1, got: %s", pval))
        }
    }

    # --- fGSEA settings ---
    if (!is.null(cfg$fgsea)) {

        valid_run_on <- c("base", "ruvseq", "both")
        if (!is.null(cfg$fgsea$run_on) && !cfg$fgsea$run_on %in% valid_run_on) {
            errors <- add_error(errors,
                sprintf("fgsea.run_on must be one of: %s. Got: '%s'",
                        paste(valid_run_on, collapse = ", "), cfg$fgsea$run_on))
        }

        if (!is.null(cfg$fgsea$gene_sets_file) && isTRUE(cfg$steps$fgsea)) {
            if (!file.exists(cfg$fgsea$gene_sets_file)) {
                errors <- add_error(errors,
                    sprintf("fgsea.gene_sets_file not found: '%s'", cfg$fgsea$gene_sets_file))
            } else {
                cat(sprintf("  [OK] gene_sets_file: %s\n", cfg$fgsea$gene_sets_file))
            }

            for (gene_list in cfg$fgsea$exclude_gene_lists) {
                if (!file.exists(gene_list)) {
                    errors <- add_error(errors,
                        sprintf("fgsea.exclude_gene_lists file not found: '%s'", gene_list))
                } else {
                    cat(sprintf("  [OK] exclude_gene_list: %s\n", gene_list))
                }
            }
        }

        min_s <- cfg$fgsea$min_gene_set_size
        max_s <- cfg$fgsea$max_gene_set_size
        if (!is.null(min_s) && !is.null(max_s) && min_s >= max_s) {
            errors <- add_error(errors,
                sprintf("fgsea.min_gene_set_size (%d) must be less than max_gene_set_size (%d)",
                        min_s, max_s))
        }
    }

    # --- fGSEA ranking stat ---
    if (!is.null(cfg$fgsea$ranking_stat)) {
        valid_stats <- c("stat", "lfc_pvalue", "lfc_shrunk")
        if (!cfg$fgsea$ranking_stat %in% valid_stats) {
            errors <- add_error(errors,
                sprintf("fgsea.ranking_stat must be one of: %s. Got: '%s'",
                        paste(valid_stats, collapse = ", "), cfg$fgsea$ranking_stat))
        }
        # Warn if using shrunk LFC but shrinkage is disabled
        if (cfg$fgsea$ranking_stat == "lfc_shrunk" && isFALSE(cfg$deseq2$lfc_shrinkage)) {
            warnings <- add_warning(warnings,
                "fgsea.ranking_stat = 'lfc_shrunk' but deseq2.lfc_shrinkage = false. Will fall back to unshrunk LFC at runtime.")
        }
    }

    # --- Donor column for paired analysis ---
    if (isTRUE(cfg$steps$paired_analysis)) {
        donor_col <- cfg$filtering$donor_col
        if (is.null(donor_col) || !is.character(donor_col) || trimws(donor_col) == "") {
            errors <- add_error(errors,
                "steps.paired_analysis = true but filtering.donor_col is missing or empty")
        }
    }

    # --- Logging level ---
    if (!is.null(cfg$logging$level)) {
        valid_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
        if (!toupper(cfg$logging$level) %in% valid_levels) {
            errors <- add_error(errors,
                sprintf("logging.level must be one of: %s. Got: '%s'",
                        paste(valid_levels, collapse = ", "), cfg$logging$level))
        }
    }

    list(errors = errors, warnings = warnings)
}


#' Validate contrasts.csv structure and cross-references
#'
#' Returns structural errors (hard-stop: broken CSV format, duplicates) separately
#' from per-contrast issues (soft: flagged and skipped downstream).
#'
#' @param contrasts_df   data.frame loaded from contrasts.csv
#' @param sample_meta    data.frame loaded from sample_metadata.csv
#' @param counts_dir     Path to counts directory (to check stratum files exist)
#' @param cfg            Parsed config list (optional, used for donor_col validation)
#' @return Named list:
#'   $structural_errors  — character vector (hard-stop)
#'   $warnings           — character vector (informational)
#'   $per_contrast_issues — named list keyed by contrast_id (soft-gate)
validate_contrasts <- function(contrasts_df, sample_meta, counts_dir, cfg = NULL) {

    structural_errors   <- character(0)
    warnings            <- character(0)
    per_contrast_issues <- list()

    cat("\n\u2500\u2500 Validating contrasts.csv \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

    # --- Required columns (structural — hard-stop) ---
    required_cols <- c("contrast_id", "contrast_var", "biosample_id_col",
                       "control_grp", "experimental_grp", "strata")
    missing_cols  <- setdiff(required_cols, colnames(contrasts_df))
    if (length(missing_cols) > 0) {
        structural_errors <- add_error(structural_errors,
            sprintf("contrasts.csv missing required columns: %s",
                    paste(missing_cols, collapse = ", ")))
        return(list(structural_errors = structural_errors,
                    warnings = warnings, per_contrast_issues = per_contrast_issues))
    }

    # --- No duplicate contrast_ids (structural — hard-stop) ---
    dupes <- contrasts_df$contrast_id[duplicated(contrasts_df$contrast_id)]
    if (length(dupes) > 0) {
        structural_errors <- add_error(structural_errors,
            sprintf("Duplicate contrast_id values: %s", paste(unique(dupes), collapse = ", ")))
    }

    # --- Row-level validation (per-contrast — soft-gate) ---
    for (i in seq_len(nrow(contrasts_df))) {

        row        <- contrasts_df[i, ]
        cid        <- row$contrast_id
        prefix     <- sprintf("contrast '%s'", cid)
        row_errors   <- character(0)
        row_warnings <- character(0)

        # control_grp and experimental_grp must be both filled or both empty
        has_ctrl <- !is.na(row$control_grp) && trimws(row$control_grp) != ""
        has_exp  <- !is.na(row$experimental_grp) && trimws(row$experimental_grp) != ""
        if (xor(has_ctrl, has_exp)) {
            row_errors <- add_error(row_errors,
                sprintf("%s: control_grp and experimental_grp must both be filled or both empty", prefix))
        }

        # biosample_id_col must exist in sample_metadata
        bid_col <- trimws(row$biosample_id_col)
        if (!bid_col %in% colnames(sample_meta)) {
            row_errors <- add_error(row_errors,
                sprintf("%s: biosample_id_col '%s' not found in sample_metadata columns",
                        prefix, bid_col))
        }

        # contrast_var must exist in sample_metadata
        cvar <- trimws(row$contrast_var)
        if (!cvar %in% colnames(sample_meta)) {
            row_errors <- add_error(row_errors,
                sprintf("%s: contrast_var '%s' not found in sample_metadata columns",
                        prefix, cvar))
        }

        # covariates must all exist in sample_metadata
        if (!is.na(row$covariates) && trimws(row$covariates) != "") {
            covs <- trimws(strsplit(row$covariates, "\\|")[[1]])
            missing_covs <- setdiff(covs, colnames(sample_meta))
            if (length(missing_covs) > 0) {
                row_errors <- add_error(row_errors,
                    sprintf("%s: covariate(s) not found in sample_metadata: %s",
                            prefix, paste(missing_covs, collapse = ", ")))
            }
        }

        # latent_corr_vars must all exist in sample_metadata (warn only — not fatal)
        if (!is.na(row$latent_corr_vars) && trimws(row$latent_corr_vars) != "") {
            corr_vars <- trimws(strsplit(row$latent_corr_vars, "\\|")[[1]])
            missing_corr <- setdiff(corr_vars, colnames(sample_meta))
            if (length(missing_corr) > 0) {
                row_warnings <- add_warning(row_warnings,
                    sprintf("%s: latent_corr_vars not found in sample_metadata (will be skipped): %s",
                            prefix, paste(missing_corr, collapse = ", ")))
            }
        }

        # --- PREFLIGHT CHECKS ---

        # Pipe-field validation: catch empty tokens (e.g. "Sex||BMI")
        if (!is.na(row$strata) && trimws(row$strata) != "") {
            pf_result <- validate_pipe_field(row$strata, "strata", context = prefix)
            row_warnings <- c(row_warnings, pf_result$warnings)
        }
        if (!is.na(row$covariates) && trimws(row$covariates) != "") {
            pf_result <- validate_pipe_field(row$covariates, "covariates", context = prefix)
            row_warnings <- c(row_warnings, pf_result$warnings)
        }
        if ("latent_corr_vars" %in% colnames(row) &&
            !is.na(row$latent_corr_vars) && trimws(row$latent_corr_vars) != "") {
            pf_result <- validate_pipe_field(row$latent_corr_vars, "latent_corr_vars", context = prefix)
            row_warnings <- c(row_warnings, pf_result$warnings)
        }

        # Sanitized name collision check
        all_model_vars <- c(cvar)
        if (!is.na(row$covariates) && trimws(row$covariates) != "") {
            all_model_vars <- c(all_model_vars, covs)
        }
        sn_result    <- validate_sanitized_names(all_model_vars, context = prefix)
        row_errors   <- c(row_errors, sn_result$errors)

        # Metadata column quality: check for all-NA columns
        if (cvar %in% colnames(sample_meta)) {
            mc_result    <- validate_metadata_column(sample_meta, cvar, context = prefix)
            row_errors   <- c(row_errors, mc_result$errors)
            row_warnings <- c(row_warnings, mc_result$warnings)
        }
        if (bid_col %in% colnames(sample_meta)) {
            mc_result    <- validate_metadata_column(sample_meta, bid_col, context = prefix)
            row_errors   <- c(row_errors, mc_result$errors)
            row_warnings <- c(row_warnings, mc_result$warnings)
        }
        if (!is.na(row$covariates) && trimws(row$covariates) != "") {
            for (cov in covs) {
                if (cov %in% colnames(sample_meta)) {
                    mc_result    <- validate_metadata_column(sample_meta, cov, context = prefix)
                    row_errors   <- c(row_errors, mc_result$errors)
                    row_warnings <- c(row_warnings, mc_result$warnings)
                }
            }
        }

        # Group membership: verify control/experimental values exist in data
        if (has_ctrl && has_exp && cvar %in% colnames(sample_meta)) {
            gm_result    <- validate_group_membership(
                sample_meta, cvar,
                trimws(row$control_grp), trimws(row$experimental_grp),
                context = prefix
            )
            row_errors   <- c(row_errors, gm_result$errors)
            row_warnings <- c(row_warnings, gm_result$warnings)
        }

        # Donor column validation
        if (!is.null(cfg) && !is.null(cfg$filtering$donor_col)) {
            dc_result    <- validate_donor_column(sample_meta, cfg$filtering$donor_col, context = prefix)
            row_errors   <- c(row_errors, dc_result$errors)
            row_warnings <- c(row_warnings, dc_result$warnings)
        }

        # filter_col_N / filter_val_N pairs must be matched (no orphaned columns)
        filter_col_names <- grep("^filter_col_", colnames(row), value = TRUE)
        for (col_name in filter_col_names) {
            n       <- sub("^filter_col_", "", col_name)
            val_col <- paste0("filter_val_", n)
            col_val <- row[[col_name]]
            val_val <- if (val_col %in% colnames(row)) row[[val_col]] else NA

            has_col <- !is.na(col_val) && trimws(col_val) != ""
            has_val <- !is.na(val_val) && trimws(val_val) != ""

            if (xor(has_col, has_val)) {
                row_errors <- add_error(row_errors,
                    sprintf("%s: filter pair %s / %s must both be filled or both empty",
                            prefix, col_name, val_col))
            }

            # filter column must exist in sample_metadata
            if (has_col && trimws(col_val) %in% colnames(sample_meta)) {
                cat(sprintf("  [OK] %s: filter_col_%s = '%s' found in metadata\n",
                            cid, n, trimws(col_val)))
            } else if (has_col) {
                row_errors <- add_error(row_errors,
                    sprintf("%s: filter_col_%s value '%s' not found in sample_metadata columns",
                            prefix, n, trimws(col_val)))
            }
        }

        # strata — each must resolve to a counts file
        strata_list <- trimws(strsplit(row$strata, "\\|")[[1]])
        if (length(strata_list) == 0 || all(strata_list == "")) {
            row_errors <- add_error(row_errors, sprintf("%s: strata field is empty", prefix))
        } else {
            for (stratum in strata_list) {
                safe_name    <- gsub("[^[:alnum:]]", "_", stratum)
                counts_file  <- file.path(counts_dir,
                                          paste0(safe_name, ".counts.csv"))
                # Also accept .csv.gz
                counts_file_gz <- paste0(counts_file, ".gz")

                if (!file.exists(counts_file) && !file.exists(counts_file_gz)) {
                    row_errors <- add_error(row_errors,
                        sprintf("%s: no count matrix found for stratum '%s'\n             Expected: %s[.gz]",
                                prefix, stratum, counts_file))
                } else {
                    # Deep structural validation of count matrix
                    actual_path <- if (file.exists(counts_file)) counts_file else counts_file_gz
                    cm_result   <- validate_count_matrix(actual_path, stratum_name = stratum)
                    row_errors  <- c(row_errors,   cm_result$errors)
                    row_warnings <- c(row_warnings, cm_result$warnings)
                    if (length(cm_result$errors) == 0) {
                        cat(sprintf("  [OK] %s: counts file valid for stratum '%s' (%d samples)\n",
                                    cid, stratum, cm_result$n_samples))
                    }
                }
            }
        }

        # Store per-contrast results
        if (length(row_errors) > 0) {
            per_contrast_issues[[cid]] <- row_errors
        }
        warnings <- c(warnings, row_warnings)
    }

    list(structural_errors   = structural_errors,
         warnings            = warnings,
         per_contrast_issues  = per_contrast_issues)
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path) {

    all_errors   <- character(0)
    all_warnings <- character(0)

    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║           EasyDE — Input Validation                         ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n")
    cat(sprintf("Config: %s\n", config_path))
    cat(sprintf("Time:   %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    # --- Load config ---
    cat("\n── Loading config.yaml ──────────────────────────────────────────\n")
    if (!file.exists(config_path)) {
        stop(sprintf("Cannot find config file: '%s'", config_path))
    }
    cfg <- suppressWarnings(yaml::read_yaml(config_path))
    cat("  [OK] config.yaml loaded\n")

    # --- Validate config ---
    config_result  <- validate_config(cfg)
    all_errors     <- c(all_errors,   config_result$errors)
    all_warnings   <- c(all_warnings, config_result$warnings)

    # Stop here if config is broken — can't validate contrasts without valid paths
    if (length(all_errors) > 0) {
        cat("\n── Validation failed — fix config errors before continuing ─────\n")
        cat(paste(all_errors, collapse = "\n"), "\n")
        stop("Validation failed. See errors above.", call. = FALSE)
    }

    # --- Load metadata and contrasts ---
    cat("\n── Loading input files ──────────────────────────────────────────\n")

    sample_meta <- tryCatch(
        fread(cfg$inputs$sample_metadata) %>% as.data.frame(),
        error = function(e) {
            stop(sprintf("Could not read sample_metadata: %s", e$message), call. = FALSE)
        }
    )
    cat(sprintf("  [OK] sample_metadata: %d rows x %d columns\n",
                nrow(sample_meta), ncol(sample_meta)))

    # Skip comment lines (lines starting with #) when reading contrasts
    raw_lines    <- suppressWarnings(readLines(cfg$inputs$contrasts_file))
    data_lines   <- raw_lines[!grepl("^#", raw_lines) & trimws(raw_lines) != ""]
    contrasts_df <- fread(text = paste(data_lines, collapse = "\n")) %>% as.data.frame()
    cat(sprintf("  [OK] contrasts.csv: %d contrasts defined\n", nrow(contrasts_df)))

    # --- Validate contrasts ---
    contrast_result <- validate_contrasts(contrasts_df, sample_meta, cfg$inputs$counts_dir, cfg = cfg)
    all_errors      <- c(all_errors,   contrast_result$structural_errors)
    all_warnings    <- c(all_warnings, contrast_result$warnings)
    flagged         <- contrast_result$per_contrast_issues   # named list by contrast_id

    # --- Final report ---
    cat("\n\u2500\u2500 Validation Report \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

    if (length(all_warnings) > 0) {
        cat(sprintf("\n%d warning(s):\n", length(all_warnings)))
        cat(paste(all_warnings, collapse = "\n"), "\n")
    }

    # Structural errors still hard-stop (broken CSV, duplicates)
    if (length(all_errors) > 0) {
        cat(sprintf("\n%d structural error(s) -- pipeline cannot run:\n", length(all_errors)))
        cat(paste(all_errors, collapse = "\n"), "\n")
        cat("\nFix all structural errors above and re-run validation.\n")
        stop("Validation failed.", call. = FALSE)
    }

    # Per-contrast issues: flagged, will be skipped downstream by step 02
    if (length(flagged) > 0) {
        cat(sprintf("\n%d contrast(s) flagged -- these will be skipped in downstream steps:\n",
                    length(flagged)))
        for (cid in names(flagged)) {
            cat(sprintf("\n  [FLAGGED] %s:\n", cid))
            for (msg in flagged[[cid]]) {
                cat(sprintf("    - %s\n", msg))
            }
        }
    }

    n_clean <- nrow(contrasts_df) - length(flagged)
    if (n_clean > 0) {
        cat(sprintf("\n  %d contrast(s) passed all checks.\n", n_clean))
    } else {
        cat("\n  WARNING: All contrasts are flagged. No clean contrasts to run.\n")
    }

    # --- Create output directory stubs ---
    dir.create(cfg$outputs$results_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(cfg$outputs$logs_dir,    recursive = TRUE, showWarnings = FALSE)

    # Write a validation summary for the audit trail
    validation_log <- file.path(cfg$outputs$logs_dir, "validation.log")
    writeLines(c(
        sprintf("EasyDE validation: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Config:    %s", config_path),
        sprintf("Contrasts: %d defined (%d clean, %d flagged)",
                nrow(contrasts_df), n_clean, length(flagged)),
        sprintf("Strata:    %s", paste(unique(unlist(
            strsplit(contrasts_df$strata, "\\|"))), collapse = ", ")),
        sprintf("Warnings:  %d", length(all_warnings)),
        if (length(flagged) > 0) sprintf("Flagged:   %s", paste(names(flagged), collapse = ", ")) else "",
        if (length(all_warnings) > 0) paste(all_warnings, collapse = "\n") else ""
    ), validation_log)

    cat(sprintf("\n  Validation log written: %s\n", validation_log))
    cat("  Ready to run.\n\n")

    invisible(list(cfg = cfg, contrasts = contrasts_df, sample_meta = sample_meta))
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive()) {

    opt_list <- list(
        make_option("--config", type = "character",
                    default = "config/config.yaml",
                    help = "Path to config.yaml [default: config/config.yaml]")
    )

    opts <- parse_args(OptionParser(option_list = opt_list))
    tryCatch({
        main(opts$config)
        quit(status = 0)
    }, error = function(e) {
        message("ERROR: ", conditionMessage(e))
        quit(status = 1)
    })
}