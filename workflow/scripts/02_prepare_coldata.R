# =============================================================================
# 02_prepare_coldata.R
# EasyDE — Sample Metadata Preparation (per stratum)
#
# PURPOSE:
#   For one contrast x stratum combination, loads the count matrix and
#   sample metadata, applies all filters, deduplicates donors, applies
#   paired analysis filtering if configured, and writes the final coldata
#   (the sample table DESeq2 will use) to disk.
#
#   This script is the "sample filtering" step — it decides which biosamples
#   survive into the DE analysis and why. Every filtering decision is logged.
#
# INPUTS:
#   - data/counts/{stratum_safe}.counts.csv     count matrix for this stratum
#   - data/sample_metadata.csv                  full sample metadata table
#   - config/config.yaml                        pipeline settings
#   - contrasts/contrasts.csv                   contrast definitions
#
# OUTPUTS:
#   - results/{contrast_id}/{stratum}/intermediates/coldata.csv
#   - results/{contrast_id}/{stratum}/intermediates/counts_filtered.csv
#   - logs/{contrast_id}/{stratum}.log
#
# USAGE:
#   Rscript workflow/scripts/02_prepare_coldata.R \
#       --config  config/config.yaml \
#       --contrast T1D_vs_ND \
#       --stratum  Beta
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(stringr)
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
source(file.path(.script_dir, "utils", "filter_utils.R"))
source(file.path(.script_dir, "utils", "logging_utils.R"))
source(file.path(.script_dir, "utils", "stats_utils.R"))
source(file.path(.script_dir, "utils", "validation_utils.R"))


# =============================================================================
# Script-specific helpers
# =============================================================================

#' Apply paired analysis filter: retain only donors present in BOTH groups
#'
#' When paired_analysis is enabled, we want matched donors — one biosample
#' in the control group and one in the experimental group for the same donor.
#' This removes donors who only appear in one arm of the comparison.
#'
#' @param meta             data.frame of sample metadata (already group-filtered)
#' @param contrast_var     Column name for the contrast variable
#' @param control_grp      Value defining the control group
#' @param experimental_grp Value defining the experimental group
#' @param donor_col        Column name for the donor/subject identifier
#' @param logger           Logger instance from make_logger()
#' @return Filtered data.frame containing only paired donors
apply_paired_filter <- function(meta,
                                contrast_var,
                                control_grp,
                                experimental_grp,
                                donor_col,
                                logger) {

    ctrl_donors <- meta[meta[[contrast_var]] == control_grp,     donor_col]
    exp_donors  <- meta[meta[[contrast_var]] == experimental_grp, donor_col]
    paired_donors <- intersect(ctrl_donors, exp_donors)

    before <- nrow(meta)
    meta   <- meta[meta[[donor_col]] %in% paired_donors, ]
    after  <- nrow(meta)

    logger$info("paired_filter",
        sprintf("donors in both groups: %d | samples: %d -> %d",
                length(paired_donors), before, after))

    if (length(paired_donors) == 0) {
        logger$error("paired_filter", "no donors found in both groups after pairing")
        return(NULL)
    }

    return(meta)
}


#' Prepare and sanitize coldata for DESeq2
#'
#' DESeq2 requires:
#'   - rownames matching colnames of the count matrix
#'   - no spaces or special characters in column names (R formula safety)
#'   - character columns that are contrast/covariate variables cast to factors
#'   - optionally: numeric covariates z-score scaled
#'
#' Returns the sanitized coldata AND the mapping from sanitized names back
#' to original names (needed for result reporting).
#'
#' @param meta             data.frame of filtered sample metadata
#' @param biosample_id_col Column name holding the biosample ID
#' @param contrast_var     Column name for the contrast variable
#' @param covariates       Character vector of covariate column names
#' @param scale_numeric    Logical: z-score scale numeric covariates?
#' @return Named list: $coldata (sanitized data.frame), $name_map (named character vector)
prepare_coldata <- function(meta,
                            biosample_id_col,
                            contrast_var,
                            covariates,
                            scale_numeric = TRUE) {

    rownames(meta) <- meta[[biosample_id_col]]

    # Build name mapping: original -> sanitized
    all_vars  <- c(contrast_var, covariates)
    safe_vars <- sanitize_names(all_vars)
    name_map  <- setNames(safe_vars, all_vars)

    # Sanitize all column names
    colnames(meta) <- sanitize_names(colnames(meta))

    safe_contrast  <- name_map[[contrast_var]]
    safe_covs      <- name_map[covariates]

    # Cast character and logical columns to factors.
    # DESeq2's results() requires the contrast variable to be a factor when using
    # contrast = c(var, exp, ctrl). Logical columns (e.g. AAB TRUE/FALSE) are not
    # automatically coerced by R model formulas, causing "X is not a factor" errors.
    for (v in c(safe_contrast, safe_covs)) {
        if (v %in% colnames(meta) && (is.character(meta[[v]]) || is.logical(meta[[v]]))) {
            meta[[v]] <- as.factor(meta[[v]])
        }
    }

    # Z-score scale numeric covariates (not the contrast variable)
    if (scale_numeric) {
        for (v in safe_covs) {
            if (v %in% colnames(meta) && is.numeric(meta[[v]])) {
                meta[[v]] <- as.numeric(scale(meta[[v]]))
            }
        }
    }

    # Drop covariates with only one unique value — they'll crash model fitting
    dropped <- character(0)
    for (v in safe_covs) {
        if (v %in% colnames(meta) && length(unique(na.omit(meta[[v]]))) < 2) {
            dropped <- c(dropped, v)
        }
    }

    return(list(
        coldata  = meta,
        name_map = name_map,
        dropped_covariates = dropped
    ))
}


#' Check whether enough samples remain for the contrast and model
#'
#' DESeq2 will fail silently or cryptically if there are fewer samples than
#' model coefficients. This check catches that before it happens.
#'
#' @param coldata          Sanitized coldata data.frame
#' @param contrast_var     Sanitized contrast variable name
#' @param control_grp      Control group value (or NULL for numeric)
#' @param experimental_grp Experimental group value (or NULL for numeric)
#' @param covariates       Character vector of sanitized covariate names (active ones only)
#' @param logger           Logger instance
#' @return Named list: $n_control, $n_experimental, $n_total, $sufficient (logical)
check_sample_counts <- function(coldata,
                                contrast_var,
                                control_grp,
                                experimental_grp,
                                covariates,
                                logger) {

    is_categorical <- is.factor(coldata[[contrast_var]]) ||
                      is.character(coldata[[contrast_var]])

    if (is_categorical) {
        n_ctrl <- sum(coldata[[contrast_var]] == control_grp,     na.rm = TRUE)
        n_exp  <- sum(coldata[[contrast_var]] == experimental_grp, na.rm = TRUE)
        n_tot  <- nrow(coldata)

        logger$info("sample_counts",
            sprintf("control=%s (n=%d) | experimental=%s (n=%d) | total=%d",
                    control_grp, n_ctrl, experimental_grp, n_exp, n_tot))

        # Need at least 3 per group and more samples than coefficients
        # Count actual expanded columns: factors with k levels contribute k-1 columns
        n_coeffs <- 1L + 1L  # intercept + contrast (binary)
        for (cv in covariates) {
            if (is.factor(coldata[[cv]]) || is.character(coldata[[cv]])) {
                n_coeffs <- n_coeffs + max(1L, nlevels(as.factor(coldata[[cv]])) - 1L)
            } else {
                n_coeffs <- n_coeffs + 1L
            }
        }
        sufficient <- (n_ctrl >= 3) && (n_exp >= 3) && (n_tot > n_coeffs)

        if (n_ctrl < 3) logger$error("sample_counts",
            sprintf("fewer than 3 control samples (n=%d)", n_ctrl))
        if (n_exp < 3)  logger$error("sample_counts",
            sprintf("fewer than 3 experimental samples (n=%d)", n_exp))
        if (n_tot <= n_coeffs) logger$error("sample_counts",
            sprintf("n_samples (%d) <= n_coefficients (%d)", n_tot, n_coeffs))

    } else {
        n_ctrl <- NA
        n_exp  <- NA
        n_tot  <- nrow(coldata)
        # Count actual expanded columns (same logic as categorical branch)
        n_coeffs <- 1L + 1L  # intercept + numeric contrast
        for (cv in covariates) {
            if (is.factor(coldata[[cv]]) || is.character(coldata[[cv]])) {
                n_coeffs <- n_coeffs + max(1L, nlevels(as.factor(coldata[[cv]])) - 1L)
            } else {
                n_coeffs <- n_coeffs + 1L
            }
        }
        sufficient <- n_tot > n_coeffs

        logger$info("sample_counts",
            sprintf("numeric contrast | total=%d | n_coefficients=%d", n_tot, n_coeffs))

        if (!sufficient) logger$error("sample_counts",
            sprintf("n_samples (%d) <= n_coefficients (%d)", n_tot, n_coeffs))
    }

    list(n_control      = n_ctrl,
         n_experimental = n_exp,
         n_total        = n_tot,
         sufficient     = sufficient)
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id, stratum) {

    # --- Load config and contrast row ---
    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df

    contrast_row <- contrasts_df[contrasts_df$contrast_id == contrast_id, ]
    if (nrow(contrast_row) == 0) {
        stop(sprintf("contrast_id '%s' not found in contrasts file", contrast_id))
    }

    # --- Set up logger ---
    stratum_safe <- sanitize_names(stratum)
    log_file     <- file.path(cfg$outputs$logs_dir, contrast_id,
                               paste0(stratum_safe, ".log"))
    logger       <- make_logger(contrast_id  = contrast_id,
                                biosample    = stratum,
                                log_file     = log_file,
                                level        = cfg$logging$level)

    logger$info("start", sprintf("script=02_prepare_coldata | stratum=%s", stratum))

    # Define output dir early and register skip sentinel
    out_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    register_skip_sentinel(out_dir, c("coldata.csv", "counts_filtered.csv", "name_map.csv"))

    # --- Define output dir early so sentinel can fire on any skip ---
    on.exit({
        for (.f in c("coldata.csv", "counts_filtered.csv", "name_map.csv")) {
            .p <- file.path(out_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
    }, add = TRUE)

    # --- Parse contrast fields ---
    biosample_id_col <- trimws(contrast_row$biosample_id_col)
    contrast_var     <- trimws(contrast_row$contrast_var)
    control_grp      <- trimws(contrast_row$control_grp)
    experimental_grp <- trimws(contrast_row$experimental_grp)
    covariates       <- parse_pipe_field(contrast_row$covariates)
    filters          <- parse_contrast_filters(contrast_row)

    is_categorical   <- !is.na(control_grp) && control_grp != ""

    # --- Load inputs ---
    logger$info("load_counts", sprintf("loading count matrix for stratum: %s", stratum))
    counts_path <- file.path(cfg$inputs$counts_dir,
                             paste0(stratum_safe, ".counts.csv"))

    # Accept .gz too
    if (!file.exists(counts_path) && file.exists(paste0(counts_path, ".gz"))) {
        counts_path <- paste0(counts_path, ".gz")
    }

    raw_mat <- try_logged(
        load_count_matrix(counts_path),
        logger  = logger,
        step    = "load_counts",
        success = sprintf("loaded %s", counts_path)
    )
    if (is.null(raw_mat)) return(invisible(NULL))

    # Separate gene column from counts
    gene_names          <- raw_mat[[1]]
    raw_mat             <- raw_mat[, -1, drop = FALSE]
    colnames(raw_mat)   <- sanitize_sample_id(colnames(raw_mat))
    raw_mat             <- as.data.frame(lapply(raw_mat, as.numeric))
    rownames(raw_mat)   <- gene_names

    logger$info("load_counts",
        sprintf("matrix dimensions: %d genes x %d samples", nrow(raw_mat), ncol(raw_mat)))

    # --- Preflight: validate count matrix data quality ---
    cm_valid <- validate_count_matrix_full(raw_mat)
    for (w in cm_valid$warnings) logger$warn("preflight_counts", w)
    if (length(cm_valid$errors) > 0) {
        for (e in cm_valid$errors) logger$error("preflight_counts", e)
        logger$skip("preflight_counts", "count matrix failed data quality checks — see errors above")
        return(invisible(NULL))
    }

    logger$info("load_metadata", "loading sample metadata")
    sample_metadata <- try_logged(
        fread(cfg$inputs$sample_metadata) %>% as.data.frame(),
        logger  = logger,
        step    = "load_metadata"
    )
    if (!is.null(sample_metadata)) {
        logger$info("load_metadata", sprintf("loaded %d rows", nrow(sample_metadata)))
    }
    if (is.null(sample_metadata)) return(invisible(NULL))

    # --- Filter: restrict to this stratum ---
    # sample_metadata contains rows for all strata — narrow to current stratum first
    # so all subsequent filters and donor dedup operate on the correct samples only.
    # normalize_name() (stats_utils.R) strips non-alnum + lowercases, so
    # "Gamma_Epsilon" matches "Gamma+Epsilon", "ActiveStellate" matches "Active Stellate", etc.
    before_stratum       <- nrow(sample_metadata)
    sample_metadata <- sample_metadata[
        normalize_name(sample_metadata$stratum) == normalize_name(stratum), ]
    logger$info("filter_stratum",
        sprintf("stratum=%s | %d -> %d samples", stratum, before_stratum, nrow(sample_metadata)))

    if (nrow(sample_metadata) == 0) {
        logger$skip("filter_stratum",
            sprintf("no samples found for stratum '%s' after cleaning — check metadata stratum column",
                    stratum))
        return(invisible(NULL))
    }

    # --- Defensive check: critical columns must exist in metadata ---
    # With step 01 soft-gating, a flagged contrast may still reach step 02.
    # Bail out gracefully rather than crashing.
    missing_critical <- character(0)
    if (!biosample_id_col %in% colnames(sample_metadata)) {
        missing_critical <- c(missing_critical,
            sprintf("biosample_id_col '%s' not found in metadata columns", biosample_id_col))
    }
    if (!contrast_var %in% colnames(sample_metadata)) {
        missing_critical <- c(missing_critical,
            sprintf("contrast_var '%s' not found in metadata columns", contrast_var))
    }
    if (length(missing_critical) > 0) {
        for (msg in missing_critical) logger$error("column_check", msg)
        logger$skip("column_check", "critical columns missing from metadata -- cannot proceed")
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        write_preflight_row(
            row_data = list(
                contrast_id = contrast_id, stratum = stratum, status = "fail",
                errors = paste(missing_critical, collapse = " | ")
            ),
            path = file.path(out_dir, "preflight.csv")
        )
        return(invisible(NULL))
    }

    # Canonicalize sample IDs to match sanitized count-matrix column names
    sample_metadata[[biosample_id_col]] <- sanitize_sample_id(sample_metadata[[biosample_id_col]])

    # --- Filter: minimum cells per sample ---
    if (!is.null(cfg$filtering$min_cells_per_sample)) {
        min_cells <- cfg$filtering$min_cells_per_sample

        # match_stratum_column() (stats_utils.R) resolves naming mismatches:
        #   "ActiveStellate" -> "Active Stellate", "Immune" -> "Immune (Macrophages)", etc.
        col_match <- match_stratum_column(stratum, colnames(sample_metadata))

        if (!is.null(col_match$col_name)) {
            if (col_match$ambiguous) {
                logger$warn("filter_min_cells",
                    sprintf("multiple columns match '%s': [%s] — using '%s' (%s match)",
                            stratum,
                            paste(col_match$candidates, collapse = ", "),
                            col_match$col_name, col_match$method))
            }
            stratum_col <- col_match$col_name
            before          <- nrow(sample_metadata)
            sample_metadata <- sample_metadata[
                !is.na(sample_metadata[[stratum_col]]) &
                sample_metadata[[stratum_col]] >= min_cells, ]
            logger$info("filter_min_cells",
                sprintf("min_cells=%d (col='%s', %s match) | %d -> %d samples",
                        min_cells, stratum_col, col_match$method, before, nrow(sample_metadata)))
        } else {
            logger$warn("filter_min_cells",
                sprintf("no cell count column found for stratum '%s' — skipping min_cells filter",
                        stratum))
        }
    }

    # --- Filter: restrict to contrast groups ---
    if (is_categorical) {
        before          <- nrow(sample_metadata)
        sample_metadata <- sample_metadata[
            sample_metadata[[contrast_var]] %in% c(control_grp, experimental_grp), ]
        logger$info("filter_groups",
            sprintf("kept groups [%s, %s] | %d -> %d samples",
                    control_grp, experimental_grp, before, nrow(sample_metadata)))
    }

    # --- Filter: user-defined subset filters ---
    if (length(filters) > 0) {
        logger$info("filter_subsets", sprintf("applying %d user-defined filter(s)", length(filters)))
        sample_metadata <- apply_subset_filters(sample_metadata, filters)
        logger$info("filter_subsets",
            sprintf("after subset filters: %d samples remain", nrow(sample_metadata)))
    }

    # --- Filter: paired analysis ---
    if (isTRUE(cfg$steps$paired_analysis) && is_categorical) {
        logger$info("paired_filter", "paired_analysis enabled — filtering to matched donors")
        sample_metadata <- apply_paired_filter(
            meta             = sample_metadata,
            contrast_var     = contrast_var,
            control_grp      = control_grp,
            experimental_grp = experimental_grp,
            donor_col        = cfg$filtering$donor_col,
            logger           = logger
        )
        if (is.null(sample_metadata)) return(invisible(NULL))
    }

    # --- Deduplicate donors ---
    # Unpaired: keep 1 sample per donor (random, deterministic via set.seed).
    # Paired:   skip dedup — multiple dishes per donor per treatment are genuine
    #           biological replicates. ~ donor + treatment handles the blocking.
    sample_metadata <- try_logged(
        dedup_donors(sample_metadata,
                     donor_col = cfg$filtering$donor_col,
                     paired    = isTRUE(cfg$steps$paired_analysis)),
        logger  = logger,
        step    = "dedup_donors",
        success = sprintf("%d unique donors retained", nrow(
            unique(sample_metadata[cfg$filtering$donor_col])))
    )
    if (is.null(sample_metadata)) return(invisible(NULL))

    # --- Align count matrix columns to surviving samples ---
    shared_samples  <- intersect(colnames(raw_mat), sample_metadata[[biosample_id_col]])
    missing_samples <- setdiff(sample_metadata[[biosample_id_col]], colnames(raw_mat))

    if (length(missing_samples) > 0) {
        logger$warn("align_samples",
            sprintf("%d samples in metadata not found in count matrix — dropping: %s",
                    length(missing_samples),
                    paste(head(missing_samples, 5), collapse = ", ")))
    }

    sample_metadata <- sample_metadata[sample_metadata[[biosample_id_col]] %in% shared_samples, ]
    raw_mat         <- raw_mat[, shared_samples, drop = FALSE]

    logger$info("align_samples",
        sprintf("aligned: %d samples | %d genes", ncol(raw_mat), nrow(raw_mat)))

    # --- Minimum sample check before gene filtering ---
    if (nrow(sample_metadata) < 3) {
        logger$skip("insufficient_samples",
            sprintf("only %d samples remain after filtering — minimum is 3", nrow(sample_metadata)))
        return(invisible(NULL))
    }

    # --- Gene filter ---
    logger$info("gene_filter", "applying group-aware gene filter")
    genes_to_keep <- group_aware_gene_filter(
        raw_mat          = as.matrix(raw_mat),
        sample_metadata  = sample_metadata,
        biosample_id_col = biosample_id_col,
        contrast_var     = contrast_var,
        control_grp      = if (is_categorical) control_grp else NULL,
        experimental_grp = if (is_categorical) experimental_grp else NULL,
        min_reads        = cfg$filtering$min_gene_counts
    )
    raw_mat <- raw_mat[genes_to_keep, , drop = FALSE]

    logger$info("gene_filter",
        sprintf("genes retained: %d / %d", length(genes_to_keep),
                length(genes_to_keep) + sum(!rownames(raw_mat) %in% genes_to_keep)))

    # --- Preflight: zero-inflation check ---
    zi_result <- check_zero_inflation(as.matrix(raw_mat))
    if (zi_result$zero_inflated) {
        logger$warn("preflight_zero_inflation",
            sprintf("zero-inflated: 100%% of genes have >= 1 zero — DESeq2 will use poscounts"))
    } else {
        logger$info("preflight_zero_inflation",
            sprintf("%.1f%% of genes have >= 1 zero", zi_result$pct_genes_with_zero))
    }

    # --- Prepare coldata ---
    coldata_result <- prepare_coldata(
        meta             = sample_metadata,
        biosample_id_col = biosample_id_col,
        contrast_var     = contrast_var,
        covariates       = covariates,
        scale_numeric    = isTRUE(cfg$deseq2$scale_numeric_covariates)
    )

    coldata          <- coldata_result$coldata
    name_map         <- coldata_result$name_map
    dropped_covs     <- coldata_result$dropped_covariates

    if (length(dropped_covs) > 0) {
        logger$warn("prepare_coldata",
            sprintf("covariates dropped (only 1 unique value): %s",
                    paste(dropped_covs, collapse = ", ")))
    }

    # --- Second-pass covariate drop: factor levels after NA removal ---
    # The first pass (above) drops covariates with < 2 unique raw values.
    # But a factor with 2 levels may effectively have only 1 usable level if
    # the minority level has 0 non-NA observations for other covariates.
    # DESeq2 calls contrasts() internally and crashes with:
    #   "contrasts can be applied only to factors with 2 or more levels"
    # Catch this BEFORE running DESeq2.
    active_covs_2 <- setdiff(sanitize_names(covariates), dropped_covs)
    for (v in active_covs_2) {
        if (v %in% colnames(coldata) && is.factor(coldata[[v]])) {
            # droplevels to see actual observed levels
            observed_levels <- levels(droplevels(coldata[[v]]))
            if (length(observed_levels) < 2) {
                dropped_covs <- c(dropped_covs, v)
                logger$warn("drop_covariate",
                    sprintf("%s has < 2 observed factor levels after filtering — dropping", v))
            }
        }
    }

    # --- Preflight: design matrix rank and near-singularity check ---
    safe_contrast_pf <- sanitize_names(contrast_var)
    active_covs_pf   <- setdiff(sanitize_names(covariates), dropped_covs)
    preflight_warnings <- character(0)
    preflight_errors   <- character(0)
    preflight_rank_ok  <- TRUE
    preflight_near_singular <- character(0)

    dm_result <- validate_design_matrix(
        contrast_var = safe_contrast_pf,
        covariates   = active_covs_pf,
        meta         = coldata,
        control_grp  = if (is_categorical) control_grp else NULL,
        experimental_grp = if (is_categorical) experimental_grp else NULL
    )

    preflight_rank_ok      <- dm_result$rank_ok
    preflight_near_singular <- dm_result$near_singular_covariates

    for (w in dm_result$warnings) {
        logger$warn("preflight_design", w)
        preflight_warnings <- c(preflight_warnings, w)
    }
    for (e in dm_result$errors) {
        logger$error("preflight_design", e)
        preflight_errors <- c(preflight_errors, e)
    }

    # --- Final sample count check ---
    safe_contrast    <- sanitize_names(contrast_var)
    active_covs      <- setdiff(sanitize_names(covariates), dropped_covs)
    counts_result    <- check_sample_counts(
        coldata          = coldata,
        contrast_var     = safe_contrast,
        control_grp      = control_grp,
        experimental_grp = experimental_grp,
        covariates       = active_covs,
        logger           = logger
    )

    if (!counts_result$sufficient) {
        logger$skip("insufficient_samples",
            "not enough samples for the model — see errors above")
        return(invisible(NULL))
    }

    # --- Write outputs ---
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(cfg$write_outputs$coldata_csv)) {
        coldata_path <- file.path(out_dir, "coldata.csv")
        write_result(coldata, coldata_path, sep = ",")
        logger$info("write_outputs", sprintf("coldata: %s", coldata_path))
    }

    counts_path_out <- file.path(out_dir, "counts_filtered.csv")
    write_result(
        cbind(gene = rownames(raw_mat), raw_mat),
        counts_path_out,
        sep = ","
    )
    logger$info("write_outputs", sprintf("counts_filtered: %s", counts_path_out))

    # Save name_map as a small lookup CSV for downstream scripts
    name_map_df   <- data.frame(original = names(name_map), sanitized = name_map,
                                stringsAsFactors = FALSE)
    write_result(name_map_df, file.path(out_dir, "name_map.csv"), sep = ",")

    # --- Write preflight report ---
    preflight_status <- if (length(preflight_errors) > 0) {
        "fail"
    } else if (length(preflight_warnings) > 0 || length(dropped_covs) > 0) {
        "warn"
    } else {
        "pass"
    }

    write_preflight_row(
        row_data = list(
            contrast_id         = contrast_id,
            stratum             = stratum,
            status              = preflight_status,
            n_genes             = nrow(raw_mat),
            n_samples           = counts_result$n_total,
            n_control           = counts_result$n_control,
            n_experimental      = counts_result$n_experimental,
            n_covariates_active = length(active_covs),
            dropped_covariates  = if (length(dropped_covs) > 0) paste(dropped_covs, collapse = "|") else NA_character_,
            rank_ok             = preflight_rank_ok,
            near_singular_vars  = if (length(preflight_near_singular) > 0) paste(preflight_near_singular, collapse = "|") else NA_character_,
            zero_inflated       = zi_result$zero_inflated,
            pct_genes_with_zero = zi_result$pct_genes_with_zero,
            warnings            = if (length(preflight_warnings) > 0) paste(preflight_warnings, collapse = "|") else NA_character_,
            errors              = if (length(preflight_errors) > 0) paste(preflight_errors, collapse = "|") else NA_character_
        ),
        path = file.path(out_dir, "preflight.csv")
    )
    logger$info("write_outputs", sprintf("preflight.csv: status=%s", preflight_status))

    logger$info("complete",
        sprintf("n_control=%s | n_experimental=%s | n_total=%d | n_genes=%d",
                counts_result$n_control, counts_result$n_experimental,
                counts_result$n_total, nrow(raw_mat)))

    invisible(list(
        coldata     = coldata,
        counts      = raw_mat,
        name_map    = name_map,
        n_control   = counts_result$n_control,
        n_experimental = counts_result$n_experimental
    ))
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {

    opt_list <- list(
        make_option("--config",   type = "character", default = "config/config.yaml"),
        make_option("--contrast", type = "character", help = "contrast_id from contrasts.csv"),
        make_option("--stratum",  type = "character", help = "stratum name (e.g. Beta)")
    )

    opts <- parse_args(OptionParser(option_list = opt_list))
    if (is.null(opts$contrast)) stop("--contrast is required")
    if (is.null(opts$stratum))  stop("--stratum is required")

    main(opts$config, opts$contrast, opts$stratum)
}