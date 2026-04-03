# =============================================================================
# validation_utils.R
# EasyDE — Reusable Preflight Validation Functions
#
# Pure validation functions with no side effects. Each returns a named list
# with $errors (character vector) and $warnings (character vector), plus any
# additional metadata. Callers decide how to log or act on the results.
#
# These are used by 01_prepare_inputs.R (structural preflight) and
# 02_prepare_coldata.R (data-level preflight).
# =============================================================================

# Needs sanitize_names() from stats_utils.R — sourced by the calling scripts


# =============================================================================
# Structural validators (used in 01_prepare_inputs.R — before any computation)
# =============================================================================

#' Quick structural validation of a count matrix file
#'
#' Reads only the first few rows to check structure without loading the full
#' matrix. Validates: readable, has genes, has samples, no duplicate column
#' names, first column is character (gene names), remaining columns are numeric.
#'
#' @param path  Path to the count matrix file (.csv or .csv.gz)
#' @param stratum_name  Name of the stratum (for error messages)
#' @return Named list: $errors, $warnings, $n_samples (from header)
validate_count_matrix <- function(path, stratum_name = "") {

    errors   <- character(0)
    warnings <- character(0)
    n_samples <- NA_integer_

    # Try to read first 5 rows
    mat <- tryCatch(
        suppressWarnings(data.table::fread(path, nrows = 5, check.names = FALSE)),
        error = function(e) NULL
    )

    if (is.null(mat)) {
        errors <- c(errors, sprintf("stratum '%s': cannot read count matrix: %s",
                                    stratum_name, path))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    if (ncol(mat) < 2) {
        errors <- c(errors, sprintf("stratum '%s': count matrix has fewer than 2 columns (need gene names + at least 1 sample)",
                                    stratum_name))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    if (nrow(mat) == 0) {
        errors <- c(errors, sprintf("stratum '%s': count matrix has 0 rows (no genes)",
                                    stratum_name))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    # Check for duplicate sample column names (excluding gene column)
    sample_cols <- colnames(mat)[-1]
    n_samples   <- length(sample_cols)
    dup_cols    <- sample_cols[duplicated(sample_cols)]
    if (length(dup_cols) > 0) {
        errors <- c(errors, sprintf("stratum '%s': duplicate column names in count matrix: %s",
                                    stratum_name, paste(unique(dup_cols), collapse = ", ")))
    }

    # Check that first column looks like gene names (character)
    if (is.numeric(mat[[1]])) {
        warnings <- c(warnings, sprintf("stratum '%s': first column of count matrix is numeric (expected gene names)",
                                        stratum_name))
    }

    # Check that remaining columns are numeric
    non_numeric <- character(0)
    for (j in 2:ncol(mat)) {
        if (!is.numeric(mat[[j]])) {
            non_numeric <- c(non_numeric, colnames(mat)[j])
        }
    }
    if (length(non_numeric) > 0) {
        errors <- c(errors, sprintf("stratum '%s': non-numeric data columns: %s",
                                    stratum_name, paste(head(non_numeric, 5), collapse = ", ")))
    }

    list(errors = errors, warnings = warnings, n_samples = n_samples)
}


#' Validate a metadata column: existence, completeness, and type
#'
#' @param meta           data.frame of sample metadata
#' @param col_name       Column name to validate
#' @param expected_type  Optional: "numeric", "character", or NULL for any
#' @param context        String for error messages (e.g. "contrast 'T1D_vs_ND'")
#' @return Named list: $errors, $warnings, $n_na, $n_unique
validate_metadata_column <- function(meta, col_name, expected_type = NULL, context = "") {

    errors   <- character(0)
    warnings <- character(0)
    n_na     <- NA_integer_
    n_unique <- NA_integer_
    prefix   <- if (nchar(context) > 0) paste0(context, ": ") else ""

    if (!col_name %in% colnames(meta)) {
        errors <- c(errors, sprintf("%scolumn '%s' not found in metadata", prefix, col_name))
        return(list(errors = errors, warnings = warnings, n_na = n_na, n_unique = n_unique))
    }

    vals     <- meta[[col_name]]
    n_na     <- sum(is.na(vals))
    n_unique <- length(unique(na.omit(vals)))

    if (all(is.na(vals))) {
        errors <- c(errors, sprintf("%scolumn '%s' is entirely NA", prefix, col_name))
    } else if (n_na > 0) {
        na_pct <- round(100 * n_na / length(vals), 1)
        if (na_pct > 50) {
            warnings <- c(warnings, sprintf("%scolumn '%s' has %s%% NA values (%d/%d)",
                                            prefix, col_name, na_pct, n_na, length(vals)))
        }
    }

    if (!is.null(expected_type) && !all(is.na(vals))) {
        actual_type <- if (is.numeric(vals)) "numeric" else "character"
        if (actual_type != expected_type) {
            warnings <- c(warnings, sprintf("%scolumn '%s' expected %s but is %s",
                                            prefix, col_name, expected_type, actual_type))
        }
    }

    list(errors = errors, warnings = warnings, n_na = n_na, n_unique = n_unique)
}


#' Validate that contrast group values exist in the data
#'
#' @param meta          data.frame of sample metadata
#' @param contrast_var  Column name of the contrast variable
#' @param control_grp   Value defining the control group
#' @param experimental_grp Value defining the experimental group
#' @param context       String for error messages
#' @return Named list: $errors, $warnings, $n_control, $n_experimental
validate_group_membership <- function(meta, contrast_var, control_grp, experimental_grp, context = "") {

    errors   <- character(0)
    warnings <- character(0)
    n_control <- NA_integer_
    n_experimental <- NA_integer_
    prefix <- if (nchar(context) > 0) paste0(context, ": ") else ""

    if (!contrast_var %in% colnames(meta)) {
        errors <- c(errors, sprintf("%scontrast_var '%s' not in metadata", prefix, contrast_var))
        return(list(errors = errors, warnings = warnings,
                    n_control = n_control, n_experimental = n_experimental))
    }

    vals <- meta[[contrast_var]]
    n_control      <- sum(vals == control_grp, na.rm = TRUE)
    n_experimental <- sum(vals == experimental_grp, na.rm = TRUE)

    if (n_control == 0) {
        errors <- c(errors, sprintf("%scontrol_grp '%s' not found in column '%s' (0 matches)",
                                    prefix, control_grp, contrast_var))
    }
    if (n_experimental == 0) {
        errors <- c(errors, sprintf("%sexperimental_grp '%s' not found in column '%s' (0 matches)",
                                    prefix, experimental_grp, contrast_var))
    }

    list(errors = errors, warnings = warnings,
         n_control = n_control, n_experimental = n_experimental)
}


#' Validate a pipe-separated field for empty tokens
#'
#' Catches malformed fields like "Sex||Age" or "|Sex" or "Sex|"
#'
#' @param field       String value from contrasts CSV
#' @param field_name  Name of the field (for error messages)
#' @param context     String for error messages
#' @return Named list: $errors, $warnings, $parsed_values
validate_pipe_field <- function(field, field_name, context = "") {

    errors   <- character(0)
    warnings <- character(0)
    prefix   <- if (nchar(context) > 0) paste0(context, ": ") else ""

    if (is.na(field) || trimws(field) == "") {
        return(list(errors = errors, warnings = warnings, parsed_values = character(0)))
    }

    tokens <- strsplit(field, "\\|")[[1]]
    trimmed <- trimws(tokens)

    empty_count <- sum(trimmed == "")
    if (empty_count > 0) {
        warnings <- c(warnings, sprintf("%s%s has %d empty token(s) after pipe-splitting: '%s'",
                                        prefix, field_name, empty_count, field))
    }

    parsed <- trimmed[trimmed != ""]
    list(errors = errors, warnings = warnings, parsed_values = parsed)
}


#' Check for sanitized name collisions
#'
#' When variable names are sanitized for R formulas, distinct names might
#' collide (e.g. "T1D-status" and "T1D status" both become "T1D_status").
#'
#' @param names_vec  Character vector of original names
#' @param context    String for error messages
#' @return Named list: $errors, $warnings, $name_map (original -> sanitized)
validate_sanitized_names <- function(names_vec, context = "") {

    errors   <- character(0)
    warnings <- character(0)
    prefix   <- if (nchar(context) > 0) paste0(context, ": ") else ""

    if (length(names_vec) == 0) {
        return(list(errors = errors, warnings = warnings, name_map = character(0)))
    }

    sanitized <- sanitize_names(names_vec)
    name_map  <- setNames(sanitized, names_vec)

    # Check for collisions
    dup_sanitized <- sanitized[duplicated(sanitized)]
    if (length(dup_sanitized) > 0) {
        for (d in unique(dup_sanitized)) {
            originals <- names_vec[sanitized == d]
            errors <- c(errors, sprintf(
                "%snames collide after sanitization: [%s] all map to '%s'",
                prefix, paste(originals, collapse = ", "), d))
        }
    }

    list(errors = errors, warnings = warnings, name_map = name_map)
}


#' Validate donor column for NAs and existence
#'
#' @param meta       data.frame of sample metadata
#' @param donor_col  Column name for the donor identifier
#' @param context    String for error messages
#' @return Named list: $errors, $warnings, $n_na_donors
validate_donor_column <- function(meta, donor_col, context = "") {

    errors   <- character(0)
    warnings <- character(0)
    n_na     <- NA_integer_
    prefix   <- if (nchar(context) > 0) paste0(context, ": ") else ""

    if (is.null(donor_col) || is.na(donor_col) || trimws(donor_col) == "") {
        return(list(errors = errors, warnings = warnings, n_na_donors = n_na))
    }

    if (!donor_col %in% colnames(meta)) {
        errors <- c(errors, sprintf("%sdonor_col '%s' not found in metadata", prefix, donor_col))
        return(list(errors = errors, warnings = warnings, n_na_donors = n_na))
    }

    n_na <- sum(is.na(meta[[donor_col]]))
    if (n_na > 0) {
        warnings <- c(warnings, sprintf(
            "%sdonor_col '%s' has %d NA value(s) — these samples may be grouped incorrectly during deduplication",
            prefix, donor_col, n_na))
    }

    list(errors = errors, warnings = warnings, n_na_donors = n_na)
}


# =============================================================================
# Data-level validators (used in 02_prepare_coldata.R — on filtered data)
# =============================================================================

#' Validate design matrix rank and detect near-singular covariates
#'
#' Builds a trial model.matrix from the formula and checks:
#'   1. Rank deficiency (rank < ncol) — error
#'   2. Near-singular categorical covariates (>=95% one level) — warning
#'   3. Near-singular numeric covariates (CV < 0.01) — warning
#'
#' @param contrast_var  Sanitized contrast variable name
#' @param covariates    Character vector of sanitized covariate names (active ones)
#' @param meta          data.frame (coldata) with sanitized column names
#' @param control_grp   Control group value (or NULL for numeric contrasts)
#' @param experimental_grp Experimental group value (or NULL)
#' @return Named list: $errors, $warnings, $rank, $n_coefficients,
#'         $near_singular_covariates, $rank_ok
validate_design_matrix <- function(contrast_var, covariates, meta,
                                   control_grp = NULL, experimental_grp = NULL) {

    errors   <- character(0)
    warnings <- character(0)
    near_singular <- character(0)
    rank_val <- NA_integer_
    n_coefs  <- NA_integer_
    rank_ok  <- NA

    # Build formula
    terms <- c(covariates, contrast_var)
    formula_str <- paste0("~ ", paste(terms, collapse = " + "))

    mm <- tryCatch(
        stats::model.matrix(as.formula(formula_str), data = meta),
        error = function(e) NULL
    )

    if (is.null(mm)) {
        errors <- c(errors, sprintf("cannot build design matrix from formula: %s", formula_str))
        return(list(errors = errors, warnings = warnings, rank = rank_val,
                    n_coefficients = n_coefs, near_singular_covariates = near_singular,
                    rank_ok = FALSE))
    }

    n_coefs  <- ncol(mm)
    rank_val <- qr(mm)$rank
    rank_ok  <- (rank_val == n_coefs)

    if (!rank_ok) {
        errors <- c(errors, sprintf(
            "design matrix is rank-deficient (rank=%d, n_coefficients=%d, formula=%s) — covariates may be collinear",
            rank_val, n_coefs, formula_str))
    }

    # Check for near-singular covariates
    for (v in covariates) {
        if (!v %in% colnames(meta)) next
        col_vals <- meta[[v]]

        if (is.factor(col_vals) || is.character(col_vals)) {
            # Categorical: check if one level has >= 95% of observations
            tbl <- table(col_vals, useNA = "no")
            if (length(tbl) > 0) {
                max_prop <- max(tbl) / sum(tbl)
                if (max_prop >= 0.95) {
                    dominant <- names(tbl)[which.max(tbl)]
                    near_singular <- c(near_singular, v)
                    warnings <- c(warnings, sprintf(
                        "covariate '%s' is near-singular: level '%s' has %.1f%% of observations",
                        v, dominant, max_prop * 100))
                }
            }
        } else if (is.numeric(col_vals)) {
            # Numeric: check coefficient of variation
            col_clean <- na.omit(col_vals)
            if (length(col_clean) > 1 && mean(col_clean) != 0) {
                cv <- sd(col_clean) / abs(mean(col_clean))
                if (cv < 0.01) {
                    near_singular <- c(near_singular, v)
                    warnings <- c(warnings, sprintf(
                        "covariate '%s' has near-zero variation (CV=%.4f)",
                        v, cv))
                }
            }
        }
    }

    list(errors = errors, warnings = warnings, rank = rank_val,
         n_coefficients = n_coefs, near_singular_covariates = near_singular,
         rank_ok = rank_ok)
}


#' Validate an already-loaded count matrix for data quality issues
#'
#' Checks for negative values (error) and all-zero gene rows (warning).
#' Called in 02_prepare_coldata.R after the full matrix is loaded.
#'
#' @param mat  Numeric matrix (genes x samples), already loaded
#' @return Named list: $errors, $warnings, $n_negative, $n_zero_rows
validate_count_matrix_full <- function(mat) {

    errors   <- character(0)
    warnings <- character(0)
    n_negative  <- 0L
    n_zero_rows <- 0L

    if (!is.matrix(mat) && !is.data.frame(mat)) {
        return(list(errors = errors, warnings = warnings,
                    n_negative = n_negative, n_zero_rows = n_zero_rows))
    }

    # Check for negative values
    num_mat <- as.matrix(mat)
    storage.mode(num_mat) <- "double"
    neg_mask <- num_mat < 0 & !is.na(num_mat)
    n_negative <- sum(neg_mask)
    if (n_negative > 0) {
        errors <- c(errors, sprintf(
            "count matrix contains %d negative value(s) — DESeq2 requires non-negative integers",
            n_negative))
    }

    # Check for all-zero gene rows
    row_sums <- rowSums(num_mat, na.rm = TRUE)
    n_zero_rows <- sum(row_sums == 0)
    if (n_zero_rows > 0) {
        warnings <- c(warnings, sprintf(
            "%d gene(s) have all-zero counts across all samples — these will have no power",
            n_zero_rows))
    }

    list(errors = errors, warnings = warnings,
         n_negative = n_negative, n_zero_rows = n_zero_rows)
}


#' Check whether a count matrix is zero-inflated for DESeq2 size factor estimation
#'
#' DESeq2's default estimateSizeFactors(type="ratio") computes a per-gene
#' geometric mean, which requires at least one gene where NO sample is zero.
#' If every gene has >= 1 zero, the geometric mean is 0 for every gene and
#' DESeq2 fails. The standard fix is type="poscounts".
#'
#' @param mat  Numeric matrix (genes x samples)
#' @return Named list: $zero_inflated (logical), $pct_genes_with_zero (numeric)
check_zero_inflation <- function(mat) {

    if (!is.matrix(mat) && !is.data.frame(mat)) {
        return(list(zero_inflated = FALSE, pct_genes_with_zero = NA_real_))
    }

    num_mat <- as.matrix(mat)
    genes_with_zero <- apply(num_mat, 1, function(x) any(x == 0, na.rm = TRUE))
    pct <- 100 * sum(genes_with_zero) / length(genes_with_zero)

    list(
        zero_inflated       = all(genes_with_zero),
        pct_genes_with_zero = round(pct, 1)
    )
}


# =============================================================================
# Preflight report I/O
# =============================================================================

#' Write a single preflight report row for one stratum
#'
#' Writes a one-row CSV to the stratum's intermediates directory.
#' Each stratum gets its own file (parallel-safe under Snakemake).
#'
#' @param row_data  Named list with preflight result fields
#' @param path      Output path for the preflight CSV
write_preflight_row <- function(row_data, path) {

    df <- data.frame(
        contrast_id           = row_data$contrast_id           %||% NA_character_,
        stratum               = row_data$stratum                %||% NA_character_,
        status                = row_data$status                 %||% "unknown",
        n_genes               = row_data$n_genes                %||% NA_integer_,
        n_samples             = row_data$n_samples              %||% NA_integer_,
        n_control             = row_data$n_control              %||% NA_integer_,
        n_experimental        = row_data$n_experimental         %||% NA_integer_,
        n_covariates_active   = row_data$n_covariates_active    %||% NA_integer_,
        dropped_covariates    = row_data$dropped_covariates     %||% NA_character_,
        rank_ok               = row_data$rank_ok                %||% NA,
        near_singular_vars    = row_data$near_singular_vars     %||% NA_character_,
        zero_inflated         = row_data$zero_inflated          %||% NA,
        pct_genes_with_zero   = row_data$pct_genes_with_zero    %||% NA_real_,
        warnings              = row_data$warnings               %||% NA_character_,
        errors                = row_data$errors                 %||% NA_character_,
        timestamp             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        stringsAsFactors      = FALSE
    )

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(df, path, row.names = FALSE, quote = TRUE)
}
