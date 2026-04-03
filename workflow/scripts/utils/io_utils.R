# =============================================================================
# io_utils.R
# Input/output helpers for EasyDE
#
# General-purpose functions for reading, caching, and writing data.
# No project-specific logic lives here — this works for any dataset.
# =============================================================================

suppressMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
})


#' Load a file from a local path OR a URL, with optional local caching
#'
#' If `cache_path` is provided and the file already exists locally,
#' it is read from disk (fast). If not, the file is fetched from `source`
#' (URL or path) and saved to `cache_path` for future runs.
#'
#' This prevents re-downloading metadata on every pipeline run.
#'
#' @param source      File path or URL to read from
#' @param cache_path  Local path to cache the file (NULL = no caching)
#' @param ...         Additional arguments passed to fread()
#' @return data.frame
load_with_cache <- function(source, cache_path = NULL, ...) {

    if (!is.null(cache_path) && file.exists(cache_path)) {
        message(sprintf("  [cache hit]  Reading from local cache: %s", cache_path))
        return(fread(cache_path, ...) %>% as.data.frame())
    }

    message(sprintf("  [fetch]      Reading: %s", source))
    df <- fread(source, ...) %>% as.data.frame()

    if (!is.null(cache_path)) {
        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        fwrite(df, cache_path)
        message(sprintf("  [cache save] Saved to: %s", cache_path))
    }

    return(df)
}


#' Load a count matrix from a local file, return NULL on failure
#'
#' Wraps fread() in tryCatch so callers get NULL instead of a crash.
#' The calling script decides what to do (log the skip, continue the loop).
#'
#' @param path  Path or URL to the count matrix file
#' @return data.frame or NULL if loading failed
load_count_matrix <- function(path) {

    # The PanKbase .fixed.txt files have no column header for the gene-name
    # column, so fread warns "Detected N column names but data has N+1 columns"
    # and auto-names it "V1". This is correct behavior — we suppress that
    # specific warning because it is expected and harmless.
    result <- tryCatch(
        withCallingHandlers(
            fread(path, check.names = FALSE) %>% as.data.frame(),
            warning = function(w) {
                if (grepl("column names but the data has", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        error = function(e) {
            message(sprintf("  [skip] Could not load matrix from: %s\n  Reason: %s",
                            path, e$message))
            return(NULL)
        }
    )

    return(result)
}


#' Build a contrast_summary row for one biosample (success or failure)
#'
#' Centralizes the repeated data.frame() construction that previously appeared
#' 6+ times in the original script. Pass NA for fields not yet determined.
#'
#' @param biosample            The stratum being analyzed — one biosample group
#'                             such as a cell type or tissue (e.g. "Beta", "islet").
#'                             This is the strata loop variable from contrasts.csv,
#'                             NOT the biosample_id_col (donor/subject unit).
#' @param contrast_id          Contrast identifier string
#' @param status               One of: "success", "skip", "error"
#' @param error_message        Human-readable reason for skip/error (or NA)
#' @param n_control_samps      Number of control samples (or NA)
#' @param n_experimental_samps Number of experimental samples (or NA)
#' @param n_samples_total      Total samples used (or NA)
#' @param base_formula         DESeq2 base formula string (or NA)
#' @param n_degs_base          Number of DEGs from base model (or NA)
#' @param ruv_formula          Final RUV formula string (or NA)
#' @param n_degs_ruv           Number of DEGs from RUV model (or NA)
#' @param n_pathways_base      Number of significant fGSEA pathways on base results (or NA)
#' @param n_pathways_ruv       Number of significant fGSEA pathways on RUV results (or NA)
#' @return Single-row data.frame
make_summary_row <- function(biosample,
                             contrast_id,
                             status,
                             error_message          = NA,
                             n_control_samps        = NA,
                             n_experimental_samps   = NA,
                             n_samples_total        = NA,
                             base_formula           = NA,
                             n_degs_base            = NA,
                             ruv_formula            = NA,
                             n_degs_ruv             = NA,
                             n_pathways_base        = NA,
                             n_pathways_ruv         = NA) {
    data.frame(
        biosample              = biosample,
        contrast               = contrast_id,
        status                 = status,
        n_control_samps        = n_control_samps,
        n_experimental_samps   = n_experimental_samps,
        n_samples_total        = n_samples_total,
        base_formula           = base_formula,
        n_degs_base            = n_degs_base,
        ruv_formula            = ruv_formula,
        n_degs_ruv             = n_degs_ruv,
        n_pathways_base        = n_pathways_base,
        n_pathways_ruv         = n_pathways_ruv,
        error_message          = error_message,
        stringsAsFactors       = FALSE
    )
}


#' Write a TSV or CSV result file, creating the output directory if needed
#'
#' @param df        data.frame to write
#' @param path      Output file path
#' @param sep       Separator: "\t" for TSV (default), "," for CSV
#' @param row_names Whether to include row names (default FALSE)
write_result <- function(df, path, sep = "\t", row_names = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(df, path, sep = sep, quote = TRUE, row.names = row_names)
    message(sprintf("  [write] %s", path))
}


#' Parse a pipe-separated field from the contrasts CSV into a character vector
#'
#' Handles empty/NA fields gracefully — returns an empty character vector.
#' Used for covariates, strata, filter values, etc.
#'
#' @param field  Single string, potentially pipe-separated (e.g. "Sex|BMI|Age")
#' @return Character vector, or character(0) if field is empty/NA
parse_pipe_field <- function(field) {
    if (is.na(field) || trimws(field) == "") return(character(0))
    trimws(strsplit(field, "\\|")[[1]])
}


#' Parse filter columns and values from a contrast row into a named list
#'
#' Reads all filter_col_N / filter_val_N column pairs from a contrast row
#' and returns a named list suitable for apply_subset_filters().
#'
#' The pipeline discovers filter pairs dynamically using grep() — the user
#' simply adds filter_col_3/filter_val_3, filter_col_4/filter_val_4, etc.
#' columns to the contrasts CSV and they are picked up automatically.
#' Pairs where either the column name or value is empty/NA are silently skipped.
#'
#' @param contrast_row  Single-row data.frame from the contrasts CSV
#' @return Named list: column_name -> character vector of allowed values
#'         Returns empty list if no filter columns exist or all are empty.
parse_contrast_filters <- function(contrast_row) {

    filters <- list()

    # Discover all filter_col_N columns dynamically — works for any N
    filter_col_names <- grep("^filter_col_", colnames(contrast_row), value = TRUE)

    for (col_name in filter_col_names) {
        n       <- sub("^filter_col_", "", col_name)
        val_col <- paste0("filter_val_", n)

        col_val <- contrast_row[[col_name]]
        val_val <- contrast_row[[val_col]]

        # Skip if either side is empty/NA
        if (is.na(col_val) || trimws(col_val) == "") next
        if (is.na(val_val) || trimws(val_val) == "") next

        filters[[trimws(col_val)]] <- parse_pipe_field(val_val)
    }

    return(filters)
}


#' Load pipeline config and contrasts table
#'
#' Centralizes the repeated config + contrasts loading that appears in every
#' pipeline script. Suppresses the "incomplete final line" warning that
#' readLines and yaml::read_yaml emit when files lack a trailing newline.
#'
#' @param config_path  Path to config.yaml
#' @return Named list: \$cfg (parsed config), \$contrasts_df (data.frame)
load_config <- function(config_path) {

    cfg <- suppressWarnings(yaml::read_yaml(config_path))

    raw_lines  <- suppressWarnings(readLines(cfg$inputs$contrasts_file))
    data_lines <- raw_lines[!grepl("^#", raw_lines) & trimws(raw_lines) != ""]
    contrasts_df <- fread(text = paste(data_lines, collapse = "\n")) %>% as.data.frame()

    list(cfg = cfg, contrasts_df = contrasts_df)
}

# =============================================================================
# Preflight status reader
# =============================================================================

#' Read the preflight status for a given stratum
#'
#' Reads the per-stratum preflight.csv written by 02_prepare_coldata.R.
#' Returns "unknown" if the file doesn't exist (backward compatibility
#' with results generated before preflight reporting was added).
#'
#' @param inter_dir  Path to the stratum's intermediates directory
#' @return Character: "pass", "warn", "fail", or "unknown"
read_preflight_status <- function(inter_dir) {
    path <- file.path(inter_dir, "preflight.csv")
    if (!file.exists(path)) return("unknown")

    tryCatch({
        df <- read.csv(path, stringsAsFactors = FALSE)
        if (nrow(df) > 0 && "status" %in% colnames(df)) {
            return(df$status[1])
        }
        "unknown"
    }, error = function(e) "unknown")
}


#' Read the full preflight report for a given stratum
#'
#' Returns the complete preflight data.frame, or NULL if not available.
#'
#' @param inter_dir  Path to the stratum's intermediates directory
#' @return data.frame or NULL
read_preflight_report <- function(inter_dir) {
    path <- file.path(inter_dir, "preflight.csv")
    if (!file.exists(path)) return(NULL)

    tryCatch(
        read.csv(path, stringsAsFactors = FALSE),
        error = function(e) NULL
    )
}


# =============================================================================
# Snakemake skip sentinel helpers
# =============================================================================
# When a stratum is skipped (too few samples, upstream failure, etc.),
# Snakemake still expects all declared output files to exist.
# These two functions handle that contract:
#
#   write_skip_sentinel(dir, filenames)
#       Called via on.exit() at the top of each script's main().
#       If any declared output file is missing when the script exits,
#       it writes a one-line placeholder so Snakemake can proceed.
#
#   is_skip_sentinel(path)
#       Called by downstream scripts before reading an upstream file.
#       Returns TRUE if the file is a placeholder, so the script can
#       bail out gracefully instead of crashing on malformed input.
# =============================================================================


#' Register an on.exit handler that writes skip sentinels for missing outputs
#'
#' Call this once at the top of main(), right after inter_dir / out_dir is
#' defined. If the script exits without writing one of the declared output
#' files (due to a skip or error), a one-line placeholder is written so
#' Snakemake's output contract is satisfied.
#'
#' @param output_dir  Directory where outputs are written
#' @param filenames   Character vector of expected output filenames
#' @param env         Environment to register on.exit in (default: parent frame)
register_skip_sentinel <- function(output_dir, filenames, env = parent.frame()) {
    do.call(on.exit, list(substitute({
        for (.sentinel_f in .sentinel_files) {
            .sentinel_p <- file.path(.sentinel_dir, .sentinel_f)
            if (!file.exists(.sentinel_p)) {
                dir.create(dirname(.sentinel_p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .sentinel_p)
            }
        }
    }, list(.sentinel_dir = output_dir, .sentinel_files = filenames)),
    add = TRUE), envir = env)
}


#' Check whether a file is a skip sentinel written by register_skip_sentinel
#'
#' @param path  Path to file to check
#' @return TRUE if the file is a sentinel placeholder, FALSE otherwise
is_skip_sentinel <- function(path) {
    if (!file.exists(path)) return(FALSE)
    tryCatch({
        first_line <- readLines(path, n = 1, warn = FALSE)
        identical(trimws(first_line), "skipped=TRUE")
    }, error = function(e) FALSE)
}


# =============================================================================
# Null-coalescing operator
# =============================================================================
# Defined once here; all scripts source io_utils.R so it's always available.
# Handles both NULL and length-0 (e.g. character(0) from failed lookups).
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a