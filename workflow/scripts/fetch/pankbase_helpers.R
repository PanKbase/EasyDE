# =============================================================================
# pankbase_helpers.R
# PanKbase-specific data fetching and formatting for EasyDE
#
# PURPOSE:
#   This is a DATA SOURCE ADAPTER — not core pipeline logic.
#   It fetches raw data from PanKbase and writes the two standard input files
#   that the EasyDE pipeline expects:
#
#     data/counts/{stratum_safe}.counts.csv   one per cell type / stratum
#     data/sample_metadata.csv               merged donor + biosample metadata
#
#   Once these files exist, the core pipeline is completely data-source agnostic.
#   Users with non-PanKbase data skip this script entirely and provide their
#   own files in the same format.
#
# USAGE:
#   Rscript workflow/scripts/fetch/pankbase_helpers.R \
#       --config config/config.yaml \
#       --pankbase-config workflow/scripts/fetch/pankbase_config.yaml \
#       --strata "Beta|Alpha|Delta"
#
# LOGGING:
#   Structured logs are written to logs/fetch/ via logging_utils.R.
#   Each stratum gets its own log file: logs/fetch/{stratum}.log
#   Warnings/errors also written to: logs/fetch/{stratum}.errors.log
#   A cross-stratum summary is written to: logs/fetch/fetch_warnings_summary.tsv
#
# DEPENDENCIES:
#   Pipeline settings (output paths, donor_col) are read from config/config.yaml.
#   PanKbase-specific settings (URLs, cache dir) are read from pankbase_config.yaml.
#   See that file for full documentation.
# =============================================================================

suppressMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
    library(stringr)
    library(yaml)
    library(R.utils)
})

# Resolve utils directory relative to this file.
# commandArgs() is the reliable method under Rscript; sys.frame() only works
# when the file is source()-d interactively.
.script_dir <- tryCatch({
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("--file=", "", file_arg[1])))
    } else {
        dirname(normalizePath(sys.frame(1)$ofile))
    }
}, error = function(e) getwd())

source(file.path(.script_dir, "..", "utils", "io_utils.R"))
source(file.path(.script_dir, "..", "utils", "logging_utils.R"))


# -----------------------------------------------------------------------------
# URL and name helpers
# -----------------------------------------------------------------------------

#' Convert a PanKbase cell type display name to its URL-encoded form
#'
#' PanKbase uses specific URL encodings for cell type names that differ
#' from simply removing spaces. This function centralizes all known
#' special cases so they don't scatter across the codebase.
#'
#' @param cell_type  Cell type display name (e.g. "MUC5B+ Ductal", "Beta")
#' @return URL-encoded string suitable for embedding in a PanKbase S3 URL
encode_celltype_for_url <- function(cell_type) {

    # Known special cases where the URL name differs from the display name
    url_name <- switch(cell_type,
        "MUC5B+ Ductal"        = "MUC5b+Ductal",
        "Immune (Macrophages)" = "Immune",
        gsub(" ", "", cell_type)    # default: remove spaces only
    )

    URLencode(url_name, reserved = TRUE)
}


#' Build the full PanKbase S3 URL for a cell type count matrix
#'
#' @param cell_type    Cell type display name
#' @param base_url     Base S3 URL (from config pankbase.base_counts_url)
#' @param file_pattern File suffix pattern (from config pankbase.file_pattern)
#' @return Full URL string
build_pankbase_counts_url <- function(cell_type, base_url, file_pattern) {
    encoded <- encode_celltype_for_url(cell_type)
    paste0(base_url, encoded, file_pattern)
}


#' Convert a cell type display name to a safe filename stem
#'
#' Uses underscores (not dots) for safety in R formulas and file systems.
#' This safe name is used for the counts file and as the stratum identifier
#' throughout the pipeline.
#'
#' @param cell_type  Cell type display name (e.g. "MUC5B+ Ductal")
#' @return Safe string (e.g. "MUC5B__Ductal")
celltype_to_safe_name <- function(cell_type) {
    gsub("[^[:alnum:]]", "_", cell_type)
}


# -----------------------------------------------------------------------------
# Metadata fetching
# -----------------------------------------------------------------------------

#' Fetch PanKbase donor and biosample metadata
#'
#' Downloads donor metadata (tar.gz) and biosample metadata (plain TSV) from
#' PanKbase S3, caches them locally to avoid re-downloading on subsequent runs,
#' and returns them as separate tables.
#'
#' The two tables are returned separately because the join key (donor_accession)
#' comes from the count matrix header rows, not from the biosample table.
#' The chained join happens per-stratum in fetch_pankbase_counts().
#'
#' @param donor_url     URL to donor metadata tar.gz
#' @param biosample_url URL to biosample metadata TSV
#' @param cache_dir     Local directory to cache downloaded files
#' @param logger        Logger instance from make_logger()
#' @return Named list: list(donor = data.frame, biosample = data.frame)
fetch_pankbase_metadata <- function(donor_url,
                                    biosample_url,
                                    cache_dir = "data/cache",
                                    logger    = NULL) {

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # Thin wrappers so we log to the logger if present, or fall back to message()
    log_info  <- function(step, msg) {
        if (!is.null(logger)) logger$info(step, msg)
        else message(sprintf("[metadata] %s | %s", step, msg))
    }
    log_warn  <- function(step, msg) {
        if (!is.null(logger)) logger$warn(step, msg)
        else message(sprintf("[metadata] WARN | %s | %s", step, msg))
    }
    log_error <- function(step, msg) {
        if (!is.null(logger)) logger$error(step, msg)
        else message(sprintf("[metadata] ERROR | %s | %s", step, msg))
    }

    # ---- Donor metadata: delivered as tar.gz, must extract before reading ----
    donor_tar   <- file.path(cache_dir, "pankbase_donor.tar.gz")
    donor_cache <- file.path(cache_dir, "pankbase_donor_metadata.tsv")

    if (!file.exists(donor_cache)) {

        log_info("donor_download", sprintf("url=%s", donor_url))

        withCallingHandlers(
            download.file(donor_url, destfile = donor_tar, quiet = TRUE, mode = "wb"),
            warning = function(w) {
                log_warn("donor_download", conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )

        archive_files <- untar(donor_tar, list = TRUE)
        txt_file      <- archive_files[grepl("\\.txt$", archive_files)]
        if (length(txt_file) == 0) {
            log_error("donor_extract", "No .txt file found inside donor metadata tar.gz")
            stop("No .txt file found inside donor metadata tar.gz")
        }

        log_info("donor_extract", sprintf("file=%s", txt_file[1]))
        untar(donor_tar, files = txt_file[1], exdir = cache_dir)

        extracted_path <- file.path(cache_dir, txt_file[1])
        file.copy(extracted_path, donor_cache, overwrite = TRUE)
        log_info("donor_cache_save", sprintf("path=%s", donor_cache))

    } else {
        log_info("donor_cache_hit", sprintf("path=%s", donor_cache))
    }

    donor_meta <- withCallingHandlers(
        fread(donor_cache) %>% as.data.frame(),
        warning = function(w) {
            log_warn("donor_read", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    donor_meta$donor_accession <- donor_meta$Accession

    # ---- Biosample metadata: plain TSV, fread handles directly ----
    biosample_cache <- file.path(cache_dir, "pankbase_biosample_metadata.tsv")

    log_info("biosample_fetch", sprintf("url=%s", biosample_url))

    biosample_meta <- withCallingHandlers(
        load_with_cache(biosample_url, cache_path = biosample_cache),
        warning = function(w) {
            log_warn("biosample_read", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    biosample_meta$biosample_accession <- biosample_meta$Accession

    log_info("metadata_ready",
             sprintf("donor=%d x %d | biosample=%d x %d",
                     nrow(donor_meta),     ncol(donor_meta),
                     nrow(biosample_meta), ncol(biosample_meta)))

    return(list(donor = donor_meta, biosample = biosample_meta))
}


# -----------------------------------------------------------------------------
# Count matrix fetching
# -----------------------------------------------------------------------------

#' Fetch one PanKbase cell type count matrix and extract its sample ID table
#'
#' PanKbase pseudobulk matrices have 4 metadata rows at the top before gene
#' counts begin. This function:
#'   1. Downloads the count matrix from PanKbase S3
#'   2. Extracts the 4-row header into a per-sample ID table
#'   3. Joins sample IDs -> donor_meta -> biosample_meta (original pipeline order)
#'   4. Strips the header rows to produce a clean genes x samples count matrix
#'   5. Writes the clean count matrix to data/counts/{safe_name}.counts.csv
#'
#' All warnings are captured via withCallingHandlers and routed through
#' logger$warn rather than being printed raw to stderr.
#'
#' @param cell_type       Cell type display name (e.g. "Beta", "MUC5B+ Ductal")
#' @param base_url        Base PanKbase S3 URL for count matrices
#' @param file_pattern    File suffix (e.g. "_sample_gex_total_counts.fixed.txt")
#' @param donor_meta      Donor metadata data.frame (from fetch_pankbase_metadata)
#' @param biosample_meta  Biosample metadata data.frame (from fetch_pankbase_metadata)
#' @param counts_dir      Output directory for clean count matrix CSVs
#' @param biosample_id_col  Column name to assign to the sample ID column
#' @param logger          Logger instance from make_logger()
#' @return data.frame of per-sample IDs with donor + biosample columns joined,
#'         or NULL on failure
fetch_pankbase_counts <- function(cell_type,
                                  base_url,
                                  file_pattern,
                                  donor_meta,
                                  biosample_meta,
                                  counts_dir       = "data/counts",
                                  biosample_id_col = "sample_accession",
                                  logger           = NULL) {

    safe_name   <- celltype_to_safe_name(cell_type)
    output_path <- file.path(counts_dir, paste0(safe_name, ".counts.csv"))

    log_info  <- function(step, msg) if (!is.null(logger)) logger$info(step, msg)
    log_warn  <- function(step, msg) if (!is.null(logger)) logger$warn(step, msg)
    log_error <- function(step, msg) if (!is.null(logger)) logger$error(step, msg)
    log_skip  <- function(reason, msg) if (!is.null(logger)) logger$skip(reason, msg)

    # ---- Download raw matrix ----
    url <- build_pankbase_counts_url(cell_type, base_url, file_pattern)
    log_info("counts_fetch", sprintf("url=%s", url))

    raw_mat <- withCallingHandlers(
        tryCatch(
            load_count_matrix(url),
            error = function(e) {
                log_error("counts_fetch", sprintf("msg=%s", conditionMessage(e)))
                NULL
            }
        ),
        warning = function(w) {
            log_warn("counts_fetch", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    if (is.null(raw_mat)) {
        log_skip("download_failed", sprintf("cell_type=%s | url=%s", cell_type, url))
        return(NULL)
    }

    log_info("counts_fetch", sprintf("raw_mat=%d rows x %d cols", nrow(raw_mat), ncol(raw_mat)))



    # ---- Extract 4-row PanKbase header into sample ID table ----
    # Row order: center_donor_id, donor_accession, biosample_accession, treatment
    sample_ids <- withCallingHandlers(
        {
            ids           <- t(raw_mat[1:4, ]) %>% as.data.frame()
            colnames(ids) <- ids[1, ]
            ids           <- ids[-1, ]
            ids <- tibble::rownames_to_column(ids, var = biosample_id_col)
            # Normalize sample IDs to match fread column name sanitization:
            # fread and as.data.frame both convert "-" to "." in column names,
            # so HPAP-019 becomes HPAP.019 in the count matrix.
            # We apply the same transformation to sample_accession so joins work.
            ids[[biosample_id_col]] <- gsub("-", ".", ids[[biosample_id_col]])
            ids
        },
        warning = function(w) {
            log_warn("header_parse", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    # ---- Chain joins: pankbase_ids -> donor -> biosample (matches original) ----
    sample_ids <- withCallingHandlers(
        {
            ids <- left_join(sample_ids, donor_meta,    by = "donor_accession")
            left_join(ids, biosample_meta, by = "biosample_accession")
        },
        warning = function(w) {
            log_warn("metadata_join", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    log_info("metadata_join",
             sprintf("n_samples=%d | n_meta_cols=%d", nrow(sample_ids), ncol(sample_ids)))

    # ---- Strip header rows to get clean genes x samples count matrix ----
    count_mat <- withCallingHandlers(
        {
            mat        <- as.data.frame(raw_mat[-c(1:4), ])
            gene_names <- mat$V1
            mat        <- mat[, colnames(mat) != "V1", drop = FALSE]

            # PanKbase uses "+" in sample column names — replace for R safety
            colnames(mat) <- gsub("\\+", "-", colnames(mat))
            rownames(mat) <- gene_names

            mat <- as.data.frame(lapply(mat, as.numeric))
            rownames(mat) <- gene_names
            mat
        },
        warning = function(w) {
            log_warn("count_matrix_clean", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    # ---- Write clean count matrix ----
    dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
    out_mat <- cbind(gene = rownames(count_mat), count_mat)

    withCallingHandlers(
        fwrite(out_mat, output_path),
        warning = function(w) {
            log_warn("counts_write", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    log_info("counts_write",
             sprintf("path=%s | genes=%d | samples=%d",
                     output_path, nrow(count_mat), ncol(count_mat)))

    return(sample_ids)
}


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

#' Fetch all PanKbase data for a given list of strata and write pipeline inputs
#'
#' Orchestrates the full PanKbase fetch:
#'   1. Fetches and caches donor + biosample metadata once
#'   2. For each stratum: fetches count matrix, extracts + joins sample IDs
#'   3. Writes one count matrix file per stratum to counts_dir
#'   4. Writes the merged sample_metadata.csv
#'   5. Writes structured logs to logs/fetch/ (one per stratum + warnings summary)
#'
#' After this function completes, data/ contains everything the core EasyDE
#' pipeline needs — no further PanKbase interaction is required.
#'
#' @param cfg      Parsed main config.yaml as a named list (pipeline settings)
#' @param pb_cfg   Parsed pankbase_config.yaml as a named list (source settings)
#' @param strata   Character vector of cell type / stratum names to fetch
#' @return Invisibly returns a data.frame summarizing fetch outcomes per stratum
fetch_all_pankbase_data <- function(cfg, pb_cfg, strata) {

    if (is.null(pb_cfg)) {
        stop("pankbase_config.yaml could not be parsed. ",
             "Check workflow/scripts/fetch/pankbase_config.yaml.")
    }

    logs_dir <- file.path(cfg$logs_dir %||% "logs", "fetch")
    dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

    # Step 1: metadata — one logger for the metadata step
    meta_logger <- make_logger(
        contrast_id = "fetch",
        biosample   = "metadata",
        log_file    = file.path(logs_dir, "metadata.log")
    )

    meta <- tryCatch(
        fetch_pankbase_metadata(
            donor_url     = pb_cfg$donor_metadata_url,
            biosample_url = pb_cfg$biosample_metadata_url,
            cache_dir     = pb_cfg$cache_dir,
            logger        = meta_logger
        ),
        error = function(e) {
            meta_logger$error("fetch_metadata", sprintf("msg=%s", conditionMessage(e)))
            stop(sprintf("Metadata fetch failed: %s\nCheck: %s",
                         conditionMessage(e),
                         file.path(logs_dir, "metadata.log")))
        }
    )

    # Step 2: count matrices — one logger per stratum
    fetch_summary  <- data.frame()
    all_sample_ids <- list()

    for (stratum in strata) {

        safe_name <- celltype_to_safe_name(stratum)
        logger    <- make_logger(
            contrast_id = "fetch",
            biosample   = safe_name,
            log_file    = file.path(logs_dir, paste0(safe_name, ".log"))
        )

        message(sprintf("[fetch] %-25s ...", stratum))
        logger$info("start", sprintf("cell_type=%s", stratum))

        sample_ids <- fetch_pankbase_counts(
            cell_type        = stratum,
            base_url         = pb_cfg$base_counts_url,
            file_pattern     = pb_cfg$file_pattern,
            donor_meta       = meta$donor,
            biosample_meta   = meta$biosample,
            counts_dir       = cfg$inputs$counts_dir,
            biosample_id_col = "sample_accession",
            logger           = logger
        )

        status <- if (is.null(sample_ids)) "failed" else "success"
        logger$info("done", sprintf("status=%s", status))
        message(sprintf("[fetch] %-25s done (%s)", stratum,
                        if (is.null(sample_ids)) "FAILED" else sprintf("%d samples", nrow(sample_ids))))

        fetch_summary <- rbind(fetch_summary, data.frame(
            stratum   = stratum,
            status    = status,
            n_samples = if (is.null(sample_ids)) NA_integer_ else nrow(sample_ids),
            log_file  = file.path(logs_dir, paste0(safe_name, ".log")),
            stringsAsFactors = FALSE
        ))

        if (!is.null(sample_ids)) {
            sample_ids$stratum <- stratum
            all_sample_ids[[stratum]] <- sample_ids
        }
    }

    # Step 3: write merged sample metadata
    if (length(all_sample_ids) > 0) {

        all_ids_df <- do.call(rbind, all_sample_ids)

        # Join supplemental local metadata files defined in pankbase_config.yaml.
        # These are PanKbase-project-specific resources (chemistry, cell counts)
        # that live in resources/ and are committed to the repository.
        if (!is.null(pb_cfg$supplemental_metadata)) {
            for (supp in pb_cfg$supplemental_metadata) {
                supp_path <- supp$path
                join_col  <- supp$join_col
                if (!file.exists(supp_path)) {
                    message(sprintf("[pankbase] [WARN] supplemental file not found, skipping: %s", supp_path))
                    next
                }
                supp_df <- fread(supp_path) %>% as.data.frame()
                if (!join_col %in% colnames(supp_df)) {
                    message(sprintf("[pankbase] [WARN] join_col '%s' not found in %s, skipping", join_col, supp_path))
                    next
                }
                n_before <- ncol(all_ids_df)
                colnames(supp_df)[colnames(supp_df) == join_col] <- "sample_accession"
                # Normalize to dots to match count matrix column name convention
                supp_df[["sample_accession"]] <- gsub("-", ".", supp_df[["sample_accession"]])
                all_ids_df <- left_join(all_ids_df, supp_df, by = "sample_accession")
                message(sprintf("[pankbase] [join] %-35s +%d column(s) -> metadata now %d cols",
                                basename(supp_path),
                                ncol(all_ids_df) - n_before,
                                ncol(all_ids_df)))
            }
        }

        metadata_path <- cfg$inputs$sample_metadata
        dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
        fwrite(all_ids_df, metadata_path)
        message(sprintf("[pankbase] [write] sample_metadata: %s  (%d rows x %d cols)",
                        metadata_path, nrow(all_ids_df), ncol(all_ids_df)))
    }

    # Step 4: print fetch summary to console
    message("\n[pankbase] Fetch summary:")
    print(fetch_summary[, c("stratum", "status", "n_samples")])

    n_ok   <- sum(fetch_summary$status == "success")
    n_fail <- sum(fetch_summary$status == "failed")
    message(sprintf("[pankbase] %d/%d strata fetched successfully, %d failed.",
                    n_ok, nrow(fetch_summary), n_fail))

    if (n_fail > 0) {
        failed <- fetch_summary$stratum[fetch_summary$status == "failed"]
        message(sprintf("[pankbase] Failed strata: %s", paste(failed, collapse = ", ")))
        message(sprintf("[pankbase] Check logs in: %s", logs_dir))
    }

    # Step 5: write cross-stratum warnings/errors summary
    summarize_logs(
        logs_dir    = logs_dir,
        output_file = file.path(logs_dir, "fetch_warnings_summary.tsv")
    )

    return(invisible(fetch_summary))
}


# -----------------------------------------------------------------------------
# CLI entry point — only runs when called as Rscript, not when sourced
# -----------------------------------------------------------------------------

if (!interactive() && identical(sys.nframe(), 0L)) {

    library(optparse)

    opt_list <- list(
        make_option("--config", type = "character",
                    default = "config/config.yaml",
                    help = "Path to main pipeline config.yaml [default: config/config.yaml]"),
        make_option("--pankbase-config", type = "character",
                    default = "workflow/scripts/fetch/pankbase_config.yaml",
                    help = "Path to pankbase_config.yaml [default: workflow/scripts/fetch/pankbase_config.yaml]"),
        make_option("--strata", type = "character",
                    help = "Pipe-separated list of strata to fetch (e.g. 'Beta|Alpha|Delta')")
    )

    opts <- parse_args(OptionParser(option_list = opt_list))

    if (is.null(opts$strata)) stop("--strata is required")

    cfg    <- suppressWarnings(yaml::read_yaml(opts$config))
    pb_cfg <- suppressWarnings(yaml::read_yaml(opts$`pankbase-config`))
    strata <- trimws(strsplit(opts$strata, "\\|")[[1]])

    fetch_all_pankbase_data(cfg, pb_cfg, strata)
}