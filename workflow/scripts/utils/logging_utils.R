# =============================================================================
# logging_utils.R
# Structured logging for EasyDE
#
# Every log line is machine-readable AND human-readable:
#   [2026-02-26 14:32:01] [INFO]  contrast=T1D_vs_ND | biosample=Beta | step=filter_samples | n_samples=47->31
#   [2026-02-26 14:32:03] [ERROR] contrast=T1D_vs_ND | biosample=Acinar | step=deseq_base | msg=every gene contains at least one zero
#
# Two log files are written per job:
#   {biosample}.log        — full structured log (all levels)
#   {biosample}.errors.log — ONLY created if WARN/ERROR/SKIP lines exist
#                            Absence of this file = clean run for that biosample
#
# Usage:
#   logger <- make_logger(contrast_id = "T1D_vs_ND", biosample = "Beta",
#                         log_file = "logs/T1D_vs_ND/Beta.log", level = "INFO")
#   logger$info("filter_samples", "n_samples=47->31")
#   logger$warn("donor_dedup", "4 samples removed")
#   logger$error("deseq_base", "every gene contains at least one zero")
#   logger$skip("deseq_base_failed", "insufficient samples after filtering")
# =============================================================================


# Log level hierarchy — messages below the configured level are suppressed
LOG_LEVELS <- c(DEBUG = 0, INFO = 1, WARN = 2, ERROR = 3, SKIP = 3)

#' Append a formatted error to a running error list
#'
#' Used by any validation or multi-step script that collects ALL errors
#' before stopping — so the user sees every problem at once, not one per run.
#'
#' @param errors  Character vector of accumulated error strings
#' @param msg     Human-readable description of the error
#' @return Updated character vector with the new error appended
add_error <- function(errors, msg) {
    c(errors, paste0("  [ERROR] ", msg))
}


#' Append a formatted warning to a running warning list
#'
#' Warnings are non-fatal — the pipeline can proceed but should inform the user.
#' Example: a latent_corr_var not found in metadata (will be skipped at runtime).
#'
#' @param warnings  Character vector of accumulated warning strings
#' @param msg       Human-readable description of the warning
#' @return Updated character vector with the new warning appended
add_warning <- function(warnings, msg) {
    c(warnings, paste0("  [WARN]  ", msg))
}




#' Create a logger instance bound to a specific contrast + biosample
#'
#' Returns a named list of logging functions (info, warn, error, skip, debug).
#' Each function writes a structured line to the full log file and stderr.
#' If any WARN/ERROR/SKIP lines are emitted, they are ALSO written to a
#' separate .errors.log file — which only exists when something went wrong.
#'
#' @param contrast_id  Contrast identifier string (e.g. "T1D_vs_ND")
#' @param biosample    Biosample identifier (e.g. "Beta", "islet_donor_42")
#' @param log_file     Path to the full log file. Created if missing.
#'                     The error log is written to the same dir as {basename}.errors.log
#' @param level        Minimum log level to emit: "DEBUG" | "INFO" | "WARN" | "ERROR"
#' @return Named list with functions: $debug, $info, $warn, $error, $skip
make_logger <- function(contrast_id, biosample, log_file, level = "INFO") {

    min_level  <- LOG_LEVELS[[toupper(level)]]
    if (is.null(min_level)) {
        stop(sprintf("Invalid log level: '%s'. Choose from: DEBUG, INFO, WARN, ERROR", level))
    }

    # Derive error-only log path from the full log path
    error_log_file <- sub("\\.log$", ".errors.log", log_file)

    # Ensure log directory exists
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

    # Core write function
    write_log <- function(log_level, step, message) {

        if (LOG_LEVELS[[log_level]] < min_level) return(invisible(NULL))

        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

        line <- sprintf(
            "[%s] [%-5s] contrast=%s | biosample=%s | step=%s | %s",
            timestamp,
            log_level,
            contrast_id,
            biosample,
            step,
            message
        )

        # Always write to full log
        cat(line, "\n", file = log_file, append = TRUE)

        # Write WARN/ERROR/SKIP to error-only log too
        if (log_level %in% c("WARN", "ERROR", "SKIP")) {
            cat(line, "\n", file = error_log_file, append = TRUE)
        }

        # Mirror WARN/ERROR/SKIP to stderr so SLURM/terminal shows problems in real time.
        # INFO and DEBUG stay file-only to keep the console quiet.
        if (log_level %in% c("WARN", "ERROR", "SKIP")) {
            message(line)
        }

        invisible(NULL)
    }

    list(
        debug = function(step, msg) write_log("DEBUG", step, msg),
        info  = function(step, msg) write_log("INFO",  step, msg),
        warn  = function(step, msg) write_log("WARN",  step, msg),
        error = function(step, msg) write_log("ERROR", step, msg),
        skip  = function(reason, msg) write_log("SKIP",
                    step    = "SKIPPED",
                    message = sprintf("reason=%s | %s", reason, msg))
    )
}


#' Summarize all log files into a single error/skip report
#'
#' Reads all .log files under logs_dir, extracts WARN/ERROR/SKIP lines,
#' and writes a consolidated report. Useful for post-run review.
#'
#' @param logs_dir   Root logs directory
#' @param output_file Path to write the summary report (TSV)
#' @return data.frame of flagged log lines (invisibly)
summarize_logs <- function(logs_dir, output_file) {

    log_files <- list.files(logs_dir, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)

    if (length(log_files) == 0) {
        message("No log files found in: ", logs_dir)
        return(invisible(data.frame()))
    }

    flagged_lines <- lapply(log_files, function(f) {
        lines <- readLines(f, warn = FALSE)
        flagged <- lines[grepl("\\[(WARN |ERROR|SKIP )\\]", lines)]
        if (length(flagged) == 0) return(NULL)
        data.frame(log_file = f, line = flagged, stringsAsFactors = FALSE)
    })

    flagged_df <- do.call(rbind, Filter(Negate(is.null), flagged_lines))

    if (is.null(flagged_df) || nrow(flagged_df) == 0) {
        message("No warnings, errors or skips found across all logs.")
        return(invisible(data.frame()))
    }

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write.table(flagged_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[log summary] %d flagged lines written to: %s", nrow(flagged_df), output_file))

    return(invisible(flagged_df))
}


#' Wrap a code block in tryCatch with automatic structured logging
#'
#' On success: logs an INFO message and returns the result.
#' On failure: logs an ERROR message and returns NULL (caller handles the skip).
#'
#' This replaces the repeated skip <- FALSE / tryCatch / if (skip) { next }
#' pattern that appears 6+ times in the original script.
#'
#' Example:
#'   dds <- try_logged(
#'       expr    = DESeq(dds),
#'       logger  = logger,
#'       step    = "deseq_base",
#'       success = "DESeq2 completed"
#'   )
#'   if (is.null(dds)) next   # one clean exit point per step
#'
#' @param expr     Expression to evaluate (use {} for multi-line blocks)
#' @param logger   Logger instance from make_logger()
#' @param step     Step name string (for log line)
#' @param success  Message to log on success (default: "completed")
#' @return Result of expr on success, NULL on failure
try_logged <- function(expr, logger, step, success = "completed") {

    result <- tryCatch(
        {
            out <- expr
            logger$info(step, success)
            out
        },
        error = function(e) {
            msg <- conditionMessage(e)
            if (is.null(msg) || !nzchar(trimws(msg))) {
                msg <- paste0("unknown error (class: ", paste(class(e), collapse = "/"), ")")
            }
            logger$error(step, sprintf("msg=%s", msg))
            NULL
        }
    )

    return(result)
}