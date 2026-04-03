# =============================================================================
# 09_pipeline_summary.R
# EasyDE -- Pipeline-Level Summary
#
# PURPOSE:
#   Runs ONCE after all contrasts complete. Reads all per-contrast
#   contrast_summary.csv files and generates a pipeline-level status
#   heatmap (contrasts x strata).
#
# OUTPUTS:
#   {results_dir}/pipeline_status_overview.pdf
#   {results_dir}/pipeline_summary.csv
#
# USAGE:
#   Rscript workflow/scripts/09_pipeline_summary.R --config config/config.yaml
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(optparse)
    library(ggplot2)
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
source(file.path(.script_dir, "utils", "stats_utils.R"))
source(file.path(.script_dir, "utils", "plot_utils.R"))


# =============================================================================
# Main
# =============================================================================

main <- function(config_path) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df
    results_dir  <- cfg$outputs$results_dir

    contrast_ids <- unique(contrasts_df$contrast_id)
    cat(sprintf("  Pipeline summary: %d contrasts\n", length(contrast_ids)))

    # --- Collect all contrast summaries ---
    all_rows <- list()
    for (cid in contrast_ids) {
        summary_path <- file.path(results_dir, cid, "summary", "contrast_summary.csv")
        if (file.exists(summary_path)) {
            df <- tryCatch(fread(summary_path) %>% as.data.frame(), error = function(e) NULL)
            if (!is.null(df) && nrow(df) > 0) {
                all_rows[[cid]] <- df
            }
        } else {
            cat(sprintf("  [WARN] Missing summary for %s\n", cid))
        }
    }

    if (length(all_rows) == 0) {
        cat("  [WARN] No contrast summaries found — skipping pipeline summary\n")
        return(invisible(NULL))
    }

    all_summaries <- do.call(rbind, all_rows)
    cat(sprintf("  Collected %d rows across %d contrasts\n",
                nrow(all_summaries), length(all_rows)))

    # --- Write combined CSV ---
    write_result(all_summaries, file.path(results_dir, "pipeline_summary.csv"), sep = ",")

    # --- Generate pipeline-level status heatmap ---
    tryCatch({
        p <- plot_pipeline_status_heatmap(all_summaries)
        if (!is.null(p)) {
            n_contrasts <- length(unique(all_summaries$contrast_id))
            n_strata    <- length(unique(all_summaries$biosample))
            pdf_w <- max(10, n_strata * 0.8 + 3)
            pdf_h <- max(6,  n_contrasts * 0.4 + 3)

            pdf(file.path(results_dir, "pipeline_status_overview.pdf"),
                width = pdf_w, height = pdf_h)
            print(p)
            dev.off()
            cat(sprintf("  Pipeline status heatmap: %s\n",
                        file.path(results_dir, "pipeline_status_overview.pdf")))
        }
    }, error = function(e) {
        cat(sprintf("  [WARN] Could not generate pipeline status heatmap: %s\n", e$message))
    })
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {
    opt_list <- list(
        make_option("--config", type = "character", default = "config/config.yaml",
                    help = "Path to pipeline config YAML")
    )
    opts <- parse_args(OptionParser(option_list = opt_list))
    main(opts$config)
}
