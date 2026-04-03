# =============================================================================
# plot_utils.R
# All plotting functions for EasyDE
#
# Every function returns a ggplot object - callers decide where to send it
# (pdf, png, print to screen). No dev.open/dev.off here.
# =============================================================================

suppressMessages({
    library(ggplot2)
    library(ggrepel)
    library(ggcorrplot)
    library(viridis)
    library(reshape2)
    library(dplyr)
    library(ggh4x)
    library(scales)
})

# Resolve utils directory relative to this file's location, robustly
# stats_utils.R is a sibling - sourced via .script_dir set by the calling script
source(file.path(.script_dir, "utils", "stats_utils.R"))


#' Volcano plot from DESeq2 results
#'
#' @param deseq_df   data.frame with columns: log2FoldChange, pvalue, padj, gene
#' @param fdr        FDR threshold for coloring significant genes (default 0.05)
#' @param title      Plot title (e.g. the formula used)
#' @param subtitle   Plot subtitle (e.g. "number DEGs: 42")
#' @param n_label    Number of top significant genes to label (default 20)
#' @return ggplot object
plot_volcano <- function(deseq_df, fdr = 0.05, title = "", subtitle = "", n_label = 20) {

    # Separate significant genes for labeling
    sig_df <- deseq_df[!is.na(deseq_df$padj) & deseq_df$padj < fdr, ]
    top_df <- sig_df[order(sig_df$pvalue), ]
    top_df <- head(top_df, n_label)

    plot_df         <- deseq_df[!is.na(deseq_df$padj), ]
    plot_df$signif  <- ifelse(plot_df$padj < fdr, "Signif.", "N.S.")

    p <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue))) +
        geom_point(aes(col = signif), size = 0.5) +
        scale_color_manual(values = c("N.S." = "gray", "Signif." = "firebrick")) +
        labs(col = "", title = title, subtitle = subtitle) +
        theme_bw() +
        theme(plot.title = element_text(size = 8))

    if (nrow(top_df) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_df,
            aes(x = log2FoldChange, y = -log10(pvalue), label = gene),
            size = 3
        )
    }

    return(p)
}


#' Correlation heatmap: RUV latent vars vs known covariates
#'
#' Displays t-statistics as cell values with BH-adjusted significance stars.
#' W factors correlated with the contrast variable (raw p < 0.05) are excluded
#' from the RUVseq model and marked with "(excl.)" on the y-axis.
#' Covariate correlations are shown for information only — they do not drive
#' exclusion since covariates are already in the DESeq2 design formula.
#'
#' @param t_mat          Numeric matrix of t-statistics (covariates x latent vars)
#' @param p_mat_display  Numeric matrix of BH-adjusted p-values (per-W, for stars)
#' @param disease_ws     Character vector of excluded W names (default NULL)
#' @param contrast_var   Sanitized contrast variable name (for caption, default NULL)
#' @param title          Plot title string
#' @param t_limits       Symmetric limits for color scale (default c(-15, 15))
#' @return ggplot object
plot_latent_correlation <- function(t_mat, p_mat_display,
                                     disease_ws   = NULL,
                                     contrast_var = NULL,
                                     title = "", t_limits = c(-15, 15)) {

    # Get significance asterisks from the display p-value matrix
    star_mat <- apply(p_mat_display, c(1, 2), pvalue_to_stars)

    # Rename excluded W columns so the y-axis labels show "(excl.)"
    if (!is.null(disease_ws) && length(disease_ws) > 0) {
        rename_w <- function(w) ifelse(w %in% disease_ws, paste0(w, " (excl.)"), w)
        colnames(t_mat)         <- rename_w(colnames(t_mat))
        colnames(p_mat_display) <- rename_w(colnames(p_mat_display))
        colnames(star_mat)      <- rename_w(colnames(star_mat))
    }

    cor_plot <- ggcorrplot(t_mat,
                           hc.order = FALSE,
                           lab      = TRUE,
                           title    = title,
                           lab_size = 2.25) +
        scale_fill_gradient2(
            low      = "#3A86FF",
            mid      = "white",
            high     = "#FD2244",
            limits   = t_limits,
            midpoint = 0,
            oob      = scales::squish
        )

    contrast_label <- if (!is.null(contrast_var)) contrast_var else "contrast var"
    caption_text <- paste0(
        "* padj<0.05  ** padj<0.01  *** padj<0.001 (BH-adjusted per-W)\n",
        "W exclusion: ", contrast_label, " at raw p<0.05 | covariate correlations shown for info only"
    )
    cor_plot <- cor_plot + labs(fill = "t-statistic", caption = caption_text)

    # Add significance stars — use ggcorrplot's internal data coordinates
    # directly to guarantee alignment (avoids Var1/Var2 mismatch issues)
    plot_data <- cor_plot[["data"]]
    plot_data$stars <- mapply(function(v1, v2) {
        # ggcorrplot maps: Var1 = rownames (covariates), Var2 = colnames (W factors)
        if (as.character(v1) %in% rownames(star_mat) &&
            as.character(v2) %in% colnames(star_mat)) {
            star_mat[as.character(v1), as.character(v2)]
        } else if (as.character(v2) %in% rownames(star_mat) &&
                   as.character(v1) %in% colnames(star_mat)) {
            # ggcorrplot may swap axes for non-square matrices
            star_mat[as.character(v2), as.character(v1)]
        } else {
            ""
        }
    }, plot_data$Var1, plot_data$Var2)
    star_data <- plot_data[plot_data$stars != "" & !is.na(plot_data$stars), ]

    if (nrow(star_data) > 0) {
        cor_plot <- cor_plot +
            geom_text(
                data    = star_data,
                aes(x = Var1, y = Var2, label = stars),
                nudge_y = 0.25,
                size    = 5,
                inherit.aes = FALSE
            )
    }
    cor_plot <- cor_plot +
        guides(fill = guide_legend(title = "t-statistic"))

    return(cor_plot)
}


#' RLE elbow plot: F-statistic vs number of RUV factors (k)
#'
#' Visualizes the ANOVA-based k selection. The chosen best_k is annotated.
#'
#' @param anova_df  data.frame with columns: k (factor), k_num, f_statistic, residual_variance
#' @param best_k    Integer: the selected best k value
#' @param title     Plot title string
#' @return ggplot object
plot_rle_elbow <- function(anova_df, best_k, title = "") {

    p <- ggplot(anova_df, aes(x = k, y = f_statistic, group = 1)) +
        geom_line() +
        geom_point(
            aes(size = residual_variance, fill = residual_variance),
            shape = 21, color = "black"
        ) +
        scale_fill_viridis(name = "Residual variance", option = "D") +
        scale_size_continuous(name = "Residual variance") +
        theme_classic() +
        labs(
            x     = "Number of RUV factors (k)",
            y     = "F-statistic\n(ANOVA: RLE ~ sample)",
            title = paste0(title, "\nBest-k -> ", best_k)
        ) +
        theme(
            axis.title  = element_text(size = 14),
            axis.text   = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )

    return(p)
}


#' fGSEA barplot: top N significant pathways by NES
#'
#' @param fgsea_sig_df  data.frame of significant fGSEA results (already filtered by FDR)
#'                      Must have columns: pathway, NES, padj
#' @param top_n         Max pathways to show (default 20, split top/bottom by NES)
#' @param title         Plot title
#' @return ggplot object, or NULL if no significant pathways
plot_fgsea_barplot <- function(fgsea_sig_df, top_n = 20, title = "Top significant pathways") {

    if (is.null(fgsea_sig_df) || nrow(fgsea_sig_df) == 0) {
        return(NULL)
    }

    n_show   <- min(top_n, nrow(fgsea_sig_df))
    plot_df  <- fgsea_sig_df[order(fgsea_sig_df$NES, decreasing = TRUE), ]
    plot_df  <- head(plot_df, n_show)

    p <- ggplot(plot_df, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = padj)) +
        coord_flip() +
        scale_fill_viridis(option = "C", direction = -1, name = "FDR") +
        labs(x = "Pathway", y = "Normalized Enrichment Score", title = title) +
        theme_minimal() +
        theme(text = element_text(size = 6))

    return(p)
}


# =============================================================================
# Status heatmap shared constants
# =============================================================================

# Canonical ordering and color palette for pipeline statuses
.status_levels <- c("success", "success_base_only",
                    "skipped_no_samples", "skipped_min_group",
                    "skipped_no_pairs", "skipped_preflight", "skipped",
                    "error", "not_run")

.status_colors <- c(
    "success"            = "#2ca02c",   # green
    "success_base_only"  = "#98df8a",   # light green
    "skipped_no_samples" = "#ff7f0e",   # orange
    "skipped_min_group"  = "#ffbb78",   # light orange
    "skipped_no_pairs"   = "#f0b27a",   # peach
    "skipped_preflight"  = "#e59866",   # tan
    "skipped"            = "#fad7a0",   # pale orange
    "error"              = "#d62728",   # red
    "not_run"            = "grey80"     # grey
)

#' Add DEG display label and factor-order status for heatmap data
#'
#' @param df  data.frame with columns: status, n_degs_base, n_degs_ruv
#' @return df with added columns: n_degs_display, deg_label, status (factored)
.prepare_status_df <- function(df) {
    df$n_degs_display <- ifelse(
        !is.na(df$n_degs_ruv), df$n_degs_ruv,
        ifelse(!is.na(df$n_degs_base), df$n_degs_base, NA)
    )
    df$deg_label <- ifelse(is.na(df$n_degs_display), "",
                           as.character(df$n_degs_display))
    df$status <- factor(df$status,
                         levels = intersect(.status_levels, unique(df$status)))
    df
}


#' Pipeline status overview heatmap
#'
#' Tile plot showing the status of each stratum (cell type) across the contrast.
#' Success statuses are shown in green shades, skips in orange shades,
#' errors in red, and not_run in grey. Cell labels show DEG counts where
#' available (RUV if present, else base).
#'
#' @param summary_df     data.frame: contrast_summary with columns biosample,
#'                       status, n_degs_base, n_degs_ruv
#' @param contrast_id    Character: contrast name for the title
#' @return ggplot object
plot_status_heatmap <- function(summary_df, contrast_id = "") {

    if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

    summary_df <- .prepare_status_df(summary_df)
    summary_df$biosample <- factor(summary_df$biosample,
                                    levels = sort(unique(summary_df$biosample)))
    used_colors <- .status_colors[names(.status_colors) %in% levels(summary_df$status)]

    p <- ggplot(summary_df, aes(x = biosample, y = 1, fill = status)) +
        geom_tile(color = "white", linewidth = 0.8) +
        geom_text(aes(label = deg_label), size = 2.5, color = "black") +
        scale_fill_manual(values = used_colors, name = "Status", drop = FALSE) +
        labs(
            title    = paste0("Pipeline Status: ", contrast_id),
            subtitle = sprintf("%d strata | %d success | %d skipped | %d errors",
                               nrow(summary_df),
                               sum(grepl("^success", summary_df$status)),
                               sum(grepl("^skipped", summary_df$status)),
                               sum(summary_df$status == "error")),
            x = "Cell type",
            y = ""
        ) +
        theme_minimal() +
        theme(
            axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            panel.grid       = element_blank(),
            legend.position  = "bottom",
            legend.text      = element_text(size = 7),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 9, color = "grey40")
        ) +
        guides(fill = guide_legend(nrow = 2))

    return(p)
}


# =============================================================================
# Pipeline-level status heatmap (all contrasts x strata)
# =============================================================================

#' Pipeline-level status heatmap across all contrasts and strata
#'
#' 2D heatmap: contrasts on y-axis, cell types on x-axis, fill = status,
#' DEG count labels in cells.
#'
#' @param all_summaries  data.frame with columns: contrast_id, biosample, status,
#'                       n_degs_base, n_degs_ruv (rbind of all contrast_summary.csv)
#' @param title          Plot title
#' @return ggplot object
plot_pipeline_status_heatmap <- function(all_summaries, title = "Pipeline Status Overview") {

    if (is.null(all_summaries) || nrow(all_summaries) == 0) return(NULL)

    all_summaries <- .prepare_status_df(all_summaries)
    all_summaries$biosample <- factor(all_summaries$biosample,
                                       levels = sort(unique(all_summaries$biosample)))
    all_summaries$contrast  <- factor(all_summaries$contrast,
                                       levels = rev(sort(unique(all_summaries$contrast))))
    used_colors <- .status_colors[names(.status_colors) %in% levels(all_summaries$status)]

    n_contrasts <- length(unique(all_summaries$contrast))
    n_strata    <- length(unique(all_summaries$biosample))
    n_success   <- sum(grepl("^success", all_summaries$status))
    n_skipped   <- sum(grepl("^skipped", all_summaries$status))
    n_errors    <- sum(all_summaries$status == "error")

    p <- ggplot(all_summaries, aes(x = biosample, y = contrast, fill = status)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = deg_label), size = 2, color = "black") +
        scale_fill_manual(values = used_colors, name = "Status", drop = FALSE) +
        labs(
            title    = title,
            subtitle = sprintf("%d contrasts x %d strata | %d success | %d skipped | %d errors",
                               n_contrasts, n_strata, n_success, n_skipped, n_errors),
            x = "Cell type",
            y = "Contrast"
        ) +
        theme_minimal() +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
            axis.text.y      = element_text(size = 7),
            panel.grid       = element_blank(),
            legend.position  = "bottom",
            legend.text      = element_text(size = 7),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 9, color = "grey40")
        ) +
        guides(fill = guide_legend(nrow = 2))

    return(p)
}


# =============================================================================
# Benchmark signature plots
# =============================================================================

#' Benchmark heatmap: fGSEA NES + FORA significance overlay
#'
#' For base and RUV panels, the fill is fGSEA NES. For the ruv_only panel,
#' the fill is the max FORA overlap fraction (up or down).
#' FORA significance is overlaid as directional markers:
#'   gold triangle-up = significant enrichment in upregulated DEGs
#'   cyan triangle-down = significant enrichment in downregulated DEGs
#'
#' @param bench_df   data.frame: consolidated benchmark results for one contrast
#' @param gene_set_name  "base", "ruv", or "ruv_only"
#' @param control_type_filter  "positive" or "negative" (new schema); used when
#'        control_type column exists. Takes precedence over categories param.
#' @param categories Character vector of class values to include (new schema)
#'        or old-style subcategory values (legacy fallback)
#' @param fora_fdr   FDR threshold for FORA significance markers (default 0.05)
#' @param title      Plot title
#' @return ggplot object, or NULL if no data
#' Save a heatmap plot result to PDF using its recommended dimensions
#'
#' Used by scripts 07 and 08 for benchmark heatmaps. The plot function returns
#' a list with $plot, $width, $height; this function saves it to a PDF.
#'
#' @param result   List with $plot (ggplot), $width, $height — as returned by
#'                 plot_fora_heatmap() or plot_benchmark_heatmap()
#' @param filepath Full path for the output PDF
#' @param logger   Optional logger instance for info messages
save_heatmap <- function(result, filepath, logger = NULL) {
    if (is.null(result)) return(invisible(NULL))
    pdf(filepath, width = result$width, height = result$height)
    print(result$plot)
    dev.off()
    if (!is.null(logger)) {
        logger$info("plot", basename(filepath))
    } else {
        cat(sprintf("  Benchmark: %s\n", filepath))
    }
}


plot_benchmark_heatmap <- function(bench_df, gene_set_name,
                                    control_type_filter = NULL,
                                    categories = NULL,
                                    fora_fdr = 0.05,
                                    title = NULL) {

    if (is.null(bench_df) || nrow(bench_df) == 0) return(NULL)

    has_class <- "class" %in% colnames(bench_df)

    # Filter to requested gene set
    df <- bench_df[bench_df$gene_set == gene_set_name, ]
    if (nrow(df) == 0) return(NULL)

    # Filter by control_type (new schema) or categories (legacy)
    if (!is.null(control_type_filter) && "control_type" %in% colnames(df)) {
        df <- df[df$control_type == control_type_filter, ]
    }
    if (!is.null(categories)) {
        if (has_class) {
            df <- df[df$class %in% categories, ]
        } else {
            df <- df[df$category %in% categories, ]
        }
    }
    if (nrow(df) == 0) return(NULL)

    # Determine facet column: use class (new) or category (legacy)
    facet_col <- if (has_class) "class" else "category"

    # --- Build fill matrix from fGSEA NES (base/ruv) or FORA fraction (ruv_only) ---
    if (gene_set_name == "ruv_only") {
        fora_df <- df[df$test == "fora", ]
        if (nrow(fora_df) == 0) return(NULL)

        fill_df <- fora_df %>%
            group_by(stratum, signature) %>%
            slice_min(pvalue, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            mutate(fill_value = ifelse(direction == "down",
                                        -statistic, statistic))
        fill_name <- "Overlap\nfraction"
    } else {
        fgsea_df <- df[df$test == "fgsea", ]
        if (nrow(fgsea_df) == 0) return(NULL)

        fill_df <- fgsea_df %>%
            mutate(fill_value = NES)
        fill_name <- "NES"
    }

    # Get FORA results for star overlay
    fora_sig <- df[df$test == "fora" &
                    !is.na(df$padj) & df$padj < fora_fdr, ]

    # --- Build tile data ---
    select_cols <- c("stratum", "signature", facet_col, "fill_value")
    if (has_class) select_cols <- c(select_cols, "category")
    tile_data <- fill_df %>% select(all_of(select_cols))

    if (nrow(tile_data) == 0) return(NULL)

    # Order signatures by facet column, then alphabetically
    sig_order <- tile_data %>%
        distinct(signature, .data[[facet_col]]) %>%
        arrange(.data[[facet_col]], signature)
    tile_data$signature <- factor(tile_data$signature, levels = sig_order$signature)

    # Default title
    if (is.null(title)) {
        label <- switch(gene_set_name,
                         base = "Base model", ruv = "RUV model",
                         ruv_only = "RUV-only DEGs")
        title <- paste0("Signature enrichment \u2014 ", label)
    }

    # --- Build plot ---
    max_abs <- max(abs(tile_data$fill_value), na.rm = TRUE)
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1

    p <- ggplot(tile_data, aes(x = stratum, y = signature)) +
        geom_tile(aes(fill = fill_value), color = "grey80", linewidth = 0.3) +
        scale_fill_gradient2(
            low      = "#3A86FF",
            mid      = "white",
            high     = "#FD2244",
            midpoint = 0,
            limits   = c(-max_abs, max_abs),
            name     = fill_name,
            oob      = scales::squish
        ) +
        labs(title = title, x = "", y = "") +
        theme_minimal() +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y      = element_text(size = 5.5),
            panel.grid       = element_blank(),
            plot.title       = element_text(size = 10, face = "bold"),
            legend.key.width = unit(0.4, "cm")
        )

    # Faceting: nested (category + class) when new schema, else flat
    n_facet_levels <- length(unique(sig_order[[facet_col]]))
    if (has_class && length(unique(tile_data$category)) > 1) {
        # Nested facets: category (outer) → class (inner)
        category_labels <- c(
            cell_identity = "Cell identity",
            artifact = "Artifact",
            subject_confounder = "Subject confounder"
        )
        tile_data$category_label <- factor(
            ifelse(tile_data$category %in% names(category_labels),
                   category_labels[tile_data$category],
                   tile_data$category),
            levels = unname(category_labels[
                intersect(names(category_labels), unique(tile_data$category))
            ])
        )
        tile_data$class_facet <- tile_data[[facet_col]]

        # Build strip colors
        cat_colors <- list(
            "Cell identity"      = c(outer = "#4CAF50", inner = "#C8E6C9"),
            "Artifact"           = c(outer = "#E05050", inner = "#FFCDD2"),
            "Subject confounder" = c(outer = "#F5A623", inner = "#FFE0B2")
        )

        outer_levels <- levels(tile_data$category_label)
        strip_bg_list <- list()
        for (ol in outer_levels) {
            col <- cat_colors[[ol]]
            if (!is.null(col)) {
                strip_bg_list <- c(strip_bg_list,
                    list(element_rect(fill = col[["outer"]])))
            } else {
                strip_bg_list <- c(strip_bg_list,
                    list(element_rect(fill = "grey90")))
            }
        }
        for (ol in outer_levels) {
            inner_classes <- tile_data %>%
                filter(category_label == ol) %>%
                pull(class_facet) %>% unique() %>% sort()
            col <- cat_colors[[ol]]
            for (ic in inner_classes) {
                if (!is.null(col)) {
                    strip_bg_list <- c(strip_bg_list,
                        list(element_rect(fill = col[["inner"]])))
                } else {
                    strip_bg_list <- c(strip_bg_list,
                        list(element_rect(fill = "grey95")))
                }
            }
        }

        themed_strip <- ggh4x::strip_nested(background_y = strip_bg_list)
        p <- p + ggh4x::facet_nested(
            category_label + class_facet ~ .,
            scales = "free_y", space = "free_y",
            switch = "y",
            nest_line = element_line(linewidth = 0.3),
            strip = themed_strip
        ) + theme(
            strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
            strip.placement = "outside"
        )
    } else if (n_facet_levels > 1) {
        tile_data$cat_facet <- factor(tile_data[[facet_col]],
                                       levels = unique(sig_order[[facet_col]]))
        p <- p + facet_grid(cat_facet ~ ., scales = "free_y", space = "free_y") +
            theme(strip.text.y = element_text(angle = 0, size = 7, face = "bold"))
    }

    # --- Overlay FORA significance markers ---
    if (nrow(fora_sig) > 0) {
        fora_sig$signature <- factor(fora_sig$signature,
                                      levels = levels(tile_data$signature))
        if (has_class && "category_label" %in% colnames(tile_data)) {
            fora_sig$category_label <- factor(
                ifelse(fora_sig$category %in% names(category_labels),
                       category_labels[fora_sig$category],
                       fora_sig$category),
                levels = levels(tile_data$category_label))
            fora_sig$class_facet <- fora_sig[[facet_col]]
        } else if (n_facet_levels > 1) {
            fora_sig$cat_facet <- factor(fora_sig[[facet_col]],
                                          levels = unique(sig_order[[facet_col]]))
        }

        fora_up   <- fora_sig[fora_sig$direction == "up", ]
        fora_down <- fora_sig[fora_sig$direction == "down", ]

        if (nrow(fora_up) > 0) {
            p <- p + geom_point(
                data = fora_up,
                aes(x = stratum, y = signature),
                shape = 24, size = 2, fill = "#FFD700", color = "black",
                stroke = 0.4, inherit.aes = FALSE
            )
        }
        if (nrow(fora_down) > 0) {
            p <- p + geom_point(
                data = fora_down,
                aes(x = stratum, y = signature),
                shape = 25, size = 2, fill = "#00CED1", color = "black",
                stroke = 0.4, inherit.aes = FALSE
            )
        }
    }

    p <- p + labs(caption = paste0(
        "\u25b2 gold = FORA enriched in upregulated DEGs (padj<",
        fora_fdr, ")   ",
        "\u25bc cyan = FORA enriched in downregulated DEGs"
    ))

    return(p)
}


#' Benchmark scatter: Base NES vs RUV NES per signature
#'
#' @param bench_df    Consolidated benchmark data
#' @param control_type_filter  "positive" or "negative" (new schema)
#' @param categories  Class values to include (new schema) or old category values
#' @param fgsea_fdr   FDR threshold for coloring (default 0.05)
#' @param title       Plot title
#' @return ggplot object, or NULL
plot_benchmark_scatter <- function(bench_df,
                                    control_type_filter = NULL,
                                    categories = NULL,
                                    fgsea_fdr = 0.05,
                                    title = "Base vs RUV enrichment") {

    if (is.null(bench_df) || nrow(bench_df) == 0) return(NULL)

    has_class <- "class" %in% colnames(bench_df)

    fgsea_df <- bench_df[bench_df$test == "fgsea", ]
    if (nrow(fgsea_df) == 0) return(NULL)

    # Filter by control_type (new schema) or categories (legacy)
    if (!is.null(control_type_filter) && "control_type" %in% colnames(fgsea_df)) {
        fgsea_df <- fgsea_df[fgsea_df$control_type == control_type_filter, ]
    }
    if (!is.null(categories)) {
        if (has_class) {
            fgsea_df <- fgsea_df[fgsea_df$class %in% categories, ]
        } else {
            fgsea_df <- fgsea_df[fgsea_df$category %in% categories, ]
        }
    }

    # Use class for shape aesthetic (new) or category (legacy)
    shape_col <- if (has_class) "class" else "category"

    # Pivot: one row per signature x stratum with base and ruv NES
    base_df <- fgsea_df[fgsea_df$gene_set == "base",
                         c("stratum", "signature", shape_col, "NES", "padj")]
    colnames(base_df)[4:5] <- c("NES_base", "padj_base")

    ruv_df <- fgsea_df[fgsea_df$gene_set == "ruv",
                        c("stratum", "signature", "NES", "padj")]
    colnames(ruv_df)[3:4] <- c("NES_ruv", "padj_ruv")

    merged <- merge(base_df, ruv_df, by = c("stratum", "signature"), all = FALSE)
    if (nrow(merged) == 0) return(NULL)

    # Significance classification
    merged$sig_class <- "Neither"
    merged$sig_class[merged$padj_base < fgsea_fdr & merged$padj_ruv >= fgsea_fdr] <- "Base only"
    merged$sig_class[merged$padj_base >= fgsea_fdr & merged$padj_ruv < fgsea_fdr] <- "RUV only"
    merged$sig_class[merged$padj_base < fgsea_fdr & merged$padj_ruv < fgsea_fdr] <- "Both"

    sig_colors <- c(
        "Both"      = "#2ca02c",
        "Base only" = "#ff7f0e",
        "RUV only"  = "#1f77b4",
        "Neither"   = "grey70"
    )

    # Axis limits
    lim <- max(abs(c(merged$NES_base, merged$NES_ruv)), na.rm = TRUE) * 1.1

    p <- ggplot(merged, aes(x = NES_base, y = NES_ruv)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_point(aes(color = sig_class, shape = .data[[shape_col]]),
                   size = 2.5, alpha = 0.8) +
        scale_color_manual(values = sig_colors, name = "Significant in") +
        coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
        labs(
            title    = title,
            x        = "Base model NES",
            y        = "RUV model NES",
            shape    = if (has_class) "Class" else "Category",
            caption  = paste0("Dashed line = identity | FDR < ", fgsea_fdr)
        ) +
        theme_bw() +
        theme(
            plot.title    = element_text(size = 10, face = "bold"),
            legend.text   = element_text(size = 7),
            legend.title  = element_text(size = 8)
        )

    # Label points that differ substantially
    merged$label <- ifelse(
        abs(merged$NES_base - merged$NES_ruv) > 1 |
        merged$sig_class %in% c("Base only", "RUV only"),
        paste0(merged$signature, "\n(", merged$stratum, ")"),
        NA
    )
    labeled <- merged[!is.na(merged$label), ]
    if (nrow(labeled) > 0 && nrow(labeled) <= 20) {
        p <- p + ggrepel::geom_text_repel(
            data = labeled,
            aes(label = label),
            size = 2, max.overlaps = 15,
            segment.color = "grey60"
        )
    }

    return(p)
}


# =============================================================================
# FORA-based specificity heatmap functions
# =============================================================================

# Color palette for benchmark categories
benchmark_category_colors <- list(
    cell_identity = list(
        fill = "#2E8B57", outer = "#4CAF50", inner = "#C8E6C9",
        label = "Cell identity"),
    artifact = list(
        fill = "#D94040", outer = "#E05050", inner = "#FFCDD2",
        label = "Artifact"),
    subject_confounder = list(
        fill = "#E8870E", outer = "#F5A623", inner = "#FFE0B2",
        label = "Subject confounder")
)


#' Convert snake_case class names to readable labels for facet strips
#' @param class_name Character vector of class names
#' @return Human-readable labels
prettify_class <- function(class_name) {
    lookup <- c(
        cell_type_marker     = "Cell type markers",
        disease_marker_T1D   = "T1D markers",
        disease_marker_T2D   = "T2D markers",
        disease_marker_PDAC  = "PDAC markers",
        functional_pathway   = "Functional pathways",
        ambient_RNA          = "Ambient RNA",
        cell_cycle           = "Cell cycle",
        dissociation_stress  = "Dissociation stress",
        doublet              = "Doublet",
        ischemia_procurement = "Ischemia / procurement",
        platform_QC          = "Platform QC",
        age                  = "Age",
        sex                  = "Sex",
        BMI                  = "BMI",
        height               = "Height",
        ancestry             = "Ancestry",
        smoking              = "Smoking",
        alcohol              = "Alcohol",
        caffeine             = "Caffeine",
        metabolic_state      = "Metabolic state",
        T1D                  = "T1D",
        T2D                  = "T2D",
        metformin            = "Metformin",
        statin               = "Statin",
        GLP1_agonist         = "GLP-1 agonist",
        immunosuppressant    = "Immunosuppressant",
        chemotherapy         = "Chemotherapy",
        other_antidiabetic   = "Other antidiabetic"
    )
    ifelse(class_name %in% names(lookup), lookup[class_name],
           gsub("_", " ", class_name))
}


#' Shorten signature names for display by stripping common prefixes
#' @param sig_name Character vector of signature descriptions
#' @return Shortened names (max 65 chars)
shorten_sig <- function(sig_name) {
    sig_name %>%
        gsub("^Identity markers - ", "", .) %>%
        gsub("^Age - ", "", .) %>%
        gsub("^Sex - ", "", .) %>%
        gsub("^BMI - ", "", .) %>%
        gsub("^Height - ", "", .) %>%
        gsub("^Weight_BMI - ", "", .) %>%
        gsub("^Ancestry - ", "", .) %>%
        gsub("^Ambient RNA - ", "", .) %>%
        gsub("^Dissociation - ", "", .) %>%
        gsub("^Doublet - ", "", .) %>%
        gsub("^Cell Cycle - ", "", .) %>%
        gsub("^Ischemia - ", "", .) %>%
        gsub("^Technical - ", "", .) %>%
        gsub("^Disease - ", "", .) %>%
        gsub("^Medication - ", "", .) %>%
        gsub("^Metabolic - ", "", .) %>%
        gsub("^Lifestyle - ", "", .) %>%
        gsub("^PDAC ", "", .) %>%
        gsub("T1D upregulated - ", "T1D up: ", .) %>%
        gsub("T1D downregulated - ", "T1D dn: ", .) %>%
        gsub("T2D upregulated - ", "T2D up: ", .) %>%
        gsub("T2D downregulated - ", "T2D dn: ", .) %>%
        gsub("T2D altered - ", "T2D: ", .) %>%
        gsub("T2D - ", "T2D: ", .) %>%
        gsub("T1D signature - ", "T1D: ", .) %>%
        gsub("T2D signature - ", "T2D: ", .) %>%
        gsub("PDAC signature - ", "PDAC: ", .) %>%
        gsub("^Alcohol: ", "", .) %>%
        gsub("^Smoking: ", "", .) %>%
        gsub("^Caffeine: ", "", .) %>%
        substr(1, 65)
}


#' Collapse per-signature FORA results into per-class mega-results
#'
#' For each class x gene_set x direction x stratum, pools all individual
#' signature genes into one mega-signature and recomputes overlap counts.
#'
#' @param bench_fora  data.frame of FORA results with columns: class, gene_set,
#'                    direction, stratum, overlap_genes, n_overlap, etc.
#' @return data.frame with collapsed rows (one per class x gene_set x direction x stratum)
collapse_benchmark_by_class <- function(bench_fora) {
    if (is.null(bench_fora) || nrow(bench_fora) == 0) return(NULL)

    collapsed <- bench_fora %>%
        group_by(contrast, stratum, class, category, control_type,
                 gene_set, direction) %>%
        summarise(
            signature         = first(class),
            subcategory       = "collapsed",
            confidence        = "High",
            n_sigs            = n(),
            # Pool all overlapping genes across signatures in this class
            overlap_genes_raw = paste(overlap_genes, collapse = "|"),
            n_genes_sig       = sum(n_genes_sig, na.rm = TRUE),
            n_sig_in_universe = sum(n_sig_in_universe, na.rm = TRUE),
            pvalue            = min(pvalue, na.rm = TRUE),
            padj              = min(padj, na.rm = TRUE),
            n_degs            = first(n_degs),
            .groups           = "drop"
        ) %>%
        mutate(
            # Deduplicate pooled genes
            overlap_genes  = sapply(overlap_genes_raw, function(g) {
                genes <- unique(unlist(strsplit(g, "\\|")))
                genes <- genes[genes != ""]
                paste(genes, collapse = "|")
            }),
            n_overlap      = sapply(overlap_genes, function(g) {
                genes <- unlist(strsplit(g, "\\|"))
                genes <- genes[genes != ""]
                length(unique(genes))
            }),
            test              = "fora",
            statistic         = ifelse(n_sig_in_universe > 0,
                                        n_overlap / n_sig_in_universe, 0),
            statistic_name    = "overlap_fraction",
            NES               = NA_real_,
            is_collapsed      = TRUE
        ) %>%
        select(-overlap_genes_raw)

    collapsed
}


#' Build FORA specificity heatmap from pipeline benchmark data
#'
#' Produces the FORA overlap fraction heatmap with directional triangle
#' overlays and black diamond significance markers.
#'
#' Two modes controlled by x_var:
#'   - "gene_set" (default): x-axis = Base/RUV/RUV-only. For per-stratum plots.
#'   - "stratum":  x-axis = cell types. For summary plots across strata.
#'     Use gene_set_filter to select which model (base/ruv) to show.
#'
#' @param bench_fora       data.frame of FORA results (test == "fora")
#' @param control_type_val "positive" or "negative"
#' @param categories_filter Category values to include (e.g., "cell_identity")
#' @param title            Plot title
#' @param fill_color       Gradient high color
#' @param fill_limit       Upper limit for gradient (values above are squished)
#' @param collapsed        If TRUE, show class-level collapsed results
#' @param deg_levels       Character vector of gene_set levels for x-axis ordering
#'                          (only used when x_var = "gene_set")
#' @param x_var            "gene_set" or "stratum" — which variable on the x-axis
#' @param gene_set_filter  When x_var = "stratum", which gene_set to show
#'                          (e.g., "base" or "ruv")
#' @return list(plot, width, height) or NULL. Width/height are recommended PDF
#'         dimensions in inches, computed from fixed per-tile sizes.
plot_fora_heatmap <- function(bench_fora,
                               control_type_val = NULL,
                               categories_filter = NULL,
                               title = "Benchmark Specificity",
                               fill_color = "#2E8B57",
                               fill_limit = 0.5,
                               collapsed = FALSE,
                               deg_levels = c("base", "ruv", "ruv_only"),
                               x_var = "gene_set",
                               gene_set_filter = NULL) {

    if (is.null(bench_fora) || nrow(bench_fora) == 0) return(NULL)

    # Filter by control_type or categories
    df <- bench_fora
    if (!is.null(control_type_val) && "control_type" %in% colnames(df)) {
        df <- df[df$control_type == control_type_val, ]
    }
    if (!is.null(categories_filter)) {
        df <- df[df$category %in% categories_filter, ]
    }
    if (nrow(df) == 0) return(NULL)

    # When x_var = "stratum", filter to one gene_set
    if (x_var == "stratum" && !is.null(gene_set_filter)) {
        df <- df[df$gene_set == gene_set_filter, ]
    }

    # Only keep gene_sets that exist in data
    avail_gs <- intersect(deg_levels, unique(df$gene_set))
    if (length(avail_gs) == 0) return(NULL)
    df <- df[df$gene_set %in% avail_gs, ]

    # Compute overlap_frac from statistic column (already overlap/size)
    df$overlap_frac <- ifelse(df$statistic_name == "overlap_fraction",
                               df$statistic, df$n_overlap / df$n_sig_in_universe)
    df$overlap_frac[is.na(df$overlap_frac)] <- 0

    # Build sig_short and sig_label
    if (collapsed && "n_sigs" %in% colnames(df)) {
        # For collapsed mode, signature == class name, so use prettify_class
        df$sig_short <- prettify_class(df$signature)
        df$sig_label <- sprintf("%s (n=%d)", df$sig_short, df$n_sigs)
    } else {
        df$sig_short <- shorten_sig(df$signature)
        df$sig_label <- sprintf("%s (n=%d)", df$sig_short,
                                 ifelse(is.na(df$n_sig_in_universe), 0, df$n_sig_in_universe))
    }

    # Prepare category and class labels for faceting
    cat_labels <- c(
        cell_identity      = "Cell identity",
        artifact           = "Artifact",
        subject_confounder = "Subject confounder"
    )
    df$category_label <- factor(
        ifelse(df$category %in% names(cat_labels), cat_labels[df$category], df$category),
        levels = cat_labels[names(cat_labels) %in% df$category]
    )
    df$class_label <- prettify_class(df$class)

    # Determine x-axis variable
    if (x_var == "stratum") {
        df$x_val <- df$stratum
        x_levels <- sort(unique(df$x_val))
    } else {
        df$x_val <- df$gene_set
        x_levels <- avail_gs
    }

    # Split into all (tile fill), up (gold triangle), down (cyan triangle)
    plot_data <- df %>%
        filter(direction == "all") %>%
        mutate(x_val = factor(x_val, levels = x_levels))

    # Remove signatures with zero overlap everywhere
    nonzero_sigs <- plot_data %>%
        group_by(sig_label, class) %>%
        summarise(any_overlap = any(n_overlap > 0), .groups = "drop") %>%
        filter(any_overlap) %>%
        pull(sig_label)
    plot_data <- plot_data %>% filter(sig_label %in% nonzero_sigs)

    if (nrow(plot_data) == 0) return(NULL)

    # Directional data (padj < 0.05 only)
    make_dir_data <- function(dir_name) {
        df %>%
            filter(direction == dir_name, padj < 0.05) %>%
            mutate(x_val = factor(x_val, levels = x_levels)) %>%
            filter(sig_label %in% nonzero_sigs)
    }
    dir_up   <- make_dir_data("up")
    dir_down <- make_dir_data("down")

    # Overall significance markers
    sig_overall <- plot_data %>% filter(padj < 0.05)

    # Build text labels
    plot_data <- plot_data %>%
        mutate(
            count_str = sprintf("(%d/%d)", n_overlap,
                                ifelse(is.na(n_sig_in_universe), n_overlap, n_sig_in_universe)),
            combo_label = ifelse(n_overlap == 0, "",
                                  sprintf("%.0f%% %s", overlap_frac * 100, count_str))
        )

    # Determine faceting structure
    has_multiple_cats <- length(unique(plot_data$category_label)) > 1

    # Build strip colors for nested facets
    active_cats <- unique(as.character(plot_data$category_label))
    strip_bg_list <- list()
    for (cat_label in active_cats) {
        cat_key <- names(cat_labels)[cat_labels == cat_label]
        if (length(cat_key) > 0 && cat_key %in% names(benchmark_category_colors)) {
            col <- benchmark_category_colors[[cat_key]]
            strip_bg_list <- c(strip_bg_list, list(
                element_rect(fill = col[["outer"]])))
        } else {
            strip_bg_list <- c(strip_bg_list, list(
                element_rect(fill = "grey80")))
        }
    }
    # Inner strip backgrounds (class_label level)
    for (cat_label in active_cats) {
        cat_key <- names(cat_labels)[cat_labels == cat_label]
        inner_classes <- plot_data %>%
            filter(category_label == cat_label) %>%
            pull(class_label) %>% unique() %>% sort()
        if (length(cat_key) > 0 && cat_key %in% names(benchmark_category_colors)) {
            col <- benchmark_category_colors[[cat_key]]
            for (cls in inner_classes) {
                strip_bg_list <- c(strip_bg_list, list(
                    element_rect(fill = col[["inner"]])))
            }
        } else {
            for (cls in inner_classes) {
                strip_bg_list <- c(strip_bg_list, list(
                    element_rect(fill = "grey92")))
            }
        }
    }

    themed_strip <- ggh4x::strip_nested(background_y = strip_bg_list)

    # Determine text size based on number of columns
    n_cols <- length(x_levels)
    combo_text_size <- if (n_cols <= 3) 2.8 else if (n_cols <= 6) 2.2 else 1.8

    if (has_multiple_cats) {
        p <- ggplot(plot_data, aes(x = x_val,
                                    y = reorder(sig_label, desc(sig_label)),
                                    fill = overlap_frac)) +
            geom_tile(color = "white", linewidth = 0.4) +
            geom_text(aes(label = combo_label), size = combo_text_size, hjust = 0.5) +
            ggh4x::facet_nested(
                category_label + class_label ~ .,
                scales = "free_y", space = "free_y",
                switch = "y",
                nest_line = element_line(linewidth = 0.3),
                strip = themed_strip)
    } else {
        p <- ggplot(plot_data, aes(x = x_val,
                                    y = reorder(sig_label, desc(sig_label)),
                                    fill = overlap_frac)) +
            geom_tile(color = "white", linewidth = 0.4) +
            geom_text(aes(label = combo_label), size = combo_text_size, hjust = 0.5) +
            ggh4x::facet_nested(
                class_label ~ .,
                scales = "free_y", space = "free_y",
                switch = "y",
                nest_line = element_line(linewidth = 0.3),
                strip = themed_strip)
    }

    # X-axis labels
    if (x_var == "gene_set") {
        gs_labels <- c(base = "Base", ruv = "RUV", ruv_only = "RUV-only")
        p <- p + scale_x_discrete(position = "bottom", labels = gs_labels)
    } else {
        p <- p + scale_x_discrete(position = "bottom")
    }

    # Compute recommended PDF dimensions from data size
    n_cols <- length(x_levels)
    n_rows <- length(nonzero_sigs)
    n_facet_groups <- length(unique(plot_data$class_label))

    # Per-column width: wider for few columns (per-stratum), narrower for many (multi-stratum)
    tile_w   <- if (n_cols <= 3) 1.8 else if (n_cols <= 6) 1.3 else 1.0
    tile_h   <- 0.35  # per signature row
    # Strip width: use longest class_label + category_label for nested facets
    max_class_chars <- max(nchar(as.character(plot_data$class_label)), na.rm = TRUE)
    max_cat_chars   <- if (has_multiple_cats) max(nchar(as.character(plot_data$category_label)), na.rm = TRUE) else 0
    strip_w  <- max(max_class_chars * 0.09, 1.2) + max(max_cat_chars * 0.09, 0) + 0.6
    strip_w  <- min(strip_w, 4.0)  # cap at 4 inches
    legend_w <- 1.2   # right legend
    title_h  <- 0.8   # title + top margin
    caption_h <- 0.6  # caption + bottom margin
    xlab_h   <- 1.0   # x-axis labels (angled)
    facet_gap <- 0.15 # gap per facet group boundary

    rec_w <- strip_w + n_cols * tile_w + legend_w
    rec_h <- title_h + n_rows * tile_h + n_facet_groups * facet_gap + xlab_h + caption_h

    # Apply minimums
    rec_w <- max(rec_w, 8)
    rec_h <- max(rec_h, 4)

    # Caption: two lines to avoid clipping
    caption_text <- paste0(
        "FORA (fgsea::fora)  |  Diamond = padj<0.05 (all DEGs)\n",
        "Gold/cyan triangles = up/down enrichment (padj<0.05)  |  (overlap / sig size)")

    p <- p +
        scale_fill_gradient(low = "white", high = fill_color,
                            limits = c(0, fill_limit),
                            oob = scales::squish, name = "Overlap\nfraction",
                            labels = scales::percent) +
        labs(title = title, x = "", y = "", caption = caption_text) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x  = element_text(angle = 35, hjust = 1, face = "bold", size = 9),
            axis.text.y  = element_text(size = 7.5),
            strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
            strip.placement = "outside",
            plot.title   = element_text(face = "bold", size = 13),
            plot.caption = element_text(size = 7, hjust = 0),
            legend.position = "right",
            panel.spacing = unit(0.15, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin  = margin(8, 12, 8, 8)
        )

    # Add significance markers — diamond at left edge of tile
    if (nrow(sig_overall) > 0) {
        p <- p + geom_point(
            data = sig_overall,
            aes(x = x_val, y = sig_label),
            shape = 18, size = 2.5, color = "black",
            inherit.aes = FALSE,
            position = position_nudge(x = -0.38)
        )
    }

    # Handle both-direction offset
    if (nrow(dir_up) > 0 && nrow(dir_down) > 0) {
        both_keys <- inner_join(
            dir_up %>% select(x_val, sig_label),
            dir_down %>% select(x_val, sig_label),
            by = c("x_val", "sig_label")
        )
        dir_up$has_both <- paste(dir_up$x_val, dir_up$sig_label) %in%
            paste(both_keys$x_val, both_keys$sig_label)
        dir_down$has_both <- paste(dir_down$x_val, dir_down$sig_label) %in%
            paste(both_keys$x_val, both_keys$sig_label)
    } else {
        if (nrow(dir_up) > 0) dir_up$has_both <- FALSE
        if (nrow(dir_down) > 0) dir_down$has_both <- FALSE
    }

    # Directional triangles at right edge of tile
    if (nrow(dir_up) > 0) {
        p <- p + geom_point(
            data = dir_up,
            aes(x = x_val, y = sig_label),
            shape = 24, size = 2.0, fill = "#FFD700", color = "black",
            stroke = 0.4, inherit.aes = FALSE,
            position = position_nudge(x = 0.35,
                                       y = ifelse(dir_up$has_both, 0.12, 0))
        )
    }
    if (nrow(dir_down) > 0) {
        p <- p + geom_point(
            data = dir_down,
            aes(x = x_val, y = sig_label),
            shape = 25, size = 2.0, fill = "#00CED1", color = "black",
            stroke = 0.4, inherit.aes = FALSE,
            position = position_nudge(x = 0.35,
                                       y = ifelse(dir_down$has_both, -0.12, 0))
        )
    }

    return(list(plot = p, width = rec_w, height = rec_h))
}