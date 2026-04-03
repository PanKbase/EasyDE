# =============================================================================
# filter_utils.R
# Gene and sample filtering functions for EasyDE
#
# All functions here are PURE: they take inputs, return outputs, no side effects
# (no reading files, no writing files, no printing — callers handle that)
# =============================================================================


# -----------------------------------------------------------------------------
# Gene-level filters
# These decide which genes are "expressed enough" to be worth testing
# -----------------------------------------------------------------------------

#' Proportion-based gene filter
#'
#' Keeps genes with >= min_reads counts in >= min_prop FRACTION of samples.
#' Applied row-wise via apply() on the count matrix.
#'
#' @param row       Numeric vector of counts for one gene across samples
#' @param min_reads Minimum count threshold (default 5)
#' @param min_prop  Minimum proportion of samples that must meet threshold (default 0.25)
#' @return Logical scalar — TRUE means keep this gene
proportion_filter <- function(row, min_reads = 5, min_prop = 0.25) {
    mean(row >= min_reads) >= min_prop
}


#' Count-based gene filter
#'
#' Keeps genes with >= min_reads counts in >= min_count NUMBER of samples.
#' Differs from proportion_filter: threshold is an absolute sample count, not a fraction.
#' Applied row-wise via apply() on the count matrix.
#'
#' @param row       Numeric vector of counts for one gene across samples
#' @param min_reads Minimum count threshold (default 5)
#' @param min_count Minimum NUMBER of samples that must meet threshold (default 3)
#' @return Logical scalar — TRUE means keep this gene
count_filter <- function(row, min_reads = 5, min_count = 3) {
    sum(row >= min_reads) >= min_count
}


#' Group-aware gene filter
#'
#' For categorical contrasts: keeps genes expressed in >= floor(n/2) samples
#' in EITHER the control OR experimental group (union logic).
#' For numeric contrasts: keeps genes expressed in >= floor(n/2) samples overall.
#'
#' Uses rowSums() instead of a for loop — operates on the whole matrix at once
#' in C, making it orders of magnitude faster on large gene x sample matrices.
#'
#' @param raw_mat          Count matrix (genes x samples), numeric
#' @param sample_metadata  Sample metadata data.frame
#' @param biosample_id_col Column name in sample_metadata holding biosample IDs
#'                         (the donor/subject unit — must match colnames of raw_mat)
#' @param contrast_var     Column name in sample_metadata for the contrast variable
#' @param control_grp      Value in contrast_var that defines the control group (or NULL for numeric)
#' @param experimental_grp Value in contrast_var that defines the experimental group (or NULL for numeric)
#' @param min_reads        Minimum count to consider a gene "expressed" in a sample (default 5)
#' @return Character vector of gene names that pass the filter
group_aware_gene_filter <- function(raw_mat,
                                    sample_metadata,
                                    biosample_id_col,
                                    contrast_var,
                                    control_grp      = NULL,
                                    experimental_grp = NULL,
                                    min_reads        = 5) {

    is_categorical <- is.factor(sample_metadata[[contrast_var]]) ||
                      is.character(sample_metadata[[contrast_var]])

    if (is_categorical) {

        ctrl_samps <- sample_metadata[sample_metadata[[contrast_var]] == control_grp, biosample_id_col]
        exp_samps  <- sample_metadata[sample_metadata[[contrast_var]] == experimental_grp, biosample_id_col]

        ctrl_cutoff <- floor(length(ctrl_samps) / 2)
        exp_cutoff  <- floor(length(exp_samps)  / 2)

        # rowSums on a logical matrix counts TRUEs per row (i.e. per gene)
        # This replaces the original for-loop — same logic, ~100x faster
        ctrl_pass <- rowSums(raw_mat[, ctrl_samps, drop = FALSE] >= min_reads) >= ctrl_cutoff
        exp_pass  <- rowSums(raw_mat[, exp_samps,  drop = FALSE] >= min_reads) >= exp_cutoff

        # Union: keep gene if it passes in EITHER group
        genes_to_keep <- rownames(raw_mat)[ctrl_pass | exp_pass]

    } else {

        all_cutoff    <- floor(nrow(sample_metadata) / 2)
        all_pass      <- rowSums(raw_mat >= min_reads) >= all_cutoff
        genes_to_keep <- rownames(raw_mat)[all_pass]

    }

    return(unique(genes_to_keep))
}


# -----------------------------------------------------------------------------
# Sample-level filters
# These decide which samples survive into the analysis
# -----------------------------------------------------------------------------

#' Apply an arbitrary list of subset filters to a metadata table
#'
#' This replaces the hardcoded subset_on1/subset_on2 pattern.
#' Instead, pass a named list of column -> allowed values pairs.
#' Any number of filters can be applied.
#'
#' Example:
#'   filters <- list(
#'       treatment      = "no_treatment",
#'       tissue         = c("islet", "exocrine"),
#'       disease_status = "ND"
#'   )
#'   sample_metadata <- apply_subset_filters(sample_metadata, filters)
#'
#' @param meta    data.frame of sample metadata
#' @param filters Named list: names are column names, values are allowed value(s)
#'                Pass an empty list (list()) to skip filtering entirely.
#'                Column names must exist in meta — unknown columns emit a warning
#'                and are skipped rather than crashing.
#' @return Filtered data.frame
apply_subset_filters <- function(meta, filters) {

    if (length(filters) == 0) return(meta)

    for (col_name in names(filters)) {

        allowed_values <- filters[[col_name]]

        if (!col_name %in% colnames(meta)) {
            warning(sprintf(
                "Filter column '%s' not found in metadata — skipping this filter.",
                col_name
            ))
            next
        }

        before <- nrow(meta)
        meta   <- meta[meta[[col_name]] %in% allowed_values, , drop = FALSE]
        after  <- nrow(meta)

        message(sprintf(
            "  Filter [%s = %s]: %d → %d samples",
            col_name, paste(allowed_values, collapse = " | "), before, after
        ))
    }

    return(meta)
}


#' Deduplicate donors: keep one biosample per donor
#'
#' When a donor has contributed multiple biosamples (e.g. different tissues,
#' runs, or time points), randomly retain one per donor before DE analysis.
#' Must be called with the donor column name from config (filtering.donor_col)
#' — there is no default, because assuming a column name would silently fail
#' on data where that column doesn't exist.
#'
#' Paired mode: Deduplication is SKIPPED. In paired treatment designs, each
#' donor has samples in both treatment arms — and may have multiple dishes
#' per treatment (genuine biological replicates from the same islet isolation).
#' The DESeq2 formula ~ donor + treatment blocks on donor while using
#' dish-to-dish variability as the error term. Dedup would discard power.
#'
#' @param meta       data.frame of sample metadata
#' @param donor_col  Column name holding donor/subject IDs (required, from config
#'                   filtering.donor_col). Distinct from biosample_id_col — this
#'                   is the biological individual, not the unit of analysis.
#' @param paired     Logical: if TRUE, skip dedup (paired design preserves all replicates)
#' @return Deduplicated data.frame (unpaired) or original data.frame (paired)
dedup_donors <- function(meta, donor_col, paired = FALSE) {

    if (!donor_col %in% colnames(meta)) {
        stop(sprintf(
            "donor_col '%s' not found in metadata. Check filtering.donor_col in config.yaml.",
            donor_col
        ))
    }

    # In paired mode, each donor has samples in BOTH treatment groups.
    # Multiple samples per donor per group are legitimate biological replicates
    # (e.g. separate dishes of islets from the same donor, treated independently).
    # The formula ~ donor_accession + treatment correctly blocks on donor while
    # using dish-to-dish variability as the error term. Do NOT dedup.
    if (paired) {
        message(sprintf("  Paired mode: skipping donor dedup (%d samples retained)", nrow(meta)))
        return(meta)
    }

    # Unpaired mode: keep 1 sample per donor regardless of group
    donor_counts <- table(meta[[donor_col]])

    # Donors with exactly 1 sample — keep as-is
    single <- meta[meta[[donor_col]] %in% names(donor_counts[donor_counts == 1]), ]

    # Donors with >1 sample — sample one randomly
    multi_donors <- names(donor_counts[donor_counts > 1])
    set.seed(123) # Ensure deterministic reproducible sampling
    multi_rows   <- lapply(multi_donors, function(d) {
        rows <- meta[meta[[donor_col]] == d, ]
        rows[sample(nrow(rows), 1), ]
    })

    result <- do.call(rbind, c(list(single), multi_rows))

    message(sprintf(
        "  Donor dedup: %d samples → %d donors retained",
        nrow(meta), nrow(result)
    ))

    return(result)
}
