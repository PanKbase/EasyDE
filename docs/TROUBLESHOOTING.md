# Troubleshooting

## Common Errors

**`skipped_preflight` with "100% of genes have >= 1 zero"**
The count matrix is extremely sparse — every gene has at least one sample
with zero counts. When detected at preflight, the stratum's `zero_inflated`
flag is set to TRUE, and steps 03/05 automatically switch to
`DESeq(sfType = "poscounts")` instead of the default ratio estimator. The
poscounts method uses only positive counts for the geometric mean, allowing
normalization to proceed. The stratum is only rejected if even poscounts
cannot handle the sparsity level.

**`lfcShrink failed`**
apeglm requires at least one sample per group with non-zero counts. This is
a warning — `log2FoldChange_shrunk` will be NA but `log2FoldChange` is valid.

**`no count matrix found for stratum 'X'`**
The count file at `data/counts/X.counts.csv` (after safe-name conversion)
does not exist. Check that `strata` values in contrasts.csv match filenames
in `data/counts/`. Spaces and special characters become underscores.

**`W pruning removed all latent factors` / status `success_base_only`**
All W factors were excluded because every W correlated with the contrast
variable (raw p < 0.05). The stratum falls back to base DESeq2 results
automatically and is marked `success_base_only` in `contrast_summary.csv`.
This is expected when the biological signal dominates (e.g. disease-driven
W factors in highly affected cell types). Check the correlation heatmap in
`ruvseq_diagnostics.pdf` to confirm the exclusions are reasonable — excluded
W factors are labeled with "(excl.)" on the y-axis, and significance stars
(*, **, ***) mark significant associations. See
[Methods — W exclusion](METHODS.md#w-exclusion-strategy) for the rationale.

**fGSEA returns 0 significant pathways**
Not necessarily an error. Check: (1) the ranking statistic produces a
reasonable spread, (2) gene symbols in the GMT match your gene names,
(3) `fgsea.min_gene_set_size` is not too high.

**`X is not a factor` or `contrasts can be applied only to factors with 2 or more levels`**
These DESeq2 errors are now handled automatically by the pipeline. The first
is caused by logical columns (TRUE/FALSE) not being coerced to factors — step
02 now applies `as.factor()` to all logical contrast and covariate columns.
The second occurs when a covariate has only one observed level after sample
filtering — step 02 now performs a second-pass `droplevels()` check and drops
single-level covariates. If you still see these errors, verify that your
contrasts CSV and metadata are consistent.

**`error` with empty error message (typically QuiescentStellate or small strata)**
DESeq2 crashed with a multi-line error message where the first line was blank,
so the pipeline captured an empty string. The most common historical cause was
`n_samples == n_coefficients` after factor expansion — this is now caught at
preflight (step 02) and reclassified as `skipped_min_group` before DESeq2
runs. The preflight coefficient counter accounts for factor level expansion
(e.g. `Ethnicities` with 4 levels → 3 columns in `model.matrix()`). If you
still see this error, check the full log at `logs/{contrast}/{stratum}.log`.

**`skipped_min_group` for all strata in a contrast**
The experimental group has too few (or zero) positive samples across all
cell types. This is common for rare autoantibody markers (AAB_IA2, AAB_IAA,
AAB_ZNT8) and clinical measurements with high missingness (C_peptide,
HbA1c). This is a data limitation, not a pipeline bug — the contrast
definition in `contrasts.csv` is valid but the cohort doesn't have enough
positive cases.

**Snakemake shows "Nothing to be done"**
All output files already exist. Delete `results/` (or the specific contrast
directory) to force a re-run.

---

## Reading Log Files

Steps 02–07 write structured logs to `logs/{contrast_id}/{stratum}.log`:

```
[2026-03-06 14:32:01] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=start | ...
[2026-03-06 14:32:03] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=filter_samples | n_samples=89
[2026-03-06 14:32:15] [WARN]  contrast=Diabetic_vs_ND | biosample=Beta | step=lfc_shrinkage | msg=...
```

To see all warnings and errors after a run:

```bash
grep -rh "WARN\|ERROR\|SKIP" logs/
```

Step 08 also collects all flagged lines into
`results/{contrast}/summary/log_summary.tsv` for easy review.

---

## Resource Files

All resources are included in the repository — no additional downloads needed.

### Gene exclusion lists

Used by step 06 to remove ribosomal and mitochondrial genes before fGSEA ranking:

```
resources/gsea_files/rpl_genes.csv     ribosomal proteins, large subunit
resources/gsea_files/rps_genes.csv     ribosomal proteins, small subunit
resources/gsea_files/mtr_genes.csv     mitochondrially encoded genes
```

### Pathway GMT files

Three MSigDB gene set collections, pre-merged into a single default file:

```
resources/gsea_files/c2.cp.reactome.v*.symbols.gmt.txt    Reactome
resources/gsea_files/c2.cp.kegg_medicus.v*.symbols.gmt.txt KEGG Medicus
resources/gsea_files/c2.cp.kegg_legacy.v*.symbols.gmt.txt  KEGG Legacy
resources/gsea_files/reactome_kegg.gmt.txt                 merged (default)
```

To use a different collection, set `fgsea.gene_sets_file` in your config.

**Updating pathways:** Download new GMT files from
[MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp) and re-merge:

```bash
cat c2.cp.reactome.*.gmt c2.cp.kegg_medicus.*.gmt c2.cp.kegg_legacy.*.gmt \
    > resources/gsea_files/reactome_kegg.gmt.txt
```

### Benchmark signature files

Used by step 07 to benchmark DE results against curated gene signatures:

```
resources/benchmarking/master_gene_signatures.tsv    master file (277 signatures)
resources/benchmarking/Pancres_geneSig.txt            legacy positive controls
resources/benchmarking/artificats.txt                 legacy negative controls
```

The master file uses a hierarchical taxonomy (`control_type` → `category` →
`class` → `subcategory`). See [Methods — Signature Benchmarking](METHODS.md#signature-benchmarking-step-07)
for details. To use a different signature set, set `benchmarking.signature_file`
in your config.
