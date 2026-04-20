# Methods

## DESeq2 (steps 03 and 05)

Standard DESeq2 pseudobulk analysis. Step 03 fits the base model with
user-specified covariates. If RUVseq is enabled (unpaired mode only), step 05
re-fits the model with additional W latent factors to control for unwanted
technical variation. Both steps produce log2 fold changes (MLE and
apeglm-shrunk) and BH-adjusted p-values.

For **continuous** contrasts (e.g. Age, BMI, HbA1c), the pipeline detects
that `control_grp` and `experimental_grp` are blank and treats the contrast
variable as numeric. Log2 fold changes represent the effect per unit increase.

### Best-model consolidation (step 08)

Step 08 produces a consolidated file (`{contrast}_deseq_final.tsv`) that
selects the best available model per stratum:

- **RUVseq results** for strata with status `success`
- **Base DESeq2 results** as fallback for `success_base_only` strata

Each row includes a `model` column ("RUVseq" or "Base") and a `formula`
column with the exact design formula used, so downstream analyses can track
which correction was applied per cell type.

---

## Paired Analysis

When `paired_analysis: true` is set in the config, the pipeline:

1. **Matches samples by donor** across treatment and control conditions
   using the `donor_col` specified in the config.
2. **Includes donor as a blocking factor** in the DESeq2 formula
   (e.g. `~ donor_accession + treatment`), which accounts for individual
   variability.
3. **Skips RUVseq (steps 04 and 05)** because the donor blocking factor
   already captures inter-individual variation that RUVseq's latent
   factors would model. Running RUVseq on top of donor blocking is
   redundant and can absorb genuine biological signal, reducing the
   number of detected DEGs.

The paired-mode pipeline flow is: **01 → 02 → 03 → 06 → 07 → 08 → 09**.

---

## RUVseq (step 04)

> **Note:** RUVseq is automatically skipped in paired analysis mode.

RUVseq removes unwanted technical variation using empirical negative control
genes — genes that are not differentially expressed (base DESeq2 p-value
above `ruvseq.empirical_pval_threshold`, default 0.5). These genes capture
technical batch effects without biological signal.

**Automatic k selection:** The pipeline runs `RUVg` for k = 1 through
`ruvseq.max_latent_factors` (default 10). For each k, it computes the
Relative Log Expression (RLE) matrix and measures between-sample variance via
a one-way ANOVA F-statistic. The optimal k is chosen at the "elbow" of the
F-statistic curve — the point where the second derivative is maximized,
indicating the transition from meaningful variance reduction to diminishing
returns. This is analogous to the scree-plot elbow method used in PCA.

### W exclusion strategy

After selecting the best k, the pipeline identifies which W factors are
associated with the **contrast variable** and should be excluded from the final
model.

**Association testing — omnibus F-test:**
Each W is tested against each covariate using an omnibus F-test (ANOVA: full
model `lm(W ~ variable)` vs null model `lm(W ~ 1)`). This tests whether the
covariate explains *any* variance in W. For numeric covariates and 2-level
factors, this is equivalent to the standard t-test (F = t^2). For multi-level
factors (e.g., Ethnicities with 3+ levels), the F-test evaluates all levels
jointly — avoiding the pitfall of testing only one arbitrary contrast
determined by R's alphabetical factor ordering. The heatmap displays signed
t-statistics from the strongest coefficient for visual interpretation, while
exclusion decisions use the F-test p-value.

**Contrast variable exclusion (raw p < 0.05):**
Any W with F-test raw p < 0.05 against the contrast variable is excluded.
A raw (unadjusted) threshold is used because the cost of errors is asymmetric:
a false-positive exclusion (removing a W that was actually technical) only costs
a modest loss of power, while a false-negative (keeping a disease-correlated W)
lets the latent factor absorb real biological signal — a far worse outcome.

**Why covariates are NOT used for exclusion:**
Covariates (Age, Sex, BMI, Ethnicities, etc.) are already included in the
DESeq2 design formula. A W factor that correlates with a covariate is not
confounding — the covariate is already controlled for. Excluding such a W
would unnecessarily reduce the pipeline's ability to remove technical
variation. Covariate correlations are shown on the diagnostic heatmap for
informational purposes only.

If all W factors are excluded, the RUV model collapses to the base model
(status becomes `success_base_only`). The `ruvseq_summary.csv` intermediate
records `best_k`, `contrast_associated_Ws`, and `safe_Ws` for full
transparency.

### Correlation heatmap

The diagnostic correlation heatmap (`ruvseq_diagnostics.pdf`) visualizes the
relationship between W factors and all known variables:

- Excluded W factors labeled with **(excl.)** suffix on the y-axis (e.g. "W_1 (excl.)")
- Significance stars (*, **, ***) based on BH-adjusted p-values (per-W, informational)
- Caption clarifies: exclusion is based on contrast variable at raw p < 0.05; covariate correlations are informational only

---

## Sample and Gene Filtering (step 02)

Step 02 (`02_prepare_coldata.R`) applies two distinct stages of filtering before dispersion estimation: **QC filters** remove low-quality samples and lowly-expressed genes, and **analytical restrictions** define which donors enter a given comparison.

### QC filters

**Samples and donors:**

- **Minimum cells per sample** (`filtering.min_cells_per_sample`; default 20). Biosamples aggregating fewer than 20 single cells, or lacking a cell count, are discarded so that pseudobulk profiles are statistically stable.
- **Donor deduplication** (unpaired mode). If a donor contributed multiple biosamples, one is retained at random (seeded for reproducibility) to avoid pseudoreplication. Skipped in paired mode, where repeat biosamples are genuine replicates.
- **Preflight sample-size gates.** Strata with fewer than 3 surviving samples are marked `skipped_no_samples`. For categorical contrasts, fewer than 3 samples per group — or a sample count not exceeding the number of model coefficients — triggers `skipped_min_group`. In paired mode, strata with no donors present in both arms are marked `skipped_no_pairs`.

**Genes:**

A group-aware low-count filter (`filtering.min_gene_counts`; default 5) is applied before size-factor and dispersion estimation. For categorical contrasts, a gene is retained if it has at least `min_gene_counts` reads in at least half the samples of *either* group (union rule), preserving condition-specific signals such as induced genes, lineage markers, or cytokine responses that a symmetric threshold would discard. For continuous contrasts, a gene must meet the count threshold in at least half of all samples.

### Analytical restrictions

After QC, the cohort is narrowed to donors relevant to the specific comparison:

- **Contrast group restriction.** For categorical contrasts, only samples in `control_grp` or `experimental_grp` are retained.
- **User-defined metadata filters** (`filter_col_N` / `filter_val_N` in the contrasts CSV, up to three pairs). Applied sequentially to restrict the analysis to a metadata stratum — for example, untreated donors, a single sex, or a specific tissue source.
- **Paired-design matching** (paired mode). Only donors present in both arms are kept.
- **Count-matrix alignment.** The count matrix is reduced to the intersection of its columns with surviving biosample IDs; metadata-only samples are logged and dropped.

Cell-type labels are matched case- and punctuation-insensitively, so `"Gamma+Epsilon"` and `"Gamma_Epsilon"` are equivalent.

---

## Edge Cases and Workarounds

The pipeline detects and handles several edge cases that commonly arise in
pseudobulk scRNA-seq data — especially with rare cell types, small cohorts,
and sparse autoantibody markers. The goal is to **fail early and cheaply**
when a stratum lacks statistical power, and **work around recoverable issues**
rather than crashing.

### Early termination (resource-saving skips)

These checks run in step 02 (preflight) before any DESeq2 computation,
saving compute time on strata that cannot produce meaningful results:

| Check | Condition | Status assigned |
|-------|-----------|-----------------|
| No samples | Zero samples after metadata filtering | `skipped_no_samples` |
| Minimum group size | Fewer than 3 control or experimental samples | `skipped_min_group` |
| Sample–coefficient ratio | `n_samples ≤ n_coefficients` (after factor expansion) | `skipped_min_group` |
| Rank-deficient design | `model.matrix()` rank < number of columns | `skipped_preflight` |
| No donor pairs | Paired mode: no donors in both conditions | `skipped_no_pairs` |

These skips are recorded with specific error messages in `contrast_summary.csv`
so the reason is always traceable.

**Coefficient counting and factor expansion:** The sample–coefficient check
counts the actual number of columns DESeq2 will create in `model.matrix()`,
not the number of formula terms. A factor covariate with k levels expands to
k−1 indicator columns (e.g., `Ethnicities` with 4 levels → 3 columns). This
matters for small strata: a model with 5 covariates might need only 7 columns
if all are binary, but 10+ if one is a multi-level factor. The pipeline counts
expanded columns at preflight time so strata where `n_samples ≤ n_coefficients`
are caught before DESeq2 runs, rather than crashing with "estimation of
dispersion is not possible."

### Zero-inflation / poscounts workaround

**Problem:** In pseudobulk data, small or rare cell types can produce count
matrices where every gene has at least one zero across samples. DESeq2's
default size-factor estimation (`sfType = "ratio"`) computes the geometric
mean per gene — which is zero when any sample is zero, making normalization
impossible.

**Detection:** Step 02 runs `check_zero_inflation()`, which flags the stratum
as `zero_inflated = TRUE` in the preflight report when 100% of genes have
≥1 zero.

**Workaround:** Steps 03 and 05 read the preflight flag. If `zero_inflated`
is TRUE, they call `DESeq(dds, sfType = "poscounts")` instead of the default.
The `poscounts` estimator uses only positive counts for the geometric mean,
sidestepping the all-zeros problem. This is DESeq2's recommended approach
for sparse data.

### Factor coercion for logical columns

**Problem:** Autoantibody columns (AAB_GADA, AAB_ZNT8, etc.) are stored as
TRUE/FALSE in the metadata. DESeq2 requires factors for categorical variables
and crashes with "X is not a factor" on logicals.

**Fix:** Step 02 coerces all logical contrast and covariate columns to factors
with `as.factor()` before building the design matrix.

### Second-pass covariate validation

**Problem:** A covariate like Gender may have 2+ levels in the full dataset
but collapse to a single level after sample filtering (e.g., only male donors
in a small stratum). DESeq2 crashes with "contrasts can be applied only to
factors with 2 or more levels."

**Fix:** After the first-pass covariate drop, step 02 calls `droplevels()`
on each factor covariate and removes any with fewer than 2 observed levels.
The dropped covariates are recorded in `preflight.csv` and
`contrast_summary.csv`.

### Graceful RUVseq failure → base fallback

When RUVseq (step 04) or the final DESeq2-with-Ws model (step 05) crashes —
due to singular fits, rank deficiency after adding W factors, or other
numerical issues — the pipeline falls back to the base DESeq2 results from
step 03. The stratum is marked `success_base_only` and included in the
consolidated `_deseq_final.tsv` file with `model = "Base"`. No results are
lost.

### Skip sentinel system

Snakemake requires all declared output files to exist. When a stratum is
skipped or errors out, the pipeline writes placeholder "sentinel" files
(containing `skipped=TRUE`) for each expected output. This satisfies
Snakemake's file contract without generating misleading result files.
Downstream steps detect sentinels and propagate the skip cleanly.

### LFC shrinkage fallback

`lfcShrink()` (apeglm) can fail when sample counts are very small or when
the coefficient name doesn't match DESeq2's internal naming. The pipeline:

1. Reconstructs the coefficient name using R's `make.names()` logic
2. Falls back to grepping `resultsNames(dds)` if the constructed name is absent
3. Wraps the call in `tryCatch()` — on failure, `log2FoldChange_shrunk` is
   set to NA while the MLE `log2FoldChange` remains valid

---

## fGSEA (step 06)

Fast Gene Set Enrichment Analysis using the Wald test statistic as ranking
metric. Ribosomal (RPL/RPS) and mitochondrial (MT-) genes are excluded
before ranking to prevent these highly variable housekeeping genes from
dominating pathway results.

Runs on base results, RUV-corrected results, or both (set `fgsea.run_on`
in config). In paired mode, fGSEA runs on base results only (since RUVseq
is skipped). Default gene set: merged Reactome + KEGG pathways from MSigDB.

---

## Signature Benchmarking (step 07)

Step 07 tests whether DE results recover known cell-type gene signatures.
It runs two complementary statistical tests and compares three gene sets,
answering the question: **does RUV correction improve, hurt, or preserve
recovery of known biology?**

### Tests

**FORA (overrepresentation analysis)** — uses `fgsea::fora()` (hypergeometric
test) to test whether DEGs overlap with signature genes more than expected by
chance. DEGs are split by direction before testing: each gene set is divided
into upregulated (LFC > 0) and downregulated (LFC < 0) subsets, and FORA is
run on each direction separately. This directly reveals whether a signature
is enriched among up- or down-regulated genes. Universe = all genes with a
non-NA adjusted p-value in the respective model's results.

**fGSEA (ranking-based enrichment)** — runs fast gene set enrichment analysis
using the Wald statistic as ranking metric. Threshold-free: uses the full
ranked gene list, not just DEGs. Returns a normalized enrichment score (NES)
indicating both direction and magnitude of enrichment. Uses
`fgseaMultilevel()` with `minSize = 2` to accommodate small curated
signatures (many are 3–12 genes).

### Gene sets compared

| Gene set | FORA subsets | fGSEA | Purpose |
|----------|-------------|-------|---------|
| `base` | `base_up`, `base_down` | Full base ranking | Baseline without correction |
| `ruv` | `ruv_up`, `ruv_down` | Full RUV ranking | After unwanted variation removal |
| `ruv_only` | `ruv_only_up`, `ruv_only_down` | — | Genes gained by correction |

fGSEA runs on base and RUV rankings only (not applicable to `ruv_only`
since fGSEA uses the full ranked gene list, not a subset).

### Master signature file

Signatures are provided as a single master TSV file with a hierarchical
taxonomy:

```
control_type → category → class → subcategory
```

| Level | Values |
|-------|--------|
| `control_type` | `positive`, `negative` |
| `category` | `cell_identity`, `artifact`, `subject_confounder` |
| `class` | ~34 mid-level buckets (e.g. `cell_type_marker`, `dissociation_stress`, `T1D`, `age`, `BMI`, `smoking`) |
| `subcategory` | Fine-grained labels |

Master file columns: `signature_id`, `tissue`, `cell_type`, `control_type`,
`category`, `class`, `subcategory`, `description`, `genes` (pipe-separated),
`num_genes`, `confidence`, `source_file`.

Signatures are automatically filtered to match the stratum being tested
(cell_type matching or `all_cells` entries) and the tissue specified in config.

**Positive controls** (cell_identity): cell-type identity markers, functional
pathways, known disease signatures (T1D, T2D, PDAC).

**Negative controls** (artifact + subject_confounder): dissociation stress,
ambient RNA, doublets, cell cycle, ischemia/procurement artifacts, plus
demographic (age, sex, BMI, ancestry) and clinical (medications, lifestyle,
comorbidities) confounders.

### Output

Per-stratum output is written to
`results/{contrast}/{stratum}/benchmarking/benchmark_signatures.tsv`.

| Column | Description |
|--------|-------------|
| `contrast` | Contrast ID |
| `stratum` | Cell type |
| `signature` | Signature name (from description column) |
| `control_type` | `positive` or `negative` |
| `category` | `cell_identity`, `artifact`, or `subject_confounder` |
| `class` | Mid-level class (e.g. `cell_type_marker`, `dissociation_stress`) |
| `subcategory` | Fine-grained subcategory |
| `confidence` | Curated confidence level (High, Medium, Low) |
| `n_genes_sig` | Total genes in the signature |
| `gene_set` | Which DEG set: `base`, `ruv`, or `ruv_only` |
| `direction` | `up` or `down` (FORA only; NA for fGSEA) |
| `test` | `fora` or `fgsea` |
| `statistic` | Overlap fraction (FORA) or NES (fGSEA) |
| `statistic_name` | `overlap_fraction` or `NES` |
| `NES` | Normalized enrichment score (fGSEA only; NA for FORA) |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |
| `n_overlap` | Genes overlapping between DEGs and signature (FORA) or signature genes found in ranking (fGSEA) |
| `n_sig_in_universe` | Signature genes present in the tested universe (FORA only) |
| `n_degs` | Number of DEGs in the directional set (FORA only) |
| `overlap_genes` | Pipe-separated list of overlapping or leading-edge genes |

### Benchmark visualization

**Per-stratum (step 07)** — in `{stratum}/benchmarking/`:
- **`benchmark_signatures.tsv`** — full FORA + fGSEA results vs all tested signatures

**Per-contrast (step 08)** — interactive HTML in `summary/`:
- **`positive_benchmarking.html`** — fGSEA pathway results across all cell types, organized by curated mega-sets; 3-level expandable drill-down (mega-set → pathway → leading-edge genes)
- **`negative_benchmarking.html`** — fGSEA enrichment of LLM-curated artifact and confounder signatures across all cell types; 3-level drill-down (class → individual signature → leading-edge genes)

Both HTML files are self-contained (no external dependencies) and sized ~6–20 MB depending on the number of cell types and pathways. See [Output — Interactive HTML Outputs](OUTPUT.md#interactive-html-outputs-step-08) for details on the structure and toggle behavior.

### Usage

Step 07 reads the signature file path from `benchmarking.signature_file` in
the config — no additional CLI flags needed:

```bash
Rscript workflow/scripts/07_benchmark_signatures.R \
    --config config/config_traits.yaml \
    --contrast ND_vs_T1D \
    --stratum Beta
```

The step is automatically skipped (with a clean sentinel file) when no
signature file is configured or the file does not exist.

---

## Preflight Validation (steps 01–02)

### Step 01 — Structural validation

Checks that all files exist, required columns are present, contrasts are
well-formed, and groups have matching samples. Per-contrast issues are
**flagged but not fatal** — only structural problems (missing CSV columns,
duplicate contrast IDs) halt the pipeline.

### Step 02 — Data-level validation

For each stratum, step 02 performs:

1. **Count matrix checks** — no negative values, warns on all-zero gene rows
2. **Zero-inflation detection** — if 100% of genes have ≥1 zero across
   samples, the stratum is rejected (no DESeq2 normalization method can
   handle this extreme sparsity)
3. **Design matrix rank check** — detects rank-deficient models
4. **Near-singularity detection** — warns on covariates with >95% one level
5. **Factor coercion** — logical columns (e.g. autoantibody TRUE/FALSE) are
   explicitly cast to factors before building the design matrix, preventing
   DESeq2's "X is not a factor" error
6. **Second-pass covariate validation** — after initial filtering, each
   factor covariate is re-checked with `droplevels()`. Covariates that
   collapse to a single level after sample filtering (e.g. Gender in a
   stratum with only male donors) are automatically dropped, preventing
   DESeq2's "contrasts can be applied only to factors with 2 or more levels"
   error

Results are written to `intermediates/preflight.csv` with columns including
`zero_inflated`, `pct_genes_with_zero`, `rank_ok`, and `near_singular_vars`.
Steps 03–07 read the preflight status and skip rejected strata gracefully.

See the [Edge Cases and Workarounds](#edge-cases-and-workarounds) section
above for a complete list of failure modes and how the pipeline handles each.
