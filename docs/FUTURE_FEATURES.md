# Future Features

Ideas and designs for pipeline improvements — vetted but not yet implemented.

---

## 1. Progressive Covariate Dropping (Degree-of-Freedom Recovery)

**Status:** Designed, not implemented

### Problem

Small strata (e.g., QuiescentStellate with 8 samples) get skipped entirely
when `n_samples ≤ n_coefficients` after factor expansion. The covariates
`Age + Gender + BMI + Ethnicities + chemistry` expand to 7+ columns in
`model.matrix()`, leaving zero residual degrees of freedom. DESeq2 cannot
estimate dispersion and the stratum produces no results.

Currently these strata are classified as `skipped_min_group`. The biology
is lost even though a simpler model could still produce usable DE results.

### Proposed Solution

When `n_samples ≤ n_coefficients`, progressively drop covariates until the
model fits, rather than skipping the stratum entirely.

**Drop priority (highest cost → dropped first):**

1. **Most degrees of freedom consumed** — drop the covariate that frees the
   most columns. A factor with k levels consumes k−1 df (e.g., Ethnicities
   with 4 levels → 3 df). Numeric covariates consume 1 df each.
2. **Tie-breaking: near-singularity** — among covariates consuming equal df,
   drop the one closest to singular (highest % single level, or lowest
   variance for numeric covariates).
3. **Floor** — keep at minimum the contrast variable + intercept. All
   covariates are expendable if needed.

**Example — QuiescentStellate (8 samples, Prediabetic_vs_ND):**

```
Original formula: ~ Age + Gender + BMI + Ethnicities + chemistry + contrast
Expanded coefficients: 1 (intercept) + 1 + 1 + 1 + 3 (Ethnicities) + 1 + 1 = 9
n_samples = 8, n_coefficients = 9 → SKIP

After dropping Ethnicities (frees 3 df):
Reduced formula: ~ Age + Gender + BMI + chemistry + contrast
Expanded coefficients: 1 + 1 + 1 + 1 + 1 + 1 = 6
n_samples = 8, n_coefficients = 6 → 2 residual df → CAN RUN
```

### Transparency Requirements

- Log every dropped covariate with the reason (`"dropped Ethnicities:
  frees 3 df, n_samples (8) still <= n_coefficients (9)"`)
- Populate `dropped_covariates` column in `contrast_summary.csv`
- Add `reduced_formula` flag (TRUE/FALSE) to `contrast_summary.csv` so
  downstream analysis can filter or flag these strata
- Add `preflight_status = "warn"` (not "fail") — the stratum runs but
  with caveats
- The `model_info_*.csv` intermediate records the actual formula used,
  so the reduction is always traceable

### Design Considerations

- **Biological relevance ordering:** An alternative to the df-based drop
  order is a user-specified priority list in the config (e.g.,
  `covariate_priority: [Age, Gender, BMI, chemistry, Ethnicities]` where
  the last entry is dropped first). This gives domain experts control over
  which confounders are most important to retain. Could be combined with the
  df-based approach as a fallback when no priority is specified.

- **Minimum residual df:** Even after dropping, we should require at least
  `n_samples > n_coefficients + 1` (not just `>`) to ensure DESeq2 has
  enough residual df for dispersion estimation. A margin of 2–3 would be
  safer.

- **Interaction with RUVseq:** RUVseq adds W factors on top of the reduced
  formula. If covariates were already dropped to fit the base model, adding
  Ws may push the model back over the df limit. Step 04 already caps
  `max_k` based on available df, but this interaction should be tested.

- **Paired analysis:** In paired mode, the donor blocking factor consumes
  `n_donors − 1` df. Progressive dropping is even more valuable here since
  the donor term is non-negotiable.

### Where to Implement

- **`02_prepare_coldata.R` → `check_sample_counts()`**: After detecting
  `n_samples ≤ n_coefficients`, enter a loop that drops one covariate per
  iteration (highest df cost first), recalculates n_coefficients, and stops
  when `n_samples > n_coefficients + margin`. Return the reduced covariate
  list.
- **Downstream scripts (03–05)**: Already read the formula from
  `model_info_*.csv` — no changes needed if the reduced formula is written
  correctly by step 02.

---

*Add new feature proposals below this line.*
