# Configuration

EasyDE separates **what** to analyze from **how** to execute:

| File | Purpose |
|------|---------|
| `config/config_traits.yaml` | Unpaired analysis: thresholds, paths, step toggles, DESeq2/RUVseq/fGSEA parameters |
| `config/config_treatments.yaml` | Paired analysis: same structure, with `paired_analysis: true` and treatment-specific contrasts |
| `profiles/local/config.yaml` | Local execution: core count, keep-going, latency wait |
| `profiles/slurm/config.yaml` | SLURM execution: job count, sbatch template, account/partition |

The pipeline config (traits or treatments) controls the R scripts. The
profile controls Snakemake's execution engine. They are independent.

---

## Contrast Definitions

Copy the template and edit:

```bash
cp contrasts/contrasts_template.csv contrasts/contrasts.csv
```

Each row defines one contrast. Key columns:

| Column | Description | Example |
|--------|-------------|---------|
| `contrast_id` | Unique name | `Diabetic_vs_ND` |
| `contrast_var` | Metadata column to test | `Derived diabetes status` |
| `biosample_id_col` | Sample ID column | `sample_accession` |
| `control_grp` | Reference level (leave blank for continuous) | `Normal` |
| `experimental_grp` | Test level (leave blank for continuous) | `Diabetes` |
| `strata` | Cell types, pipe-separated | `Beta\|Alpha\|Delta` |
| `covariates` | Confounders, pipe-separated | `Age (years)\|Gender\|BMI` |
| `latent_corr_vars` | Variables to check for W-factor correlation | `Age (years)\|Gender\|BMI` |
| `filter_col_N` / `filter_val_N` | Subset filters (up to 3 pairs) | `treatment` / `no_treatment` |

See `contrasts/contrasts_template.csv` for full column documentation with
examples for binary, continuous, and filtered contrasts.

---

## Count Matrix Format

One CSV per stratum at `data/counts/{stratum_safe_name}.counts.csv`. The safe
name replaces spaces and special characters with underscores
(`"CD4 T cell"` becomes `CD4_T_cell.counts.csv`).

```
gene,SAMN001,SAMN002,SAMN003,...
INS,142,0,88,...
GCG,0,201,0,...
```

Genes as rows, samples as columns. First column is the gene symbol.

---

## Sample Metadata Format

`data/sample_metadata.csv` — one row per biosample, must contain:
- The `biosample_id_col` from your contrasts file (e.g. `sample_accession`)
- The contrast variable column (e.g. `Derived diabetes status`)
- All covariates listed in `covariates`
- All variables listed in `latent_corr_vars` (for RUVseq correlation heatmap)
- The donor ID column set in `filtering.donor_col` (e.g. `donor_accession`)

---

## Benchmarking Configuration

Signature benchmarking (step 07) is enabled when `benchmarking.signature_file`
points to a valid TSV file. Set it to `""` or remove the key to disable.

```yaml
benchmarking:
  # Master signature file (TSV). Columns: signature_id, tissue, cell_type,
  # control_type, category, class, subcategory, description, genes (pipe-separated),
  # num_genes, confidence, source_file
  signature_file: "resources/benchmarking/master_gene_signatures.tsv"

  # Tissue filter: only signatures matching this tissue (case-insensitive)
  # or with Tissue = "ALL" / "Pan-tissue" are used.
  tissue: "Pancreas"
```

The signature file uses a hierarchical taxonomy:
`control_type` (positive/negative) → `category` (cell_identity / artifact /
subject_confounder) → `class` (~34 mid-level buckets) → `subcategory`.

See [Methods — Signature Benchmarking](METHODS.md#signature-benchmarking-step-07)
for details on the statistical tests and output format.
