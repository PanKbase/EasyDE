# EasyDE

A modular R pipeline for pseudobulk differential expression analysis with
RUVseq correction and gene set enrichment. Designed for scRNA-seq data
(PanKbase and similar resources) but works with any pseudobulk count matrices.

Orchestrated by Snakemake for local or SLURM cluster execution.

---

## Pipeline Overview

Nine steps. Steps 02–07 run in parallel across every
**contrast × stratum** (cell type) combination:

```
01  Validate config + contrasts          (once per run — soft gate)
02  Filter samples, build coldata        (per contrast × stratum)
03  Base DESeq2                          (per contrast × stratum)
04  RUVseq normalization + k selection   (per contrast × stratum)  ← skipped in paired mode
05  Final DESeq2 with W factors          (per contrast × stratum)  ← skipped in paired mode
06  fGSEA pathway enrichment             (per contrast × stratum)
07  Benchmark against gene signatures    (per contrast × stratum)  ← optional
08  Aggregate results + status plot      (once per contrast)
09  Pipeline summary heatmap             (once after all contrasts)
```

Step 02 includes **preflight validation** — checks for design matrix rank,
zero-inflation, sample sufficiency, and data quality issues. Rejected strata
are skipped cleanly through the rest of the pipeline, saving compute time.
The pipeline also handles several edge cases automatically: `poscounts`
normalization for sparse count matrices, factor coercion for logical columns,
second-pass covariate drops, and graceful RUV-to-base fallback. See
[Methods — Edge Cases](docs/METHODS.md#edge-cases-and-workarounds) for details.

---

## Quick Start

```bash
# 1. Install and activate
micromamba env create -f installation/EasyDE_install.yml
micromamba activate EasyDE

# 2. Fetch data from PanKbase (optional — skip if you have your own data)
Rscript workflow/scripts/fetch/pankbase_helpers.R \
    --config config/config_traits.yaml \
    --pankbase-config workflow/scripts/fetch/pankbase_config.yaml \
    --strata "Beta|Alpha|Delta"

# 3. Run the full pipeline
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml

# 4. Check results
cat results/Diabetic_vs_ND/summary/contrast_summary.csv | column -t -s,
```

---

## Project Layout

```
EasyDE/
├── Snakefile                             ← Snakemake workflow (do not edit)
├── config/                               ← EDIT THESE
│   ├── config_traits.yaml                ← unpaired analysis config
│   └── config_treatments.yaml            ← paired analysis config
├── contrasts/                            ← EDIT THESE
│   ├── contrasts.csv                     ← your contrast definitions
│   └── contrasts_template.csv            ← documented template to copy from
├── profiles/                             ← Snakemake execution profiles
│   ├── local/config.yaml
│   └── slurm/config.yaml
├── data/                                 ← YOUR DATA (or fetched from PanKbase)
│   ├── counts/                           ← one CSV per cell type
│   └── sample_metadata.csv
├── resources/
│   ├── gsea_files/                      ← pathway GMT files + gene exclusion lists
│   └── benchmarking/                    ← master gene signatures for step 07
├── workflow/scripts/                     ← pipeline code (do not edit)
│   ├── 01–09_*.R                         ← pipeline steps
│   ├── fetch/                            ← PanKbase data adapter
│   └── utils/                            ← shared R utilities
├── vignette.md                           ← step-by-step walkthrough
├── notebooks/                            ← post-run exploration
│   └── explore_results.ipynb
├── docs/                                 ← detailed documentation
└── installation/
    └── EasyDE_install.yml
```

| What to edit | Where |
|-------------|-------|
| Analysis parameters | `config/config_*.yaml` |
| Contrast definitions | `contrasts/contrasts.csv` |
| Execution settings | `profiles/*/config.yaml` |
| Data files | `data/counts/` + `data/sample_metadata.csv` |

---

## Documentation

| Guide | Contents |
|-------|----------|
| **[Vignette](vignette.md)** | Step-by-step walkthrough of a single contrast x stratum |
| **[Installation](docs/INSTALLATION.md)** | Environment setup, verification, known issues |
| **[Configuration](docs/CONFIGURATION.md)** | Config files, contrast definitions, data formats |
| **[Running](docs/RUNNING.md)** | Snakemake commands, SLURM setup, step-by-step manual run |
| **[Methods](docs/METHODS.md)** | DESeq2, RUVSeq (W exclusion strategy), fGSEA, signature benchmarking, paired analysis, preflight validation |
| **[Output](docs/OUTPUT.md)** | File tree, column descriptions, best-model file, status taxonomy |
| **[Troubleshooting](docs/TROUBLESHOOTING.md)** | Common errors, log reading, resource files |
| **[Future Features](docs/FUTURE_FEATURES.md)** | Designed improvements not yet implemented |

---

## Pipeline Status Values

The `contrast_summary.csv` produced by step 08 categorizes each stratum with
a fine-grained status taxonomy:

| Status | Meaning |
|--------|---------|
| `success` | RUVseq model ran successfully; full DE results produced |
| `success_base_only` | Base DESeq2 succeeded but RUV was skipped or failed (e.g. paired mode, RUV crash) |
| `skipped_no_samples` | Zero samples remained after filtering |
| `skipped_min_group` | Fewer than 3 control or experimental samples |
| `skipped_no_pairs` | No valid donor pairs found (paired contrasts only) |
| `skipped_preflight` | Design matrix rejected by preflight checks (rank-deficient, zero-inflated) |
| `skipped` | Other skip reason — check `error_message` column for details |
| `error` | DESeq2 or downstream step crashed — `error_message` populated |
| `not_run` | No log or output data found for this stratum |

Step 08 also generates a per-contrast **status overview heatmap** (`status_overview.pdf`)
showing status and DEG counts across all strata at a glance. After all contrasts
complete, step 09 generates a **pipeline-level status heatmap**
(`pipeline_status_overview.pdf`) spanning all contrasts × cell types.

---

## Exploratory Notebook

After running the pipeline, use the Jupyter notebook to visualize results:

```bash
jupyter notebook notebooks/explore_results.ipynb
```

The notebook includes:

| Section | Content |
|---------|---------|
| Load & boxplots | Load results including `_deseq_final.tsv`, DEG boxplots with model indicator (circle = RUVseq, triangle = Base) |
| fGSEA & status | Pathway heatmaps, divergent pathway analysis, fine-grained status overview |
| Export | Per-stratum files as `{stratum}--{contrast}--{model}.DESEQ2.tsv` |
| Ensembl mapping | biomaRt gene symbol → Ensembl ID conversion |
| PankBase manifest | Submission manifest with NFS paths, S3 URLs, file aliases, MD5 placeholders |

