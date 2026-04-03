# Running the Pipeline

## Running with Snakemake

Snakemake is the recommended way to run EasyDE. It handles the full DAG of
dependencies, parallelizes steps 02–07 across all contrast × stratum
combinations, and automatically retries failed jobs.

### Run all contrasts

```bash
micromamba activate EasyDE

snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

This runs **all** contrasts defined in your contrasts CSV across **all**
their strata. For 20 contrasts × 11 cell types, that is ~1,100 parallel jobs.

### Run one contrast

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml \
    results/Diabetic_vs_ND/summary/contrast_summary.csv
```

### Dry run (see what would execute)

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml -n
```

### Run treatment (paired) analysis

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_treatments.yaml
```

In paired mode the pipeline flow is **01 → 02 → 03 → 06 → 07 → 08 → 09** (steps 04
and 05 are automatically skipped — see [Methods](METHODS.md#paired-analysis)).
Benchmarking (step 07) is optional and only runs when `benchmarking.signature_file`
is configured.

### Clean and re-run

```bash
# Remove results for one contrast
rm -rf results/Diabetic_vs_ND/

# Remove all results
rm -rf results/ logs/

# Then re-run
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

---

## Running on SLURM

For cluster execution, use the SLURM profile. Edit
`profiles/slurm/config.yaml` to set your account, partition, and QOS:

```yaml
# profiles/slurm/config.yaml (edit these)
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
  slurm_qos: your_qos
  mem_mb: 8000
  time: 240      # minutes
  cpus: 1
```

Then run:

```bash
snakemake --profile $PWD/profiles/slurm \
    --config pipeline_config=config/config_traits.yaml
```

> **Note:** Use `$PWD/profiles/slurm` (absolute path) for the SLURM profile.
> SLURM jobs inherit the working directory, so all paths resolve correctly.

### SLURM with a cluster config

For production runs with absolute output paths (e.g. writing results to a
shared project directory), create a cluster-specific config:

```bash
snakemake --profile $PWD/profiles/slurm \
    --config pipeline_config=config/TSCC_config_Run_452026_allTraits.yaml
```

---

## Running Step by Step

For debugging or educational purposes, you can run each script manually.
All commands assume you are in the pipeline root directory.

```bash
micromamba activate EasyDE

CONFIG="config/config_traits.yaml"
CONTRAST="Diabetic_vs_ND"
STRATUM="Beta"

# Step 01 — Validate (once per run)
Rscript workflow/scripts/01_prepare_inputs.R --config "$CONFIG"

# Step 02 — Prepare coldata
Rscript workflow/scripts/02_prepare_coldata.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 03 — Base DESeq2
Rscript workflow/scripts/03_run_deseq_base.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 04 — RUVseq (skipped automatically in paired mode)
Rscript workflow/scripts/04_run_ruvseq.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 05 — Final DESeq2 with W factors (skipped automatically in paired mode)
Rscript workflow/scripts/05_run_deseq_final.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 06 — fGSEA
Rscript workflow/scripts/06_run_fgsea.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 07 — Benchmark signatures (optional — requires benchmarking.signature_file in config)
Rscript workflow/scripts/07_benchmark_signatures.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 08 — Aggregate (once per contrast, no --stratum needed)
Rscript workflow/scripts/08_aggregate_results.R \
    --config "$CONFIG" --contrast "$CONTRAST"

# Step 09 — Pipeline summary (once after ALL contrasts complete)
Rscript workflow/scripts/09_pipeline_summary.R --config "$CONFIG"
```
