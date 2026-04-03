# Installation

## Option 1: Conda/micromamba (recommended)

```bash
cd /path/to/EasyDE
micromamba env create -f installation/EasyDE_install.yml
micromamba activate EasyDE
```

This installs R, DESeq2, RUVSeq, fgsea, Snakemake, and all dependencies.
Takes ~5-10 minutes on first run.

> Using conda/mamba? Replace `micromamba` with `conda` or `mamba`.

### Verify

```bash
Rscript -e "cat('R ok\n')"        # should print "R ok"
snakemake --version                 # should print 7.x
```

---

## Option 2: Install micromamba first (HPC / fresh system)

```bash
# Install micromamba (user-level, no root needed)
wget -O ~/micromamba.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest
tar -xvjf ~/micromamba.tar.bz2 -C ~/
rm ~/micromamba.tar.bz2

# Initialize shell integration
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc

# Configure channels
micromamba config append channels bioconda
micromamba config append channels conda-forge
micromamba config set channel_priority strict

# Then proceed with Option 1
micromamba env create -f installation/EasyDE_install.yml
```

---

## Post-Installation

EasyDE's conda environment includes all required packages. No manual
post-installation steps are needed under normal conditions.

### Known issues

**Conda R exit code 255**
Some conda-installed R builds return exit code 255 on completion regardless
of success. Test with `Rscript -e "cat('hello\n')"` — if exit code is not 0,
reinstall R:
```bash
micromamba install -n EasyDE -c conda-forge r-base --force-reinstall
```

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| R (r-base) | Runtime for all analysis scripts |
| DESeq2 | Differential expression (negative binomial GLM) |
| RUVSeq + EDASeq | Removal of unwanted variation |
| fgsea | Fast gene set enrichment analysis |
| apeglm | Log fold change shrinkage |
| Snakemake | Pipeline orchestration |
| r-yaml, r-optparse | Config parsing, CLI arguments |
| r-tidyverse | Data wrangling |
| r-ggplot2, r-ggrepel | Visualization |
