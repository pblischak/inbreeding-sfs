# SFS-Based Demographic Inference with Inbreeding

This repository contains code, data, and results for all simulations and empirical
data analyses conducted in the manuscript *Inferring the Demographic History of
Inbred Species From Genome-Wide SNP Frequency Data*.

### `sims/`

The `sims/` folder contains code for performing the four main simulation experiments that
were conducted in the paper. Results for these simulations are provided as well.

### `data/`

The `data/` folder contains two example data sets (cabbage and puma), as well as
scripts for fitting demographic models, estimating parameter uncertainties, and plotting
comparisons between the observed and expected SFS.

### `bbc-shiny/`

The `bbc-shiny/` folder contains R code for a small Shiny application to demonstrate the
properties of the beta-binomial convolution that is used to derive the expected SFS with
inbreeding.

#### Setup

These analyses can be recreated using the latest version of ∂a∂i (v2.0.3). 

```bash
# 0. Add the Bioconda channel if you don't already have it
conda config --append channels bioconda

# 1. Create a new Python 3 environment for dadi
conda create -n dadi-env python=3

# 2. Activate the new environment
conda activate dadi-env

# 3. Install dadi (and dependencies)
conda install dadi
```

Figures were created either using Python and the built-in plotting functions in ∂a∂i or
using `ggplot2` within the `tidyverse` package in R v3.6.