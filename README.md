# Replication package for manuscript "Conspiracy Mentality and the Erosion of Democratic Accountability"

This folder contains the limited data and R scripts needed to reproduce the empirical results reported in the manuscript and appendix for the Australia survey experiment.

## Files

- `au_replication_limited_data.csv` and `au_replication_limited_data.rds`
  - Limited submitted-sample data with 2,013 rows.
  - The variable `main_analytic` marks the 1,227 respondents used in the manuscript's main models.
- `01_build_limited_data.R`
  - Rebuilds the supplied limited data from the original raw `.sav` file.
  - This script is retained for transparency and is not required to reproduce the submitted tables and figures.
- `02_replicate_main_results.R`
  - Reproduces the main vote-intention models, the main conspiracy-moderation figure, Appendix 1 descriptive outputs, and Appendix 2 balance checks.
  - Writes outputs to `output`.
- `03_replicate_mediation_results.R`
  - Reproduces the core mediation table and figure reported in the manuscript.
  - Writes outputs to `output`.
- `04_replicate_appendix_party_alignment.R`
  - Reproduces the appendix robustness checks that use detailed party-identification and feeling-thermometer alignment codings.
  - Writes outputs to `output`.
- `05_replicate_appendix_robustness.R`
  - Reproduces the sample-flow, excluded-profile, broader-sample, ordered-logit, covariate-adjusted, factor-loading, and additive-index appendix checks.
  - Writes outputs to `output`.

## Required R packages

The scripts use `tidyverse`, `haven`, `psych`, `janitor`, `broom`, `emmeans`, `flextable`, `mediation`, and `MASS`.

## How to run the replication

From this folder, run:

```r
source("02_replicate_main_results.R")
source("03_replicate_mediation_results.R")
source("04_replicate_appendix_party_alignment.R")
source("05_replicate_appendix_robustness.R")
```

Or from a terminal:

```bash
Rscript 02_replicate_main_results.R
Rscript 03_replicate_mediation_results.R
Rscript 04_replicate_appendix_party_alignment.R
Rscript 05_replicate_appendix_robustness.R
```

To rebuild the limited data from the original raw survey export on the authors' machine, run:

```bash
Rscript 01_build_limited_data.R
```

Readers do not need the raw `.sav` file if they use the supplied limited data.
