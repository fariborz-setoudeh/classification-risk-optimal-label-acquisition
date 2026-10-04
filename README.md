# Classification-Risk-Optimal Label Acquisition

R code and reproducibility materials for **Classification-Risk-Optimal Label Acquisition**, by Fariborz Setoudehtazangi and Geoffrey J. McLachlan.

The manuscript develops a label-acquisition criterion that combines conditional label information with classification-risk curvature to minimize the leading asymptotic coefficient of expected multiclass zero-one excess risk.

## Important: final analysis workflow

The final manuscript results use the validated workflow in scripts **07–10**.

The four original analysis scripts are retained because the validated workflow reads functions and model specifications from them. Running those original scripts directly does **not** reproduce the revised manuscript results.

Do not edit the original dependency scripts. The validated workflow checks their file contents against the reviewed versions.

## Repository files

The following scripts are located in the repository root:

| File | Role |
|---|---|
| `01_simulation_scenario1.R` | Original simulation dependency |
| `02_simulation_scenario2.R` | Original simulation dependency |
| `03_landsat_analysis.R` | Original Gaussian-analysis dependency |
| `03_drybean_analysis.R` | Original Dry Bean dependency |
| `validation_07_realdata_convergence.R` | Final real-data convergence and analysis workflow |
| `validation_08_realdata_figures.R` | Real-data figures from the audited saved results |
| `validation_09_simulation.R` | Final simulation convergence and analysis workflow |
| `validation_10_simulation_figures.R` | Simulation figures from the audited saved results |

The superseded scripts `04_landsat_figure.R` and `06_drybean_figures.R` are not part of the final workflow.

The Dry Bean dependency may instead be named `03_drybean_analysis.R.R`. The validated workflow accepts either name, but exactly one copy must be present.

## Numerical studies

All studies compare:

- Random acquisition;
- Entropy sampling;
- Margin sampling;
- Fisher acquisition using the complete-classification-information target;
- Adaptive risk-optimal acquisition.

The simulations additionally include an oracle risk-optimal design using the true model. This is a benchmark unavailable in real-data applications.

### Simulations

Two three-class quadratic discriminant analysis scenarios are considered:

1. Heterogeneous class covariance matrices.
2. Weak covariance heterogeneity, providing a near-linear control configuration.

Each scenario uses:

- 1,000 observations in the acquisition pool;
- an initial pilot of 60 labels;
- total labeling budgets of 10%, 20%, and 30%, including the pilot;
- 100 Monte Carlo replications;
- 100,000 independently generated test observations per replication.

Validation 09 implements the corrected fitting, curvature, optimization, and numerical checks.

Validation 10 reproduces the simulation performance figures from the frozen 100-replication results. Its acquisition-geometry illustrations use separate feature-only pools generated from the true models. These illustrations are not selected performance replications.

### Statlog Landsat Satellite

The supplied UCI benchmark split is retained:

- 4,435 training observations;
- 2,000 test observations;
- six represented classes.

Attributes 17–20 provide the four central-pixel spectral measurements.

The experiment uses a 5% random pilot, total labeling budgets of 10%, 20%, and 30%, and 50 pilot/acquisition replications. The supplied training and test samples remain fixed.

Principal components are used only to display acquisition geometry. Model fitting and acquisition use all four selected attributes.

### Dry Bean

Exact duplicate observations are removed before analysis.

Each replication uses a stratified 70/30 training/test split. Standardization and PCA are estimated from the training observations, and the first five principal components are retained.

The experiment uses a 5% random pilot, total labeling budgets of 10%, 20%, and 30%, and 50 replications.

Validation 07 implements the final real-data analysis. Validation 08 generates figures from its audited saved results.

## Model fitting and acquisition

The Gaussian classification model is fitted by semi-supervised maximum likelihood. Observed class memberships are fixed, while memberships of unlabeled observations are latent.

The validated workflow includes checks of EM convergence, parameter stability, curvature precision, and design optimization.

Fisher and adaptive relaxed designs are optimized using Frank–Wolfe with a relative optimality-gap tolerance of 1e-5. Information matrices are inverted using ordinary Cholesky factorization.

Relaxed weights are converted into labeling sets by deterministic top-budget rounding. The criterion is recomputed at the resulting binary design to assess rounding loss.

## Data

Download the original datasets from their providers.

### Statlog Landsat Satellite

Required files:

```text
sat.trn
sat.tst
```

Dataset reference:

Srinivasan, A. (1993). *Statlog (Landsat Satellite)*.
UCI Machine Learning Repository.
DOI: 10.24432/C55887.

### Dry Bean

Required file:

```text
Dry_Bean_Dataset.xlsx
```

Dataset description:

Koklu, M. and Ozkan, I. A. (2020).
Multiclass classification of dry beans using computer vision and machine learning techniques.
*Computers and Electronics in Agriculture*, 174, 105507.
DOI: 10.1016/j.compag.2020.105507.

Place the datasets anywhere beneath the project root, for example:

```text
data/landsat/sat.trn
data/landsat/sat.tst
data/drybean/Dry_Bean_Dataset.xlsx
```

The validation scripts search recursively. Exactly one copy of each required dataset and dependency script must occur beneath the project root.

## Saved results and required caches

The manuscript workflow uses these completed output folders:

```text
realdata_final_38744d8d783f
realdata_converged_387442b8314e
simulation_validated_3874617777a9
realdata_figures_3874327b635d
simulation_figures_387423727afa
```

Obtain the corresponding archives from the repository release assets and extract them without changing their contents.

Dependencies are:

| Workflow | Required saved input |
|---|---|
| Validation 07 | Validation 06 results in `realdata_final_38744d8d783f` |
| Validation 08 | Validation 07 results in `realdata_converged_387442b8314e` |
| Validation 10 | Validation 09 results in `simulation_validated_3874617777a9` |

Preserve all `.rds`, protocol, raw-result, summary, comparison, and diagnostic files. Do not resave the audited `.rds` files: the figure scripts check their exact file hashes.

These archives are required for the frozen workflow; the scripts alone are insufficient.

## R environment

The audited runs used R 4.5.1 on Windows 11. Exact package versions and session information are recorded in the archived logs.

Packages used across the project include:

```text
ggplot2
patchwork
readxl
dplyr
tidyr
readr
scales
```

Base and recommended R packages are also used. Start from a clean R session so that objects in the global environment do not shadow numerical functions.

## Reproducing figures from the audited results

Start R in the repository root after downloading the datasets and extracting the saved-result archives.

```r
project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

source("validation_08_realdata_figures.R")

validation_08_realdata_figures(
  project_root = project_root,
  validation_07_dir = file.path(
    project_root,
    "realdata_converged_387442b8314e"
  )
)

source("validation_10_simulation_figures.R")

validation_10_simulation_figures(
  project_root = project_root,
  validation_09_dir = file.path(
    project_root,
    "simulation_validated_3874617777a9"
  )
)
```

New output folders are created. The archived results are not overwritten.

Validation 08 uses the saved real-data fits and acquisition sets. Validation 10 uses the saved simulation performance records and separately generates the specified true-model geometry illustrations.

## Running the final analysis workflows

Full analysis runs require substantial computation. Figure reproduction from the audited caches does not require rerunning the performance experiments.

### Real-data analysis

Validation 07 requires the archived validation 06 cache.

```r
source("validation_07_realdata_convergence.R")

validation_07_realdata(
  project_root = project_root,
  validation_06_dir = file.path(
    project_root,
    "realdata_final_38744d8d783f"
  ),
  n_rep = 50L
)
```

### Simulation analysis

Explicitly request 100 replications. The function's default is a two-replication diagnostic.

```r
source("validation_09_simulation.R")

validation_09_simulation(
  project_root = project_root,
  n_rep = 100L
)
```

Interrupted analysis runs can be resumed using the `resume_dir` argument and the output directory printed by the relevant script. Resumption requires matching inputs, implementation, protocol, and the checked R environment.

The supplied figure scripts intentionally accept the specific audited caches used for the manuscript. A newly generated cache may have a different file hash and is not automatically interchangeable with the archived cache.

## Interpretation and reproducibility notes

- Methods within each replication share the initial pilot.
- Budgets are constructed independently from the same pilot fit.
- Test outcomes are excluded from preprocessing, fitting, curvature estimation, and acquisition.
- Dry Bean preprocessing is estimated separately within each training split.
- Simulation oracle acquisition uses true-model quantities.
- Acquisition-geometry figures illustrate the specified designs; they are not evidence of general performance superiority.
- Monte Carlo standard errors and paired comparisons are reported in the saved outputs.
- The simulation results do not establish an advantage over Fisher acquisition.
- Real-data performance depends on the dataset and labeling budget.

## Citation

Setoudehtazangi, F. and McLachlan, G. J. (2026).
*Classification-Risk-Optimal Label Acquisition*.
Manuscript.

Publication details will be updated when available.

## License

See the repository's `LICENSE` file for the code license. External datasets remain subject to their providers' terms.

## Contact

Please open a repository issue for questions about the code or reproduction of the reported results.
