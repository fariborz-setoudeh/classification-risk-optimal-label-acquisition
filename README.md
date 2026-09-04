# Classification-Risk-Optimal Label Acquisition

This repository contains the R code and reproducibility materials for the manuscript

**Classification-Risk-Optimal Label Acquisition**

by Fariborz Setoudehtazangi and Geoffrey J. McLachlan.

The paper develops a label-acquisition framework for parametric multiclass classification in which the acquisition criterion is derived from the local geometry of multiclass zero-one classification risk.

## Repository structure

```text
classification-risk-optimal-label-acquisition/
├── README.md
├── LICENSE
├── R/
│   ├── 01_simulation_scenario1.R
│   ├── 02_simulation_scenario2.R
│   ├── 03_landsat_analysis.R
│   ├── 04_landsat_figure.R
│   ├── 05_drybean_analysis.R
│   └── 06_drybean_figures.R
├── data/
│   ├── landsat/
│   └── drybean/
└── results/
```

## Numerical studies

The code reproduces the two simulation scenarios and the two real-data analyses reported in the manuscript and Supplementary Material.

The acquisition methods considered are:

- **Random**
- **Entropy**
- **Margin**
- **Fisher**
- **Adaptive risk-optimal**

The simulation studies additionally include the **oracle risk-optimal** design as a theoretical benchmark.

### Simulation Scenario 1

`R/01_simulation_scenario1.R` reproduces the heterogeneous three-class QDA simulation. The experiment uses a sample size of 1,000, an initial pilot of 60 labels, total labeling budgets of 10%, 20%, and 30%, and 100 Monte Carlo replications.

### Simulation Scenario 2

`R/02_simulation_scenario2.R` reproduces the weak-covariance-heterogeneity scenario used as a near-linear control configuration. The sample size, pilot size, labeling budgets, and number of Monte Carlo replications are the same as in Scenario 1.

## Statlog Landsat Satellite application

The main real-data analysis uses the Statlog (Landsat Satellite) dataset from the UCI Machine Learning Repository.

The supplied benchmark split is retained:

- 4,435 training observations;
- 2,000 test observations;
- six represented classes.

The analysis uses attributes 17--20, corresponding to the four spectral measurements of the central pixel of the 3 × 3 neighborhood.

The experiment uses a 5% random pilot sample, total labeling budgets of 10%, 20%, and 30%, and 50 repeated pilot/acquisition replications. The supplied UCI test set is kept fixed across replications.

`R/03_landsat_analysis.R` reproduces the numerical analysis.

`R/04_landsat_figure.R` reproduces the acquisition visualization. The first two principal components are used only for visualization; model fitting and acquisition are performed using the four-dimensional feature representation.

## Dry Bean application

The supplementary real-data analysis uses the Dry Bean dataset from the UCI Machine Learning Repository.

After removal of exact duplicate observations, each replication uses a stratified 70/30 training/test split. Predictors are standardized using training-set quantities, and PCA is fitted to the training set. The first five principal components are retained.

The experiment uses a 5% random pilot sample, total labeling budgets of 10%, 20%, and 30%, and 50 repeated replications.

`R/05_drybean_analysis.R` reproduces the numerical analysis.

`R/06_drybean_figures.R` reproduces the corresponding figures.

## Model fitting and acquisition

For the real-data analyses, the Gaussian mixture underlying QDA is fitted by semi-supervised maximum likelihood using the EM algorithm. Observed class memberships are held fixed, while memberships of unlabeled observations are treated as latent.

For the Fisher and adaptive risk-optimal methods, the finite-pool relaxed design problem is solved using a Frank--Wolfe algorithm and converted to an exact labeling set by deterministic top-budget rounding.

## Data

The raw datasets are publicly available from the UCI Machine Learning Repository and are not included in this repository.

### Landsat Satellite

Download the original Statlog (Landsat Satellite) files

```text
sat.trn
sat.tst
```

and place them in

```text
data/landsat/
```

Reference:

> Srinivasan, A. (1993). *Statlog (Landsat Satellite)*.  
> UCI Machine Learning Repository.  
> DOI: 10.24432/C55887.

### Dry Bean

Download

```text
Dry_Bean_Dataset.xlsx
```

and place it in

```text
data/drybean/
```

Reference:

> Koklu, M. and Ozkan, I. A. (2020). Multiclass classification of dry beans using computer vision and machine learning techniques.  
> *Computers and Electronics in Agriculture*, 174, 105507.  
> DOI: 10.1016/j.compag.2020.105507.

## R packages

The analyses use the following R packages:

```text
ggplot2
dplyr
tidyr
patchwork
readxl
readr
scales
```

Additional base and recommended R packages are used where required by the individual scripts.

## Reproducing the results

Start R from the repository root directory after placing the datasets in the locations described above.

Run:

```r
source("R/01_simulation_scenario1.R")
source("R/02_simulation_scenario2.R")

source("R/03_landsat_analysis.R")
source("R/04_landsat_figure.R")

source("R/05_drybean_analysis.R")
source("R/06_drybean_figures.R")
```

Fixed random seeds are specified in the scripts.

The full simulation experiments use 100 Monte Carlo replications, while the Landsat and Dry Bean analyses use 50 replications. The simulation and design-optimization scripts can therefore require substantial computation time.

## Reproducibility notes

- Competing acquisition methods within a replication use the same initial pilot sample.
- Each labeling budget is constructed independently from the same pilot fit.
- Test-set information is not used in feature preprocessing or model fitting.
- The supplied Landsat training/test split is retained.
- Dry Bean preprocessing is re-estimated within each training split.
- The oracle risk-optimal rule is used only in the simulation studies and is not an implementable acquisition method.
- The two-dimensional Landsat representation is used only for visualization.

## Citation

If you use the code or methodology, please cite:

> Setoudehtazangi, F. and McLachlan, G. J. (2026).  
> *Classification-Risk-Optimal Label Acquisition*.  
> Manuscript.

The citation will be updated following publication.

## License

The code in this repository is released under the MIT License. The external datasets remain subject to the terms specified by their original providers.

## Contact

For questions about the code or reproducibility materials, please open an issue in this repository.
