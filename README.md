# Classification-Risk-Optimal Label Acquisition

This repository contains the R code and reproducibility materials accompanying the manuscript

**Classification-Risk-Optimal Label Acquisition**.

The paper develops a label-acquisition framework for parametric multiclass semi-supervised classification in which labels are selected according to their contribution to **classification risk**, rather than posterior uncertainty or global Fisher information alone.

The proposed acquisition criterion combines the conditional information supplied by a label with the local geometry of multiclass zero-one classification risk. The numerical studies compare the resulting classification-risk-optimal rule with random sampling, posterior-uncertainty rules, and a Fisher-information-based design.

## Repository structure

```text
classification-risk-optimal-label-acquisition/
│
├── README.md
├── LICENSE
│
├── R/
│   ├── 01_simulation_scenario1.R
│   ├── 02_simulation_scenario2.R
│   ├── 03_landsat_analysis.R
│   ├── 04_landsat_figure.R
│   ├── 05_drybean_analysis.R
│   └── 06_drybean_figures.R
│
├── data/
│   ├── landsat/
│   └── drybean/
│
└── results/
```

If the scripts are stored in the repository root rather than in an `R/` directory, the same execution order applies.

## Methods compared

The simulation and real-data studies consider the following acquisition strategies:

- **Random**: labels are acquired uniformly at random.
- **Entropy**: observations with the largest posterior entropy are selected.
- **Margin**: observations with the smallest difference between the two largest posterior class probabilities are selected.
- **Fisher**: labels are selected by minimizing a complete-classification-information criterion.
- **Adaptive risk-optimal**: labels are selected using the estimated classification-risk criterion developed in the paper.

The simulation studies additionally include an **oracle risk-optimal** design, which evaluates the acquisition criterion using the true data-generating parameters and serves only as a theoretical benchmark.

## Simulation studies

### Scenario 1: heterogeneous QDA

`01_simulation_scenario1.R` implements the main three-class Gaussian QDA simulation with heterogeneous class covariance matrices and curved Bayes decision boundaries.

The default final experiment uses:

- three classes;
- two-dimensional features;
- sample size \(n=1000\);
- an initial pilot of 60 labels;
- total labeling budgets of 10%, 20%, and 30%;
- 100 Monte Carlo replications;
- a large independent test sample for estimating classification error.

The script also evaluates the geometry of posterior uncertainty and classification value and produces the corresponding simulation figures and diagnostic files.

### Scenario 2: weak covariance heterogeneity

`02_simulation_scenario2.R` considers a control regime with covariance matrices that are much more similar across classes. The resulting Bayes geometry is therefore closer to the homogeneous-covariance or LDA setting.

The class proportions, means, sample size, pilot size, labeling budgets, and number of Monte Carlo replications are kept the same as in Scenario 1. This scenario is used to examine how the relative performance of the acquisition rules changes as covariance heterogeneity and decision-boundary curvature are reduced.

## Real-data application: Statlog Landsat Satellite

The main real-data experiment uses the **Statlog (Landsat Satellite)** data from the UCI Machine Learning Repository.

The supplied benchmark split is retained:

- 4,435 training observations;
- 2,000 test observations;
- six represented land-cover classes;
- original class codes 1, 2, 3, 4, 5, and 7.

Following the analysis protocol in the manuscript, only attributes 17--20 are used. These correspond to the four spectral measurements of the central pixel of the \(3\times3\) neighborhood.

The classification model is a six-class semi-supervised Gaussian QDA model. The benchmark training features are observed for all observations, while class labels are revealed according to the acquisition design.

The experiment uses:

- a 5% random pilot sample;
- total labeling budgets of 10%, 20%, and 30%;
- 50 repeated pilot/acquisition replications;
- the fixed UCI benchmark test set for evaluation.

The primary outcome is test classification error. Balanced classification error is also reported.

`03_landsat_analysis.R` reproduces the numerical analysis.

`04_landsat_figure.R` reproduces the geometric/acquisition visualization. The displayed two-dimensional coordinates are the first two principal components and are used **only for visualization**; model fitting, posterior probabilities, classification-value calculations, and acquisition are performed using the full four-dimensional feature representation.

The representative replication used in the figure is selected objectively: at the 20% total labeling budget, it is the replication whose paired Adaptive-risk-versus-Entropy difference is closest to the median difference over the 50 repetitions.

## Supplementary real-data application: Dry Bean

The supplementary real-data experiment uses the **Dry Bean Dataset** from the UCI Machine Learning Repository.

The analysis protocol is:

- remove exact duplicate rows;
- use repeated stratified 70/30 training/test splits;
- standardize predictors using training-set quantities only;
- perform PCA using the training set only;
- retain the first five principal components;
- fit a seven-class semi-supervised Gaussian QDA model;
- use a 5% random pilot sample;
- consider total labeling budgets of 10%, 20%, and 30%;
- repeat the experiment 50 times.

The primary outcome is test classification error and the secondary outcome is balanced classification error.

`05_drybean_analysis.R` reproduces the analysis.

`06_drybean_figures.R` reproduces the corresponding figures from the saved numerical summaries.

## Semi-supervised model fitting

For the real-data analyses, the Gaussian mixture underlying QDA is fitted by maximum likelihood using the EM algorithm.

Observed class memberships are kept fixed, while the class memberships of unlabeled observations are treated as latent. Posterior class probabilities for unlabeled observations are updated in the E-step, and class proportions, means, and covariance matrices are updated in the M-step.

Small covariance regularization terms are used only for numerical stability.

## Finite-pool design optimization

The Fisher and adaptive classification-risk criteria are optimized over the finite unlabeled pool.

For the real-data studies, the relaxed convex design problem is solved using a Frank--Wolfe algorithm. Convergence is assessed using the relative Frank--Wolfe optimality gap.

The continuous design is subsequently converted to an exact labeling set using deterministic top-\(B\) rounding. The scripts save convergence and rounding diagnostics so that the numerical accuracy of the design approximation can be inspected.

## Data

The datasets are publicly available from the UCI Machine Learning Repository.

### Landsat Satellite

Download the Statlog (Landsat Satellite) data and place

```text
sat.trn
sat.tst
```

in

```text
data/landsat/
```

The analysis expects the original UCI files with 36 predictors followed by the class label.

Reference:

> Srinivasan, A. (1993). Statlog (Landsat Satellite).  
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
> *Computers and Electronics in Agriculture*, **174**, 105507.  
> https://doi.org/10.1016/j.compag.2020.105507

The raw datasets are not required to be redistributed through this repository; they can be obtained directly from their public source.

## R requirements

The analyses were developed in R.

Packages used across the scripts include:

```text
ggplot2
dplyr
tidyr
patchwork
readxl
readr
scales
```

Missing packages are installed automatically by some scripts. For a fully controlled environment, they may instead be installed beforehand.

## Reproducing the analyses

Clone or download the repository and start R from the **repository root directory**.

Run the scripts in the following order:

```r
source("R/01_simulation_scenario1.R")
source("R/02_simulation_scenario2.R")

source("R/03_landsat_analysis.R")
source("R/04_landsat_figure.R")

source("R/05_drybean_analysis.R")
source("R/06_drybean_figures.R")
```

The simulation studies are computationally more demanding than the plotting scripts.

The full experiments use 100 Monte Carlo replications for each simulation scenario and 50 repeated experiments for each real-data application. Consequently, reproducing all results may require substantial computation time.

## Random seeds

Fixed random seeds are specified in the scripts to make the simulation studies, pilot selections, and acquisition experiments reproducible.

The principal master seeds are chosen separately for the simulation and real-data analyses.

## Main outputs

Depending on the script, the code produces:

- raw replication-level results;
- summary tables;
- paired comparisons between AdaptiveRisk and competing acquisition rules;
- class-specific error summaries;
- Frank--Wolfe convergence diagnostics;
- rounding diagnostics;
- Bayes-geometry figures;
- classification-performance figures;
- Landsat acquisition-geometry figures.

The main Landsat analysis produces files including:

```text
Landsat_raw.csv
Landsat_summary.csv
Landsat_paired_vs_AdaptiveRisk.csv
Landsat_design_diagnostics.csv
Landsat_design_diagnostics_summary.csv
```

The Dry Bean analysis produces files including:

```text
DryBean_realdata_raw.csv
DryBean_realdata_summary.csv
DryBean_realdata_paired_comparisons.csv
DryBean_realdata_design_diagnostics.csv
```

The simulation scripts analogously save raw results, summaries, figures, and design diagnostics.

## Reproducibility notes

Several details are fixed deliberately to match the experiments reported in the manuscript:

1. All competing acquisition rules within a replication use the same initial pilot sample.
2. Each labeling budget is constructed independently from the same pilot fit rather than sequentially extending the smaller budget.
3. Feature preprocessing does not use test-set information.
4. In the Landsat study, the official UCI training/test split is kept fixed.
5. In the Dry Bean study, preprocessing is re-estimated separately within every training split.
6. The oracle acquisition rule is used only in simulation and is not an implementable method.
7. The two-dimensional Landsat PCA representation is used only for visualization.

These choices are intended to keep the numerical comparisons aligned with the theoretical two-stage acquisition framework developed in the manuscript.

## Citation

If you use this code or methodology, please cite the accompanying manuscript:

> Setoudehtazangi, F.  
> *Classification-Risk-Optimal Label Acquisition*.  
> Manuscript, 2026.

A complete journal citation will be added after publication.

## Code availability

The repository is available at:

**https://github.com/fariborz-setoudeh/classification-risk-optimal-label-acquisition**

## License

The code in this repository is released under the MIT License. The external datasets remain subject to the terms and licenses specified by their original providers.

## Contact

For questions about the code or reproducibility materials, please open an issue in this repository.
