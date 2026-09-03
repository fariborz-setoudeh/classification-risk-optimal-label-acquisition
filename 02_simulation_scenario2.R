# ============================================================
# SCENARIO 2
# Weak covariance heterogeneity / near-LDA geometry
#
# PURPOSE:
# A control regime in which the three covariance matrices are
# similar. The Bayes boundaries are therefore much closer to
# the homogeneous-covariance geometry than in Scenario 1.
#
# We keep:
#   g = 3
#   p = 2
#   n = 1000
#   pilot_n = 60
#   budgets = 10%, 20%, 30%
#   R = 100
#
# All statistical and computational procedures are identical
# to Scenario 1.
# ============================================================



# ============================================================
# 2. CREATE SEPARATE OUTPUT DIRECTORY
# ============================================================

scenario2_dir <- file.path(
  old_work_dir,
  "Scenario2_near_LDA"
)

if (!dir.exists(scenario2_dir)) {
  dir.create(
    scenario2_dir,
    recursive = TRUE
  )
}

setwd(scenario2_dir)

cat(
  "\nScenario 2 output directory:\n",
  getwd(),
  "\n"
)


# ============================================================
# 3. SCENARIO 2 DATA-GENERATING MODEL
#
# Means and class proportions are kept identical to Scenario 1.
# Only covariance heterogeneity is substantially reduced.
#
# This isolates the role of decision-boundary geometry.
# ============================================================

g <- 3
p <- 2

pi0_s2 <- c(
  0.36,
  0.34,
  0.30
)

mu0_s2 <- rbind(
  c(-1.60, 0.00),
  c( 1.50, 0.20),
  c( 0.00, 1.90)
)


# ------------------------------------------------------------
# Near-common covariance structure
#
# Common reference:
#
#       [ 1.00   0.15 ]
#       [ 0.15   0.90 ]
#
# The three matrices differ only moderately around this
# reference. They remain positive definite.
# ------------------------------------------------------------

Sigma0_s2 <- array(
  NA_real_,
  dim = c(
    p,
    p,
    g
  )
)


Sigma0_s2[, , 1] <- matrix(
  c(
    1.00, 0.15,
    0.15, 0.90
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE
)


Sigma0_s2[, , 2] <- matrix(
  c(
    0.95, 0.12,
    0.12, 0.95
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE
)


Sigma0_s2[, , 3] <- matrix(
  c(
    1.05, 0.18,
    0.18, 0.85
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE
)


true_par_s2 <- list(
  pi = pi0_s2,
  mu = mu0_s2,
  Sigma = Sigma0_s2
)


# ============================================================
# 4. BASIC VALIDITY CHECKS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCENARIO 2: covariance diagnostics\n"
)

cat(
  "============================================\n"
)


for (k in seq_len(g)) {
  
  eig_k <-
    eigen(
      Sigma0_s2[, , k],
      symmetric = TRUE,
      only.values = TRUE
    )$values
  
  cat(
    "Class",
    k,
    ": eigenvalues =",
    paste(
      round(
        eig_k,
        6
      ),
      collapse = ", "
    ),
    "\n"
  )
  
  if (min(eig_k) <= 0) {
    
    stop(
      paste(
        "Scenario 2 covariance matrix",
        k,
        "is not positive definite."
      )
    )
  }
}


# ============================================================
# 5. MAKE SCENARIO 2 THE ACTIVE TRUE MODEL
#
# The validated functions use the global object true_par.
# We therefore replace only this object; no statistical
# function is changed.
# ============================================================

true_par <- true_par_s2


# ============================================================
# 6. POPULATION CHECK
#
# Estimate Bayes error independently before the simulation.
# This is only a diagnostic and does not affect the experiment.
# ============================================================

set.seed(20260902)

population_check <-
  simulate_qda(
    250000,
    true_par
  )


population_bayes_pred <-
  predict_qda(
    population_check$Y,
    true_par
  )


population_bayes_error <-
  mean(
    population_bayes_pred !=
      population_check$z
  )


cat(
  "\nApproximate Scenario 2 Bayes error =",
  round(
    population_bayes_error,
    5
  ),
  "\n"
)


# ============================================================
# 7. POPULATION GEOMETRY
#
# This produces Figure1_geometry.png INSIDE the Scenario 2
# directory, so the Scenario 1 figure is not overwritten.
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCENARIO 2: population geometry\n"
)

cat(
  "============================================\n"
)


set.seed(20260903)

fig_s2_geometry <-
  make_geometry_figure(
    true_par,
    rho = 0.20,
    Npool = 5000,
    ngrid = 180
  )


print(
  fig_s2_geometry
)


# Rename the automatically generated file so that its meaning
# is explicit within the Scenario 2 folder.

if (
  file.exists(
    "Figure1_geometry.png"
  )
) {
  
  file.rename(
    "Figure1_geometry.png",
    "Scenario2_geometry.png"
  )
}


# ============================================================
# 8. MONTE CARLO STUDY
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCENARIO 2: Monte Carlo study\n"
)

cat(
  "============================================\n"
)


# ------------------------------------------------------------
# IMPORTANT
#
# design solver. Therefore we can directly use R = 100 here.
#
# If you prefer one final safety test, replace 100 by 10 once.
# ------------------------------------------------------------

set.seed(20260904)


RES_S2 <-
  run_simulation(
    R = 100,
    n = 1000,
    pilot_n = 60,
    budgets = c(
      0.10,
      0.20,
      0.30
    ),
    test_n = 100000,
    H_ngrid = 150
  )


# ============================================================
# 9. PRINT FULL RESULTS
# ============================================================

cat(
  "\n\n============================================\n"
)

cat(
  "SCENARIO 2: CLASSIFICATION RESULTS\n"
)

cat(
  "============================================\n"
)


print(
  RES_S2$summary,
  width = Inf
)


cat(
  "\n\n============================================\n"
)

cat(
  "SCENARIO 2: DESIGN ROUNDING DIAGNOSTICS\n"
)

cat(
  "============================================\n"
)


print(
  RES_S2$design_diagnostics_summary,
  width = Inf
)


# ============================================================
# 10. RENAME OUTPUT FILES
#
# The functions retain their validated Scenario 1 filenames.
# Since we are already inside a separate directory this is not
# technically necessary, but explicit names are safer.
# ============================================================

rename_if_exists <- function(
    old,
    new) {
  
  if (file.exists(old)) {
    
    if (file.exists(new)) {
      file.remove(new)
    }
    
    file.rename(
      old,
      new
    )
  }
}


rename_if_exists(
  "Study1_raw.csv",
  "Scenario2_raw.csv"
)

rename_if_exists(
  "Study1_raw_PROGRESS.csv",
  "Scenario2_raw_PROGRESS.csv"
)

rename_if_exists(
  "Study1_summary.csv",
  "Scenario2_summary.csv"
)

rename_if_exists(
  "Study1_design_diagnostics.csv",
  "Scenario2_design_diagnostics.csv"
)

rename_if_exists(
  "Study1_design_diagnostics_PROGRESS.csv",
  "Scenario2_design_diagnostics_PROGRESS.csv"
)

rename_if_exists(
  "Study1_design_diagnostics_summary.csv",
  "Scenario2_design_diagnostics_summary.csv"
)

rename_if_exists(
  "Figure2_performance.png",
  "Scenario2_performance.png"
)


# ============================================================
# 11. SAVE MODEL PARAMETERS
# ============================================================

scenario2_parameters <-
  data.frame(
    class = 1:g,
    pi = true_par$pi,
    mu1 = true_par$mu[, 1],
    mu2 = true_par$mu[, 2],
    sigma11 = sapply(
      1:g,
      function(k)
        true_par$Sigma[1, 1, k]
    ),
    sigma12 = sapply(
      1:g,
      function(k)
        true_par$Sigma[1, 2, k]
    ),
    sigma22 = sapply(
      1:g,
      function(k)
        true_par$Sigma[2, 2, k]
    )
  )


write.csv(
  scenario2_parameters,
  "Scenario2_parameters.csv",
  row.names = FALSE
)


# ============================================================
# 12. SAVE A SHORT RUN SUMMARY
# ============================================================

scenario2_run_info <-
  data.frame(
    scenario =
      "Near-LDA / weak covariance heterogeneity",
    n =
      1000,
    pilot_n =
      60,
    replications =
      100,
    test_n =
      100000,
    approximate_population_bayes_error =
      population_bayes_error
  )


write.csv(
  scenario2_run_info,
  "Scenario2_run_info.csv",
  row.names = FALSE
)


# ============================================================
# 13. FINISHED
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "SCENARIO 2 COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "\nFiles saved in:\n",
  getwd(),
  "\n\n"
)

cat(
  "Main files:\n",
  "  Scenario2_geometry.png\n",
  "  Scenario2_performance.png\n",
  "  Scenario2_summary.csv\n",
  "  Scenario2_raw.csv\n",
  "  Scenario2_design_diagnostics_summary.csv\n",
  "  Scenario2_parameters.csv\n",
  "  Scenario2_run_info.csv\n"
)


# ============================================================
# FINAL COMBINED PERFORMANCE FIGURE
# Scenario 1 + Scenario 2
# ============================================================

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# Read results
# ------------------------------------------------------------

S1 <- read.csv(
  "Study1_summary.csv"
)

S2 <- read.csv(
  "Scenario2_near_LDA/Scenario2_summary.csv"
)

S1$Scenario <- "(a) Heterogeneous QDA"
S2$Scenario <- "(b) Weak covariance heterogeneity"

dat <- bind_rows(
  S1,
  S2
)

# ------------------------------------------------------------
# Method ordering
# ------------------------------------------------------------

dat$method <- factor(
  dat$method,
  levels = c(
    "Random",
    "Entropy",
    "Margin",
    "Fisher",
    "Adaptive risk-optimal",
    "Oracle risk-optimal"
  )
)

dat$Scenario <- factor(
  dat$Scenario,
  levels = c(
    "(a) Heterogeneous QDA",
    "(b) Weak covariance heterogeneity"
  )
)

# ------------------------------------------------------------
# Bayes errors
# ------------------------------------------------------------

bayes_df <- data.frame(
  Scenario = factor(
    c(
      "(a) Heterogeneous QDA",
      "(b) Weak covariance heterogeneity"
    ),
    levels = levels(dat$Scenario)
  ),
  bayes_error = c(
    0.14845,
    0.15882
  )
)

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

fig3 <- ggplot(
  dat,
  aes(
    x = rho,
    y = mean_error,
    group = method,
    linetype = method,
    shape = method
  )
) +
  
  geom_hline(
    data = bayes_df,
    aes(
      yintercept = bayes_error
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_error - 1.96 * se_error,
      ymax = mean_error + 1.96 * se_error
    ),
    width = 0.008,
    linewidth = 0.5
  ) +
  
  facet_wrap(
    ~ Scenario,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = c(
      0.10,
      0.20,
      0.30
    ),
    labels = c(
      "10%",
      "20%",
      "30%"
    )
  ) +
  
  labs(
    x = "Proportion of acquired labels",
    y = "Test classification error",
    linetype = "Acquisition rule",
    shape = "Acquisition rule"
  ) +
  
  theme_bw(
    base_size = 13
  ) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    )
  )

# ------------------------------------------------------------
# Show
# ------------------------------------------------------------

print(fig3)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  "Figure3_performance_combined.png",
  fig3,
  width = 12,
  height = 6.5,
  dpi = 400
)

ggsave(
  "Figure3_performance_combined.pdf",
  fig3,
  width = 12,
  height = 6.5
)