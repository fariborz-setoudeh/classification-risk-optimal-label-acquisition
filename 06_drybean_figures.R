
library(ggplot2)
library(dplyr)
library(readr)

res <- read_csv(
  "DryBean_realdata_summary.csv",
  show_col_types = FALSE
)

res <- res %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "Random",
        "Entropy",
        "Margin",
        "Fisher",
        "AdaptiveRisk"
      )
    )
  )

p <- ggplot(
  res,
  aes(
    x = Budget,
    y = MeanError,
    group = Method,
    linetype = Method,
    shape = Method
  )
) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_errorbar(
    aes(
      ymin = MeanError - 1.96 * SEError,
      ymax = MeanError + 1.96 * SEError
    ),
    width = 0.007,
    linewidth = 0.55
  ) +
  scale_x_continuous(
    breaks = c(0.10, 0.20, 0.30),
    labels = c("10%", "20%", "30%")
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.001
    )
  ) +
  labs(
    x = "Total labeling budget",
    y = "Test classification error",
    linetype = "Acquisition rule",
    shape = "Acquisition rule"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  "DryBean_realdata_error.pdf",
  p,
  width = 7.5,
  height = 5.2,
  device = cairo_pdf
)

ggsave(
  "DryBean_realdata_error.png",
  p,
  width = 7.5,
  height = 5.2,
  dpi = 500
)

print(p)
############################################################
# DRY BEAN REAL-DATA STUDY
# Classification-Risk-Optimal Label Acquisition
#
# Fixed analysis protocol (chosen BEFORE comparing methods):
#   - remove exact duplicate rows
#   - repeated stratified 70/30 train/test splits
#   - training-only centering/scaling and PCA
#   - retain first 5 PCs
#   - semi-supervised Gaussian QDA (7-component Gaussian mixture
#     with observed labels fixed and unlabelled labels latent)
#   - 5% random pilot labels
#   - total label budgets: 10%, 20%, 30%
#   - methods: Random, Entropy, Margin, Fisher, Adaptive risk-optimal
#   - each target budget is designed independently from the same pilot fit
#   - Fisher and AdaptiveRisk use a certified Frank--Wolfe solver for the
#     relaxed finite-pool convex design, followed by deterministic top-B rounding
#   - primary outcome: test classification error
#   - secondary outcome: balanced classification error
#
# IMPORTANT:
# Run first with TEST_MODE <- FALSE.
# If the test run completes cleanly, set TEST_MODE <- FALSE.
############################################################

rm(list = ls())
gc()

############################################################
# 0. Packages
############################################################

required_packages <- c("readxl", "ggplot2", "patchwork", "scales")

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if(length(missing_packages) > 0){
  install.packages(missing_packages)
}

library(readxl)
library(ggplot2)
library(patchwork)
library(scales)

############################################################
# 1. User settings
############################################################

DATA_FILE <- "Dry_Bean_Dataset.xlsx"

SEED_MASTER <- 20260901L

TEST_MODE <- TRUE

# In TEST_MODE the certified design solver is allowed up to 600 iterations.
# In the final run it is allowed up to 1200 iterations. Convergence is judged
# by the relative Frank--Wolfe optimality gap, not by iteration count alone.

# Final recommendation:
N_REP_FINAL <- 50L

N_REP <- if(TEST_MODE) 2L else N_REP_FINAL

TRAIN_FRAC <- 0.70

N_PC <- 5L

PILOT_FRAC <- 0.05

BUDGETS <- c(0.10, 0.20, 0.30)

METHODS <- c(
  "Random",
  "Entropy",
  "Margin",
  "Fisher",
  "AdaptiveRisk"
)

# Semi-supervised Gaussian-mixture EM
EM_MAXIT <- if(TEST_MODE) 100L else 300L
EM_TOL <- 1e-7

# Small numerical covariance floor.
# This is numerical stabilization, not a tuning parameter selected
# according to method performance.
COV_FLOOR_REL <- 1e-6

# Information-matrix inverse floor
INFO_EIG_FLOOR_REL <- 1e-8

# Boundary-kernel settings used to estimate H_R
BOUNDARY_BW_MIN <- 0.05
BOUNDARY_BW_MAX <- 1.00
BOUNDARY_KERNEL_CUTOFF <- 4.0

# Require at least this many near-boundary observations for a pair
MIN_BOUNDARY_POINTS <- 10L

# Console progress
VERBOSE <- TRUE

############################################################
# 2. Basic utilities
############################################################

log_sum_exp_rows <- function(M){
  m <- apply(M, 1L, max)
  m + log(rowSums(exp(M - m)))
}

softmax_rows <- function(logM){
  z <- log_sum_exp_rows(logM)
  exp(logM - z)
}

symmetrize <- function(M){
  (M + t(M)) / 2
}

regularize_cov <- function(S, rel_floor = COV_FLOOR_REL){
  S <- symmetrize(S)
  
  ee <- eigen(S, symmetric = TRUE)
  
  base <- max(mean(diag(S)), 1e-10)
  floor_val <- rel_floor * base
  
  vals <- pmax(ee$values, floor_val)
  
  S2 <- ee$vectors %*% (vals * t(ee$vectors))
  symmetrize(S2)
}

safe_inverse <- function(M, rel_floor = INFO_EIG_FLOOR_REL){
  M <- symmetrize(M)
  
  ee <- eigen(M, symmetric = TRUE)
  
  maxeig <- max(ee$values)
  floor_val <- max(rel_floor * max(maxeig, 1), 1e-12)
  
  vals <- pmax(ee$values, floor_val)
  
  out <- ee$vectors %*% ((1 / vals) * t(ee$vectors))
  symmetrize(out)
}

row_quad_bilinear <- function(A, M, B){
  # Returns row-wise A_i' M B_i.
  rowSums((A %*% M) * B)
}

balanced_error <- function(truth, pred, levels_all){
  truth <- factor(truth, levels = levels_all)
  pred  <- factor(pred,  levels = levels_all)
  
  errs <- sapply(levels_all, function(k){
    ind <- truth == k
    if(!any(ind)) return(NA_real_)
    mean(pred[ind] != truth[ind])
  })
  
  mean(errs, na.rm = TRUE)
}

class_error_vector <- function(truth, pred, levels_all){
  truth <- factor(truth, levels = levels_all)
  pred  <- factor(pred,  levels = levels_all)
  
  out <- sapply(levels_all, function(k){
    ind <- truth == k
    if(!any(ind)) return(NA_real_)
    mean(pred[ind] != truth[ind])
  })
  
  names(out) <- levels_all
  out
}

############################################################
# 3. Stratified train/test split
############################################################

stratified_split <- function(y, train_frac, seed){
  set.seed(seed)
  
  lev <- levels(y)
  
  train_idx <- integer(0)
  
  for(k in lev){
    idx <- which(y == k)
    nk_train <- floor(train_frac * length(idx))
    
    train_idx <- c(
      train_idx,
      sample(idx, size = nk_train, replace = FALSE)
    )
  }
  
  train_idx <- sort(train_idx)
  test_idx <- setdiff(seq_along(y), train_idx)
  
  list(
    train = train_idx,
    test = test_idx
  )
}

############################################################
# 4. Training-only standardization and PCA
############################################################

fit_preprocess <- function(X_train, n_pc = N_PC){
  mu <- colMeans(X_train)
  ss <- apply(X_train, 2L, sd)
  
  if(any(!is.finite(ss)) || any(ss <= 0)){
    stop("A training predictor has zero or non-finite standard deviation.")
  }
  
  Xs <- sweep(X_train, 2L, mu, "-")
  Xs <- sweep(Xs, 2L, ss, "/")
  
  pca <- prcomp(
    Xs,
    center = FALSE,
    scale. = FALSE
  )
  
  list(
    center = mu,
    scale = ss,
    rotation = pca$rotation[, seq_len(n_pc), drop = FALSE],
    sdev = pca$sdev
  )
}

apply_preprocess <- function(X, prep){
  Xs <- sweep(X, 2L, prep$center, "-")
  Xs <- sweep(Xs, 2L, prep$scale, "/")
  
  Xs %*% prep$rotation
}

############################################################
# 5. Gaussian density calculations
############################################################

log_dmvnorm_chol <- function(X, mu, Sigma){
  d <- ncol(X)
  
  Sigma <- regularize_cov(Sigma)
  
  R <- chol(Sigma)
  
  XC <- sweep(X, 2L, mu, "-")
  
  # Solve R' z = x' row-by-row through forwardsolve.
  Z <- t(
    forwardsolve(
      t(R),
      t(XC)
    )
  )
  
  mahal <- rowSums(Z^2)
  
  logdet <- 2 * sum(log(diag(R)))
  
  -0.5 * (
    d * log(2 * pi) +
      logdet +
      mahal
  )
}

compute_log_joint <- function(X, fit){
  n <- nrow(X)
  g <- length(fit$pi)
  
  out <- matrix(
    NA_real_,
    nrow = n,
    ncol = g
  )
  
  for(k in seq_len(g)){
    out[, k] <-
      log(fit$pi[k]) +
      log_dmvnorm_chol(
        X,
        fit$mu[k, ],
        fit$Sigma[[k]]
      )
  }
  
  colnames(out) <- fit$classes
  out
}

predict_ss_qda <- function(fit, X){
  lj <- compute_log_joint(X, fit)
  post <- softmax_rows(lj)
  
  cls <- fit$classes[max.col(post, ties.method = "first")]
  
  list(
    class = factor(cls, levels = fit$classes),
    posterior = post,
    log_joint = lj
  )
}

############################################################
# 6. Semi-supervised QDA / Gaussian mixture EM
############################################################

fit_ss_qda <- function(
    X,
    y,
    labelled,
    classes,
    maxit = EM_MAXIT,
    tol = EM_TOL,
    cov_floor_rel = COV_FLOOR_REL
){
  n <- nrow(X)
  d <- ncol(X)
  g <- length(classes)
  
  y <- factor(y, levels = classes)
  
  lab_idx <- which(labelled)
  
  if(length(lab_idx) == 0L){
    stop("No pilot labels supplied.")
  }
  
  pilot_counts <- table(
    factor(y[lab_idx], levels = classes)
  )
  
  if(any(pilot_counts <= d + 1L)){
    stop(
      paste0(
        "Pilot has too few observations for at least one class. Counts: ",
        paste(names(pilot_counts), pilot_counts, collapse = ", "),
        ". Increase PILOT_FRAC only if this recurs systematically."
      )
    )
  }
  
  global_mu <- colMeans(X)
  global_S <- regularize_cov(cov(X), cov_floor_rel)
  
  pi_k <- as.numeric(pilot_counts + 0.5)
  pi_k <- pi_k / sum(pi_k)
  
  mu_k <- matrix(
    NA_real_,
    nrow = g,
    ncol = d
  )
  
  Sigma_k <- vector("list", g)
  
  for(k in seq_len(g)){
    idx <- lab_idx[y[lab_idx] == classes[k]]
    
    mu_k[k, ] <- colMeans(X[idx, , drop = FALSE])
    
    Sk <- cov(X[idx, , drop = FALSE])
    
    if(any(!is.finite(Sk))){
      Sk <- global_S
    }
    
    Sigma_k[[k]] <- regularize_cov(
      Sk,
      cov_floor_rel
    )
  }
  
  old_ll <- -Inf
  
  resp <- matrix(
    0,
    nrow = n,
    ncol = g
  )
  
  for(iter in seq_len(maxit)){
    
    fit_now <- list(
      pi = pi_k,
      mu = mu_k,
      Sigma = Sigma_k,
      classes = classes
    )
    
    log_joint <- compute_log_joint(X, fit_now)
    
    # Unlabelled posterior memberships
    resp[,] <- softmax_rows(log_joint)
    
    # Labelled observations remain fixed to their observed classes
    resp[lab_idx, ] <- 0
    
    lab_class_num <- match(
      as.character(y[lab_idx]),
      classes
    )
    
    resp[cbind(lab_idx, lab_class_num)] <- 1
    
    Nk <- colSums(resp)
    
    if(any(Nk <= d + 1)){
      stop("EM produced an effectively empty class.")
    }
    
    pi_new <- Nk / n
    
    mu_new <- matrix(
      0,
      nrow = g,
      ncol = d
    )
    
    Sigma_new <- vector("list", g)
    
    for(k in seq_len(g)){
      wk <- resp[, k]
      
      mu_new[k, ] <-
        colSums(X * wk) / Nk[k]
      
      XC <- sweep(
        X,
        2L,
        mu_new[k, ],
        "-"
      )
      
      S <- crossprod(
        XC * sqrt(wk)
      ) / Nk[k]
      
      Sigma_new[[k]] <- regularize_cov(
        S,
        cov_floor_rel
      )
    }
    
    # Observed-data log likelihood:
    # labelled observations use their known class;
    # unlabelled observations use the mixture likelihood.
    unl_idx <- which(!labelled)
    
    ll_unl <- if(length(unl_idx) > 0L){
      sum(log_sum_exp_rows(log_joint[unl_idx, , drop = FALSE]))
    } else {
      0
    }
    
    ll_lab <- sum(
      log_joint[
        cbind(
          lab_idx,
          lab_class_num
        )
      ]
    )
    
    ll <- ll_unl + ll_lab
    
    pi_k <- pi_new
    mu_k <- mu_new
    Sigma_k <- Sigma_new
    
    if(iter > 1L){
      rel_change <- abs(ll - old_ll) / (1 + abs(old_ll))
      
      if(rel_change < tol){
        break
      }
    }
    
    old_ll <- ll
  }
  
  fit <- list(
    pi = pi_k,
    mu = mu_k,
    Sigma = Sigma_k,
    classes = classes,
    logLik = ll,
    iterations = iter,
    pilot_counts = pilot_counts
  )
  
  fit
}

############################################################
# 7. Parameter indexing for Gaussian QDA
############################################################

vech_pairs <- function(d){
  out <- which(lower.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
  
  # which() returns column-major order, which is a standard vech ordering.
  out
}

build_param_index <- function(g, d){
  n_alpha <- g - 1L
  n_cov <- d * (d + 1L) / 2L
  block_size <- d + n_cov
  
  alpha_idx <- seq_len(n_alpha)
  
  class_blocks <- vector("list", g)
  mu_blocks <- vector("list", g)
  cov_blocks <- vector("list", g)
  
  pos <- n_alpha + 1L
  
  for(k in seq_len(g)){
    mu_idx <- pos:(pos + d - 1L)
    pos <- pos + d
    
    cov_idx <- pos:(pos + n_cov - 1L)
    pos <- pos + n_cov
    
    mu_blocks[[k]] <- mu_idx
    cov_blocks[[k]] <- cov_idx
    class_blocks[[k]] <- c(mu_idx, cov_idx)
  }
  
  list(
    alpha = alpha_idx,
    mu = mu_blocks,
    cov = cov_blocks,
    class = class_blocks,
    n_alpha = n_alpha,
    n_cov = n_cov,
    class_block_size = block_size,
    p = pos - 1L,
    vech_pairs = vech_pairs(d)
  )
}

############################################################
# 8. Class-score components
############################################################

prior_score_vectors <- function(pi_k){
  g <- length(pi_k)
  A <- matrix(
    0,
    nrow = g,
    ncol = g - 1L
  )
  
  base <- -pi_k[seq_len(g - 1L)]
  
  for(k in seq_len(g)){
    A[k, ] <- base
    
    if(k < g){
      A[k, k] <- A[k, k] + 1
    }
  }
  
  A
}

cov_score_matrix <- function(X, mu, Sigma, pairs){
  n <- nrow(X)
  d <- ncol(X)
  m <- nrow(pairs)
  
  invS <- solve(Sigma)
  
  out <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    q <- X[i, ] - mu
    
    M <-
      invS %*%
      (tcrossprod(q) - Sigma) %*%
      invS
    
    for(j in seq_len(m)){
      a <- pairs[j, 1L]
      b <- pairs[j, 2L]
      
      out[i, j] <- if(a == b){
        0.5 * M[a, a]
      } else {
        M[a, b]
      }
    }
  }
  
  out
}

compute_score_components <- function(X, fit, index){
  n <- nrow(X)
  g <- length(fit$pi)
  d <- ncol(X)
  
  pred <- predict_ss_qda(fit, X)
  
  tau <- pred$posterior
  log_joint <- pred$log_joint
  
  A <- prior_score_vectors(fit$pi)
  
  class_scores <- vector("list", g)
  
  for(k in seq_len(g)){
    invS <- solve(fit$Sigma[[k]])
    
    XC <- sweep(
      X,
      2L,
      fit$mu[k, ],
      "-"
    )
    
    mean_score <- XC %*% invS
    
    cov_score <- cov_score_matrix(
      X,
      fit$mu[k, ],
      fit$Sigma[[k]],
      index$vech_pairs
    )
    
    class_scores[[k]] <- cbind(
      mean_score,
      cov_score
    )
  }
  
  # Marginal feature score \bar t(y).
  bar <- matrix(
    0,
    nrow = n,
    ncol = index$p
  )
  
  # Prior-score part
  bar[, index$alpha] <- tau %*% A
  
  # Class-specific blocks:
  # only t_k has nonzero class-k block, so its posterior mean
  # is tau_k * score_k.
  for(k in seq_len(g)){
    bar[, index$class[[k]]] <-
      class_scores[[k]] * tau[, k]
  }
  
  list(
    tau = tau,
    log_joint = log_joint,
    A = A,
    class_scores = class_scores,
    bar = bar,
    index = index
  )
}

############################################################
# 9. Empirical Fisher-information objects
############################################################

weighted_E_tt <- function(score_obj, weights){
  # Computes:
  # (1/n) sum_i weights_i sum_k tau_ik t_ik t_ik'
  #
  # Exploits sparsity of each class score t_k.
  
  tau <- score_obj$tau
  A <- score_obj$A
  Slist <- score_obj$class_scores
  index <- score_obj$index
  
  n <- nrow(tau)
  g <- ncol(tau)
  P <- index$p
  
  M <- matrix(
    0,
    nrow = P,
    ncol = P
  )
  
  for(k in seq_len(g)){
    wk <- weights * tau[, k]
    
    sw <- sum(wk)
    
    ak <- A[k, ]
    
    # alpha-alpha
    M[index$alpha, index$alpha] <-
      M[index$alpha, index$alpha] +
      sw * tcrossprod(ak)
    
    Sk <- Slist[[k]]
    Bk <- index$class[[k]]
    
    # alpha-class k
    cross_ab <- tcrossprod(
      ak,
      colSums(Sk * wk)
    )
    
    M[index$alpha, Bk] <-
      M[index$alpha, Bk] +
      cross_ab
    
    M[Bk, index$alpha] <-
      M[Bk, index$alpha] +
      t(cross_ab)
    
    # class k - class k
    M[Bk, Bk] <-
      M[Bk, Bk] +
      crossprod(
        Sk * sqrt(wk)
      )
  }
  
  symmetrize(M / n)
}

compute_information_objects <- function(score_obj, labelled){
  n <- nrow(score_obj$tau)
  
  # Feature information
  IY <- crossprod(score_obj$bar) / n
  IY <- symmetrize(IY)
  
  # Complete-classification information:
  # I_CC = E[t_Z t_Z']
  ICC <- weighted_E_tt(
    score_obj,
    rep(1, n)
  )
  
  # Label information retained at currently labelled points:
  w_lab <- as.numeric(labelled)
  
  Ett_lab <- weighted_E_tt(
    score_obj,
    w_lab
  )
  
  if(any(labelled)){
    B <- score_obj$bar[labelled, , drop = FALSE]
    
    barbar_lab <- crossprod(B) / n
  } else {
    barbar_lab <- matrix(
      0,
      nrow = ncol(score_obj$bar),
      ncol = ncol(score_obj$bar)
    )
  }
  
  J_lab <- Ett_lab - barbar_lab
  
  I_current <- IY + J_lab
  
  list(
    IY = symmetrize(IY),
    ICC = symmetrize(ICC),
    J_lab = symmetrize(J_lab),
    I = symmetrize(I_current)
  )
}

############################################################
# 10. Empirical active-face curvature H_R
#
# We use a kernel/coarea approximation of the boundary integral.
#
# For pair k,l, with ell_kl(y)=log r_k(y)-log r_l(y),
#
#   H_kl = integral_{F_kl}
#          r_kl(s) delta_t delta_t' / ||grad ell_kl|| dS
#
# can be represented through a delta kernel in ell_kl.
# Sampling y from the empirical feature distribution yields the
# plug-in weight sqrt(tau_k tau_l) K_h(ell_kl), localized to
# observations for which k and l are the two leading classes.
#
# Each contribution is an outer product, so the estimator is PSD.
############################################################

estimate_HR <- function(score_obj){
  tau <- score_obj$tau
  log_joint <- score_obj$log_joint
  A <- score_obj$A
  Slist <- score_obj$class_scores
  index <- score_obj$index
  
  n <- nrow(tau)
  g <- ncol(tau)
  P <- index$p
  
  top1 <- max.col(
    tau,
    ties.method = "first"
  )
  
  tau2 <- tau
  tau2[cbind(seq_len(n), top1)] <- -Inf
  
  top2 <- max.col(
    tau2,
    ties.method = "first"
  )
  
  H <- matrix(
    0,
    nrow = P,
    ncol = P
  )
  
  pair_diag <- list()
  rr <- 1L
  
  for(k in seq_len(g - 1L)){
    for(l in (k + 1L):g){
      
      active <-
        (top1 == k & top2 == l) |
        (top1 == l & top2 == k)
      
      n_active <- sum(active)
      
      if(n_active < MIN_BOUNDARY_POINTS){
        pair_diag[[rr]] <- data.frame(
          k = k,
          l = l,
          n_active = n_active,
          bandwidth = NA_real_,
          n_kernel = 0L,
          weight_sum = 0
        )
        rr <- rr + 1L
        next
      }
      
      ell <- log_joint[, k] - log_joint[, l]
      
      ell_active <- ell[active]
      
      s_ell <- sd(ell_active)
      
      if(!is.finite(s_ell) || s_ell <= 0){
        s_ell <- mad(ell_active, constant = 1)
      }
      
      h <-
        1.06 *
        s_ell *
        n_active^(-1 / 5)
      
      h <- max(
        BOUNDARY_BW_MIN,
        min(BOUNDARY_BW_MAX, h)
      )
      
      keep <-
        active &
        abs(ell) <= BOUNDARY_KERNEL_CUTOFF * h
      
      idx <- which(keep)
      
      if(length(idx) < MIN_BOUNDARY_POINTS){
        pair_diag[[rr]] <- data.frame(
          k = k,
          l = l,
          n_active = n_active,
          bandwidth = h,
          n_kernel = length(idx),
          weight_sum = 0
        )
        rr <- rr + 1L
        next
      }
      
      kern <-
        dnorm(
          ell[idx] / h
        ) / h
      
      # Symmetric approximation to r_kl / p on the boundary.
      w <-
        sqrt(
          pmax(tau[idx, k], 0) *
            pmax(tau[idx, l], 0)
        ) *
        kern
      
      good <- is.finite(w) & w > 0
      
      idx <- idx[good]
      w <- w[good]
      
      if(length(idx) < MIN_BOUNDARY_POINTS){
        pair_diag[[rr]] <- data.frame(
          k = k,
          l = l,
          n_active = n_active,
          bandwidth = h,
          n_kernel = length(idx),
          weight_sum = sum(w)
        )
        rr <- rr + 1L
        next
      }
      
      D <- matrix(
        0,
        nrow = length(idx),
        ncol = P
      )
      
      # Prior-score difference
      D[, index$alpha] <-
        matrix(
          A[k, ] - A[l, ],
          nrow = length(idx),
          ncol = index$n_alpha,
          byrow = TRUE
        )
      
      # Class k score and class l score
      D[, index$class[[k]]] <-
        Slist[[k]][idx, , drop = FALSE]
      
      D[, index$class[[l]]] <-
        -Slist[[l]][idx, , drop = FALSE]
      
      H <-
        H +
        crossprod(
          D * sqrt(w / n)
        )
      
      pair_diag[[rr]] <- data.frame(
        k = k,
        l = l,
        n_active = n_active,
        bandwidth = h,
        n_kernel = length(idx),
        weight_sum = sum(w)
      )
      
      rr <- rr + 1L
    }
  }
  
  H <- symmetrize(H)
  
  list(
    H = H,
    diagnostics = do.call(rbind, pair_diag)
  )
}

############################################################
# 11. Efficient pointwise value score
#
# For any target matrix L, let
#
#   G = I^{-1} L I^{-1}.
#
# The marginal label value is
#
#   tr{G J(y)}
#   = E[t_Z' G t_Z | Y=y] - \bar t' G \bar t.
#
# Each t_k contains only:
#   - the common prior-logit block, and
#   - the class-k Gaussian block.
#
# We exploit that sparsity rather than constructing J(y) explicitly.
############################################################

classification_value_from_G <- function(
    score_obj,
    G,
    candidate = NULL
){
  tau <- score_obj$tau
  A <- score_obj$A
  Slist <- score_obj$class_scores
  index <- score_obj$index
  
  if(is.null(candidate)){
    candidate <- seq_len(nrow(tau))
  }
  
  tau_c <- tau[candidate, , drop = FALSE]
  
  g <- ncol(tau)
  nc <- length(candidate)
  
  GAA <- G[index$alpha, index$alpha, drop = FALSE]
  
  # First term: sum_k tau_k q_kk
  term1 <- numeric(nc)
  
  # Second term: sum_{k,l} tau_k tau_l q_kl
  term2 <- numeric(nc)
  
  # Cache class score subsets
  Sc <- lapply(
    Slist,
    function(S){
      S[candidate, , drop = FALSE]
    }
  )
  
  # Diagonal q_kk
  qdiag <- vector("list", g)
  
  for(k in seq_len(g)){
    ak <- A[k, ]
    Bk <- index$class[[k]]
    Sk <- Sc[[k]]
    
    const <- as.numeric(
      t(ak) %*% GAA %*% ak
    )
    
    v_left <- as.numeric(
      G[Bk, index$alpha, drop = FALSE] %*% ak
    )
    
    lin <- 2 * as.numeric(
      Sk %*% v_left
    )
    
    bil <- row_quad_bilinear(
      Sk,
      G[Bk, Bk, drop = FALSE],
      Sk
    )
    
    qkk <- const + lin + bil
    
    qdiag[[k]] <- qkk
    
    term1 <- term1 + tau_c[, k] * qkk
    
    term2 <- term2 +
      (tau_c[, k]^2) * qkk
  }
  
  # Off-diagonal q_kl; q_lk=q_kl since G is symmetric.
  for(k in seq_len(g - 1L)){
    for(l in (k + 1L):g){
      ak <- A[k, ]
      al <- A[l, ]
      
      Bk <- index$class[[k]]
      Bl <- index$class[[l]]
      
      Sk <- Sc[[k]]
      Sl <- Sc[[l]]
      
      const <- as.numeric(
        t(ak) %*% GAA %*% al
      )
      
      # a_k' G_{A,Bl} s_l
      v_l <- as.numeric(
        G[Bl, index$alpha, drop = FALSE] %*% ak
      )
      
      part_l <- as.numeric(
        Sl %*% v_l
      )
      
      # s_k' G_{Bk,A} a_l
      v_k <- as.numeric(
        G[Bk, index$alpha, drop = FALSE] %*% al
      )
      
      part_k <- as.numeric(
        Sk %*% v_k
      )
      
      bil <- row_quad_bilinear(
        Sk,
        G[Bk, Bl, drop = FALSE],
        Sl
      )
      
      qkl <- const + part_l + part_k + bil
      
      term2 <- term2 +
        2 *
        tau_c[, k] *
        tau_c[, l] *
        qkl
    }
  }
  
  psi <- term1 - term2
  
  # Numerical roundoff can produce tiny negative values.
  psi[psi < 0 & psi > -1e-8] <- 0
  
  psi
}

############################################################
# 12. Acquisition scores
############################################################

score_entropy <- function(tau){
  -rowSums(
    tau * log(pmax(tau, 1e-15))
  )
}

score_margin <- function(tau){
  # Smaller top-two margin = more uncertain.
  # Return negative margin so that "larger score is better".
  sorted <- t(
    apply(
      tau,
      1L,
      sort,
      decreasing = TRUE
    )
  )
  
  -(sorted[, 1L] - sorted[, 2L])
}

score_fisher <- function(score_obj, info_obj, candidate){
  Iinv <- safe_inverse(info_obj$I)
  
  G <-
    Iinv %*%
    info_obj$ICC %*%
    Iinv
  
  G <- symmetrize(G)
  
  classification_value_from_G(
    score_obj,
    G,
    candidate
  )
}

score_adaptive_risk <- function(
    score_obj,
    info_obj,
    H,
    candidate
){
  Iinv <- safe_inverse(info_obj$I)
  
  G <-
    Iinv %*%
    H %*%
    Iinv
  
  G <- symmetrize(G)
  
  classification_value_from_G(
    score_obj,
    G,
    candidate
  )
}

############################################################
# 13. Exact finite-pool convex design solver
#
# IMPORTANT:
# This implements the criterion used in the paper:
#
#   min_a tr{ L I(a)^(-1) }
#
# subject to 0 <= a_i <= 1 and sum_i a_i = B,
#
# with the pilot labels fixed at acquisition probability 1.
# The fitted model, J(y), I_Y and H_R are estimated ONCE from
# the pilot, and each target budget is solved independently.
############################################################

project_capped_simplex <- function(v, B, tol = 1e-10){
  n <- length(v)
  
  B <- max(0, min(B, n))
  
  if(B <= tol){
    return(rep(0, n))
  }
  
  if(B >= n - tol){
    return(rep(1, n))
  }
  
  lo <- min(v) - 1
  hi <- max(v) + 1
  
  for(iter in seq_len(120L)){
    lambda <- 0.5 * (lo + hi)
    
    a <- pmin(
      1,
      pmax(0, v - lambda)
    )
    
    if(sum(a) > B){
      lo <- lambda
    } else {
      hi <- lambda
    }
  }
  
  pmin(
    1,
    pmax(
      0,
      v - 0.5 * (lo + hi)
    )
  )
}


weighted_J_information <- function(
    score_obj,
    weights_full
){
  n <- nrow(score_obj$tau)
  
  Ett <- weighted_E_tt(
    score_obj,
    weights_full
  )
  
  B <- score_obj$bar
  
  barbar <- crossprod(
    B * sqrt(weights_full)
  ) / n
  
  symmetrize(
    Ett - barbar
  )
}


build_information_from_design <- function(
    score_obj,
    pilot_labelled,
    candidate,
    a_candidate
){
  n <- nrow(score_obj$tau)
  
  w <- numeric(n)
  
  w[pilot_labelled] <- 1
  w[candidate] <- a_candidate
  
  IY <- crossprod(
    score_obj$bar
  ) / n
  
  IY <- symmetrize(IY)
  
  Jw <- weighted_J_information(
    score_obj,
    w
  )
  
  symmetrize(
    IY + Jw
  )
}


design_objective <- function(
    score_obj,
    pilot_labelled,
    candidate,
    a_candidate,
    L
){
  Icur <- build_information_from_design(
    score_obj,
    pilot_labelled,
    candidate,
    a_candidate
  )
  
  Iinv <- safe_inverse(Icur)
  
  sum(
    L * t(Iinv)
  )
}


solve_finite_pool_design <- function(
    score_obj,
    pilot_labelled,
    candidate,
    L,
    B,
    maxit = 1200L,
    tol = 1e-5,
    verbose = FALSE
){
  # Frank--Wolfe solver with an explicit global optimality certificate.
  #
  # Objective:
  #   Phi(a) = tr{ L I(a)^(-1) }
  #
  # Feasible set:
  #   0 <= a_i <= 1, sum_i a_i = B.
  #
  # Since I(a) is affine in a and Phi is convex on the positive-definite
  # cone, the Frank--Wolfe gap is a valid upper bound on suboptimality.
  #
  # The gradient is
  #   d Phi / d a_i = -(1/n) psi_a(Y_i),
  # so the linear minimization oracle puts unit weight on the B candidate
  # points having the largest current marginal values psi_a(Y_i).
  
  N <- length(candidate)
  
  B <- as.integer(
    round(
      min(
        max(B, 0),
        N
      )
    )
  )
  
  if(B == 0L){
    a0 <- rep(0, N)
    
    obj0 <- design_objective(
      score_obj,
      pilot_labelled,
      candidate,
      a0,
      L
    )
    
    return(
      list(
        a = a0,
        objective = obj0,
        iterations = 0L,
        converged = TRUE,
        fw_gap = 0,
        relative_fw_gap = 0,
        last_step = 0
      )
    )
  }
  
  if(B == N){
    a1 <- rep(1, N)
    
    obj1 <- design_objective(
      score_obj,
      pilot_labelled,
      candidate,
      a1,
      L
    )
    
    return(
      list(
        a = a1,
        objective = obj1,
        iterations = 0L,
        converged = TRUE,
        fw_gap = 0,
        relative_fw_gap = 0,
        last_step = 0
      )
    )
  }
  
  n_total <- nrow(
    score_obj$tau
  )
  
  # Interior feasible starting point.
  a <- rep(
    B / N,
    N
  )
  
  Icur <- build_information_from_design(
    score_obj,
    pilot_labelled,
    candidate,
    a
  )
  
  Iinv <- safe_inverse(
    Icur
  )
  
  obj <- sum(
    L * t(Iinv)
  )
  
  converged <- FALSE
  fw_gap <- Inf
  relative_fw_gap <- Inf
  last_step <- NA_real_
  
  for(iter in seq_len(maxit)){
    
    Iinv <- safe_inverse(
      Icur
    )
    
    G <-
      Iinv %*%
      L %*%
      Iinv
    
    G <- symmetrize(
      G
    )
    
    psi <- classification_value_from_G(
      score_obj,
      G,
      candidate
    )
    
    if(any(!is.finite(psi))){
      stop(
        "Non-finite Frank-Wolfe gradient."
      )
    }
    
    # Linear minimization oracle:
    # maximize sum_i s_i psi_i subject to s in the capped simplex.
    ord <- order(
      psi,
      decreasing = TRUE
    )
    
    s <- numeric(
      N
    )
    
    s[
      ord[seq_len(B)]
    ] <- 1
    
    # Exact Frank--Wolfe gap:
    # <a-s, grad Phi(a)>, where grad_i Phi = -psi_i/n_total.
    fw_gap <-
      sum(
        (s - a) * psi
      ) /
      n_total
    
    # Numerical protection only.
    if(
      fw_gap < 0 &&
      fw_gap > -1e-10
    ){
      fw_gap <- 0
    }
    
    if(!is.finite(fw_gap)){
      stop(
        "Non-finite Frank-Wolfe gap."
      )
    }
    
    relative_fw_gap <-
      fw_gap /
      max(
        abs(obj),
        1e-12
      )
    
    if(verbose && (
      iter == 1L ||
      iter %% 25L == 0L ||
      relative_fw_gap <= tol
    )){
      cat(
        "    FW iter=",
        iter,
        " objective=",
        signif(obj, 10),
        " gap=",
        signif(fw_gap, 6),
        " rel_gap=",
        signif(relative_fw_gap, 6),
        "
",
sep = ""
      )
    }
    
    if(
      fw_gap >= 0 &&
      relative_fw_gap <= tol
    ){
      converged <- TRUE
      last_step <- 0
      break
    }
    
    # Information matrix at the Frank--Wolfe vertex s.
    Is <- build_information_from_design(
      score_obj,
      pilot_labelled,
      candidate,
      s
    )
    
    # Along the segment a(gamma)=(1-gamma)a+gamma s,
    # information is exactly affine:
    #   I(gamma)=(1-gamma)Icur+gamma Is.
    line_objective <- function(gamma){
      Ig <-
        (1 - gamma) * Icur +
        gamma * Is
      
      Ig_inv <- safe_inverse(
        Ig
      )
      
      sum(
        L * t(Ig_inv)
      )
    }
    
    ls <- optimize(
      line_objective,
      interval = c(0, 1),
      tol = 1e-8
    )
    
    gamma <- ls$minimum
    new_obj <- ls$objective
    
    # Compare also with both endpoints to avoid a numerical line-search
    # artefact near gamma=0 or gamma=1.
    obj0 <- obj
    obj1 <- line_objective(1)
    
    if(obj0 <= new_obj && obj0 <= obj1){
      gamma <- 0
      new_obj <- obj0
    } else if(obj1 < new_obj && obj1 < obj0){
      gamma <- 1
      new_obj <- obj1
    }
    
    if(
      !is.finite(new_obj) ||
      new_obj > obj + 1e-9 * max(1, abs(obj))
    ){
      stop(
        "Frank-Wolfe line search failed to produce a non-increasing objective."
      )
    }
    
    last_step <- gamma
    
    if(gamma <= 1e-12){
      # If the line search says no movement is useful while the certified
      # gap is still too large, stop and report non-convergence rather than
      # pretending convergence.
      break
    }
    
    a <-
      (1 - gamma) * a +
      gamma * s
    
    # Enforce the feasible set against accumulated floating-point error.
    a <- pmin(
      1,
      pmax(
        0,
        a
      )
    )
    
    # Renormalization should be tiny; projection is used only if required.
    if(
      abs(sum(a) - B) >
      1e-8
    ){
      a <- project_capped_simplex(
        a,
        B
      )
    }
    
    Icur <-
      (1 - gamma) * Icur +
      gamma * Is
    
    Icur <- symmetrize(
      Icur
    )
    
    obj <- new_obj
  }
  
  # Recompute the certificate at the returned design, including when
  # maxit was reached immediately after an update.
  Icur <- build_information_from_design(
    score_obj,
    pilot_labelled,
    candidate,
    a
  )
  
  Iinv <- safe_inverse(
    Icur
  )
  
  obj <- sum(
    L * t(Iinv)
  )
  
  G <-
    Iinv %*%
    L %*%
    Iinv
  
  G <- symmetrize(
    G
  )
  
  psi <- classification_value_from_G(
    score_obj,
    G,
    candidate
  )
  
  ord <- order(
    psi,
    decreasing = TRUE
  )
  
  s <- numeric(
    N
  )
  
  s[
    ord[seq_len(B)]
  ] <- 1
  
  fw_gap <-
    sum(
      (s - a) * psi
    ) /
    n_total
  
  if(
    fw_gap < 0 &&
    fw_gap > -1e-10
  ){
    fw_gap <- 0
  }
  
  relative_fw_gap <-
    fw_gap /
    max(
      abs(obj),
      1e-12
    )
  
  converged <-
    is.finite(relative_fw_gap) &&
    fw_gap >= 0 &&
    relative_fw_gap <= tol
  
  list(
    a = a,
    objective = obj,
    iterations = iter,
    converged = converged,
    fw_gap = fw_gap,
    relative_fw_gap = relative_fw_gap,
    last_step = last_step
  )
}


round_design_topB <- function(
    design,
    score_obj,
    pilot_labelled,
    candidate,
    L,
    B
){
  if(B <= 0){
    selected <- integer(0)
  } else {
    selected <- order(
      design$a,
      decreasing = TRUE
    )[seq_len(B)]
  }
  
  a_round <- numeric(
    length(candidate)
  )
  
  if(length(selected) > 0){
    a_round[selected] <- 1
  }
  
  rounded_objective <- design_objective(
    score_obj,
    pilot_labelled,
    candidate,
    a_round,
    L
  )
  
  relaxed_objective <- design$objective
  
  rounding_loss_pct <-
    100 *
    (
      rounded_objective -
        relaxed_objective
    ) /
    max(
      abs(relaxed_objective),
      1e-12
    )
  
  list(
    selected = selected,
    fractional_count = sum(
      design$a > 1e-8 &
        design$a < 1 - 1e-8
    ),
    fractional_prop = mean(
      design$a > 1e-8 &
        design$a < 1 - 1e-8
    ),
    relaxed_objective = relaxed_objective,
    rounded_objective = rounded_objective,
    rounding_loss_pct = rounding_loss_pct
  )
}

############################################################
# 14. Evaluate fitted classifier
############################################################

evaluate_fit <- function(
    fit,
    X_test_pc,
    y_test,
    classes
){
  pred <- predict_ss_qda(
    fit,
    X_test_pc
  )
  
  err <- mean(
    pred$class != y_test
  )
  
  berr <- balanced_error(
    y_test,
    pred$class,
    classes
  )
  
  cerr <- class_error_vector(
    y_test,
    pred$class,
    classes
  )
  
  list(
    error = err,
    balanced_error = berr,
    class_error = cerr
  )
}

############################################################
# 15. Load and clean Dry Bean data
############################################################

if(!file.exists(DATA_FILE)){
  stop(
    "Cannot find data file: ",
    DATA_FILE
  )
}

dat <- as.data.frame(
  read_excel(DATA_FILE)
)

if(!"Class" %in% names(dat)){
  stop("Column 'Class' was not found.")
}

dat$Class <- factor(dat$Class)

predictor_names <- setdiff(
  names(dat),
  "Class"
)

if(!all(vapply(dat[predictor_names], is.numeric, logical(1)))){
  stop("All predictors must be numeric.")
}

if(anyNA(dat)){
  stop("Unexpected missing values found.")
}

n_before <- nrow(dat)

dup_full <- duplicated(dat)

dat <- dat[!dup_full, , drop = FALSE]
rownames(dat) <- NULL

n_after <- nrow(dat)

cat("\n============================================================\n")
cat("DRY BEAN DATA CLEANING\n")
cat("============================================================\n")
cat("Rows before duplicate removal:", n_before, "\n")
cat("Exact duplicate rows removed:", sum(dup_full), "\n")
cat("Rows retained:", n_after, "\n")

# Check whether identical feature vectors remain with different labels.
dup_X_remaining <- duplicated(
  dat[predictor_names]
) |
  duplicated(
    dat[predictor_names],
    fromLast = TRUE
  )

if(any(dup_X_remaining)){
  warning(
    sum(dup_X_remaining),
    " observations share an identical predictor vector after exact-row ",
    "duplicate removal. Inspect before final analysis."
  )
}

classes <- levels(dat$Class)

cat("\nClass counts after cleaning:\n")
print(table(dat$Class))

############################################################
# 16. Monte Carlo repeated-split experiment
#
# Each replication follows the two-stage procedure in the paper:
#
#   1. draw ONE common random pilot;
#   2. fit the semi-supervised QDA once from that pilot;
#   3. estimate J(y), I_Y, I_CC and H_R from the pilot fit;
#   4. for EACH total budget, construct every acquisition rule
#      from the SAME pilot fit;
#   5. for Fisher and AdaptiveRisk, solve the relaxed finite-pool
#      convex design and then use deterministic top-B rounding;
#   6. refit the classifier after revealing the selected labels.
#
# No method is updated sequentially from the 10% result to the
# 20% or 30% result. This is deliberate and matches the paper's
# two-stage adaptive theory and the simulation implementation.
############################################################



############################################################
# REAL-DATA GEOMETRY / ACQUISITION FIGURE
#
# This script reproduces ONE objectively chosen replication from
# the frozen 50-replication Dry Bean study.
#
# Representative replication:
#   among the 50 final replications, choose the one whose paired
#   10% improvement
#
#       Entropy error - AdaptiveRisk error
#
#   is closest to the median paired improvement.
#
# Hence the displayed replication is NOT selected because it makes
# AdaptiveRisk look unusually good.
#
# IMPORTANT:
#   PC1 and PC2 are used ONLY for visualization.
#   Entropy, classification value, H_R, and label acquisition are
#   all computed from the full 5-dimensional PCA-QDA model.
############################################################

FINAL_RAW_FILE <- "DryBean_realdata_raw.csv"
FIG_BUDGET <- 0.10

if(!file.exists(FINAL_RAW_FILE)){
  stop(
    paste0(
      "Cannot find ", FINAL_RAW_FILE,
      ". Put this visualization script in the same folder as the ",
      "final 50-replication CSV files."
    )
  )
}

final_raw <- read.csv(
  FINAL_RAW_FILE,
  stringsAsFactors = FALSE
)

needed_cols <- c(
  "Replication",
  "Budget",
  "Method",
  "Error"
)

if(!all(needed_cols %in% names(final_raw))){
  stop(
    paste0(
      FINAL_RAW_FILE,
      " does not contain the required columns: ",
      paste(needed_cols, collapse = ", ")
    )
  )
}

ad <- final_raw[
  abs(final_raw$Budget - FIG_BUDGET) < 1e-12 &
    final_raw$Method == "AdaptiveRisk",
  c("Replication", "Error"),
  drop = FALSE
]

names(ad)[2] <- "AdaptiveError"

en <- final_raw[
  abs(final_raw$Budget - FIG_BUDGET) < 1e-12 &
    final_raw$Method == "Entropy",
  c("Replication", "Error"),
  drop = FALSE
]

names(en)[2] <- "EntropyError"

paired_rep <- merge(
  ad,
  en,
  by = "Replication"
)

if(nrow(paired_rep) < 2L){
  stop("Could not construct paired 10% AdaptiveRisk-versus-Entropy results.")
}

paired_rep$Improvement <-
  paired_rep$EntropyError -
  paired_rep$AdaptiveError

median_improvement <- median(
  paired_rep$Improvement
)

paired_rep$DistanceToMedian <-
  abs(
    paired_rep$Improvement -
      median_improvement
  )

paired_rep <- paired_rep[
  order(
    paired_rep$DistanceToMedian,
    paired_rep$Replication
  ),
  ,
  drop = FALSE
]

TARGET_REP <- paired_rep$Replication[1L]

cat("\n============================================================\n")
cat("REPRESENTATIVE REAL-DATA VISUALIZATION\n")
cat("============================================================\n")
cat("Target budget: ", sprintf("%.0f%%", 100 * FIG_BUDGET), "\n", sep = "")
cat("Median paired improvement (Entropy - AdaptiveRisk): ",
    signif(median_improvement, 6), "\n", sep = "")
cat("Selected replication: ", TARGET_REP, "\n", sep = "")
cat("Selected-rep improvement: ",
    signif(paired_rep$Improvement[1L], 6), "\n", sep = "")
cat("AdaptiveRisk test error: ",
    signif(paired_rep$AdaptiveError[1L], 6), "\n", sep = "")
cat("Entropy test error: ",
    signif(paired_rep$EntropyError[1L], 6), "\n", sep = "")

write.csv(
  paired_rep,
  "DryBean_realdata_figure_replication_selection.csv",
  row.names = FALSE
)

############################################################
# 1. Reproduce the selected final-study replication exactly
############################################################

# These objects were created immediately before the Monte Carlo loop in
# the frozen v3 script.  The visualization script stops before that loop,
# so recreate them explicitly here.
X_all <- as.matrix(
  dat[, predictor_names, drop = FALSE]
)

y_all <- dat$Class

rep_id <- TARGET_REP

split_seed <- SEED_MASTER + 10000L * rep_id

sp <- stratified_split(
  y_all,
  TRAIN_FRAC,
  split_seed
)

X_train_raw <- X_all[
  sp$train,
  ,
  drop = FALSE
]

y_train <- factor(
  y_all[sp$train],
  levels = classes
)

X_test_raw <- X_all[
  sp$test,
  ,
  drop = FALSE
]

y_test <- factor(
  y_all[sp$test],
  levels = classes
)

prep <- fit_preprocess(
  X_train_raw,
  N_PC
)

X_train_pc <- apply_preprocess(
  X_train_raw,
  prep
)

X_test_pc <- apply_preprocess(
  X_test_raw,
  prep
)

n_train <- nrow(
  X_train_pc
)

pilot_n <- ceiling(
  PILOT_FRAC * n_train
)

total_n <- ceiling(
  FIG_BUDGET * n_train
)

additional_B <-
  total_n -
  pilot_n

set.seed(
  SEED_MASTER + 20000L * rep_id
)

pilot_idx <- sample(
  seq_len(n_train),
  size = pilot_n,
  replace = FALSE
)

pilot_labelled <- rep(
  FALSE,
  n_train
)

pilot_labelled[pilot_idx] <- TRUE

candidate <- which(
  !pilot_labelled
)

cat("Training n: ", n_train, "\n", sep = "")
cat("Pilot n: ", pilot_n, "\n", sep = "")
cat("Additional labels at 10% budget: ", additional_B, "\n", sep = "")

############################################################
# 2. Pilot fit and full 5-dimensional acquisition quantities
############################################################

pilot_fit <- fit_ss_qda(
  X_train_pc,
  y_train,
  pilot_labelled,
  classes,
  maxit = 300L,
  tol = EM_TOL
)

pilot_pred <- predict_ss_qda(
  pilot_fit,
  X_train_pc
)

entropy_all <- score_entropy(
  pilot_pred$posterior
)

index <- build_param_index(
  g = length(classes),
  d = ncol(X_train_pc)
)

score_obj <- compute_score_components(
  X_train_pc,
  pilot_fit,
  index
)

Hout <- estimate_HR(
  score_obj
)

Hfit <- Hout$H

if(
  !all(is.finite(Hfit)) ||
  sum(abs(Hfit)) <= 0
){
  stop("Estimated H_R is non-finite or numerically zero.")
}

############################################################
# 3. Reproduce Entropy and AdaptiveRisk acquisition at 10%
############################################################

sel_entropy <- candidate[
  order(
    entropy_all[candidate],
    decreasing = TRUE
  )[seq_len(additional_B)]
]

Adaptive_design <- solve_finite_pool_design(
  score_obj = score_obj,
  pilot_labelled = pilot_labelled,
  candidate = candidate,
  L = Hfit,
  B = additional_B,
  maxit = 1200L,
  tol = 1e-5,
  verbose = FALSE
)

if(!isTRUE(Adaptive_design$converged)){
  stop(
    paste0(
      "Adaptive design did not converge. Relative FW gap = ",
      signif(Adaptive_design$relative_fw_gap, 6)
    )
  )
}

Adaptive_round <- round_design_topB(
  design = Adaptive_design,
  score_obj = score_obj,
  pilot_labelled = pilot_labelled,
  candidate = candidate,
  L = Hfit,
  B = additional_B
)

sel_adaptive <- candidate[
  Adaptive_round$selected
]

cat("Adaptive design converged: ", Adaptive_design$converged, "\n", sep = "")
cat("Relative FW gap: ",
    signif(Adaptive_design$relative_fw_gap, 6), "\n", sep = "")
cat("Rounding loss (%): ",
    signif(Adaptive_round$rounding_loss_pct, 6), "\n", sep = "")

############################################################
# 4. Classification value at the fitted adaptive design
#
# This is psi_a(y) evaluated using the information matrix of
# the relaxed adaptive design. It is computed in ALL five PCs.
############################################################

I_adaptive <- build_information_from_design(
  score_obj = score_obj,
  pilot_labelled = pilot_labelled,
  candidate = candidate,
  a_candidate = Adaptive_design$a
)

I_adaptive_inv <- safe_inverse(
  I_adaptive
)

G_adaptive <-
  I_adaptive_inv %*%
  Hfit %*%
  I_adaptive_inv

G_adaptive <- symmetrize(
  G_adaptive
)

psi_all <- classification_value_from_G(
  score_obj = score_obj,
  G = G_adaptive,
  candidate = NULL
)

if(any(!is.finite(psi_all))){
  stop("Non-finite classification values.")
}

# Numerical roundoff protection.
psi_all[psi_all < 0 & psi_all > -1e-8] <- 0

# Relative value used ONLY for the colour scale.
# The acquisition itself used the unscaled psi values.
psi_nonneg <- pmax(
  psi_all,
  0
)

psi_scale <- max(
  psi_nonneg
)

if(!is.finite(psi_scale) || psi_scale <= 0){
  stop("Classification value is numerically zero.")
}

relative_psi <-
  psi_nonneg /
  psi_scale

############################################################
# 5. Observation-level plotting data
############################################################

plot_df <- data.frame(
  Row = seq_len(n_train),
  PC1 = X_train_pc[, 1L],
  PC2 = X_train_pc[, 2L],
  Class = y_train,
  Entropy = entropy_all,
  RelativeValue = relative_psi,
  Pilot = pilot_labelled,
  Candidate = !pilot_labelled,
  EntropySelected = seq_len(n_train) %in% sel_entropy,
  AdaptiveSelected = seq_len(n_train) %in% sel_adaptive
)

write.csv(
  plot_df,
  "DryBean_realdata_geometry_plot_data.csv",
  row.names = FALSE
)

selection_summary <- data.frame(
  Replication = TARGET_REP,
  Budget = FIG_BUDGET,
  PilotN = pilot_n,
  AdditionalBudget = additional_B,
  AdaptiveConverged = Adaptive_design$converged,
  AdaptiveIterations = Adaptive_design$iterations,
  AdaptiveFWGap = Adaptive_design$fw_gap,
  AdaptiveRelativeFWGap = Adaptive_design$relative_fw_gap,
  AdaptiveRoundingLossPct = Adaptive_round$rounding_loss_pct,
  MedianPairedImprovement = median_improvement,
  SelectedRepPairedImprovement = paired_rep$Improvement[1L]
)

write.csv(
  selection_summary,
  "DryBean_realdata_geometry_diagnostics.csv",
  row.names = FALSE
)

############################################################
# 6. Publication-quality four-panel figure
############################################################

# Downweight overplotting while retaining the full training pool.
POINT_ALPHA <- 0.34
POINT_SIZE <- 0.65

p1 <- ggplot(
  plot_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Class
  )
) +
  geom_point(
    alpha = 0.48,
    size = 0.72
  ) +
  labs(
    title = "(a) Dry Bean feature geometry",
    x = expression(PC[1]),
    y = expression(PC[2]),
    colour = "Class"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        size = 2.3
      )
    )
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(
      face = "plain",
      size = 12
    )
  )

p2 <- ggplot(
  plot_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Entropy
  )
) +
  geom_point(
    alpha = 0.60,
    size = POINT_SIZE
  ) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    name = "Entropy"
  ) +
  labs(
    title = "(b) Posterior uncertainty",
    x = expression(PC[1]),
    y = expression(PC[2])
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "plain",
      size = 12
    )
  )

p3 <- ggplot(
  plot_df,
  aes(
    x = PC1,
    y = PC2,
    colour = RelativeValue
  )
) +
  geom_point(
    alpha = 0.62,
    size = POINT_SIZE
  ) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    name = expression("Relative " * psi[a](y))
  ) +
  labs(
    title = "(c) Classification value",
    x = expression(PC[1]),
    y = expression(PC[2])
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "plain",
      size = 12
    )
  )

selection_df <- rbind(
  data.frame(
    PC1 = plot_df$PC1,
    PC2 = plot_df$PC2,
    Selected = plot_df$AdaptiveSelected,
    Rule = "Classification-risk optimal"
  ),
  data.frame(
    PC1 = plot_df$PC1,
    PC2 = plot_df$PC2,
    Selected = plot_df$EntropySelected,
    Rule = "Entropy"
  )
)

selection_df$Rule <- factor(
  selection_df$Rule,
  levels = c(
    "Classification-risk optimal",
    "Entropy"
  )
)

p4 <- ggplot(
  selection_df,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    colour = "grey78",
    alpha = 0.32,
    size = 0.55
  ) +
  geom_point(
    data = selection_df[
      selection_df$Selected,
      ,
      drop = FALSE
    ],
    colour = "black",
    alpha = 0.78,
    size = 0.85
  ) +
  facet_wrap(
    ~ Rule,
    nrow = 1
  ) +
  labs(
    title = "(d) Additional labels selected under a 10% total budget",
    x = expression(PC[1]),
    y = expression(PC[2])
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "grey92",
      colour = "grey50"
    ),
    strip.text = element_text(
      size = 10.5
    ),
    plot.title = element_text(
      face = "plain",
      size = 12
    )
  )

final_figure <-
  (
    p1 |
      p2
  ) /
  (
    p3 |
      p4
  ) +
  plot_annotation(
    title =
      "Posterior uncertainty and classification value target different observations",
    subtitle = paste0(
      "Dry Bean data; representative replication ",
      TARGET_REP,
      " selected by the median paired 10% AdaptiveRisk-versus-Entropy improvement"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10
      )
    )
  )

ggsave(
  "DryBean_realdata_geometry_acquisition.pdf",
  final_figure,
  width = 13.5,
  height = 9.5,
  device = cairo_pdf
)

ggsave(
  "DryBean_realdata_geometry_acquisition.png",
  final_figure,
  width = 13.5,
  height = 9.5,
  dpi = 500
)

print(
  final_figure
)

cat("\n============================================================\n")
cat("FILES WRITTEN\n")
cat("============================================================\n")
cat(
  paste(
    c(
      "DryBean_realdata_geometry_acquisition.pdf",
      "DryBean_realdata_geometry_acquisition.png",
      "DryBean_realdata_geometry_plot_data.csv",
      "DryBean_realdata_geometry_diagnostics.csv",
      "DryBean_realdata_figure_replication_selection.csv"
    ),
    collapse = "\n"
  ),
  "\n"
)

cat("\nDONE.\n")
