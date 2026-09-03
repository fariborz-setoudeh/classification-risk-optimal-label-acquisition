############################################################
# STATLOG LANDSAT SATELLITE REAL-DATA STUDY
# Classification-Risk-Optimal Label Acquisition
#
# FINAL PRE-SPECIFIED PROTOCOL
#
# Dataset:
#   sat.trn : 4435 benchmark training observations
#   sat.tst : 2000 benchmark test observations
#
# IMPORTANT:
#   The supplied UCI train/test split is kept fixed.
#   We do NOT perform cross-validation or create new train/test splits.
#
# Features:
#   Use ONLY attributes 17--20: the four spectral measurements
#   of the central pixel in each 3x3 neighbourhood.
#   This choice is explicitly suggested in the dataset documentation
#   and avoids the neighbourhood-straddling-boundary issue.
#
# Model:
#   semi-supervised 6-class Gaussian QDA
#   classes = 1,2,3,4,5,7
#   no class 6 observations are present
#
# Repeated experiment:
#   fixed benchmark training/test sets
#   50 independent random pilot/acquisition replications
#   5% random pilot labels from the benchmark training set
#   total labeling budgets: 10%, 20%, 30%
#   each target budget is designed independently from the same pilot fit
#
# Methods:
#   Random
#   Entropy
#   Margin
#   Fisher
#   AdaptiveRisk
#
# Fisher and AdaptiveRisk:
#   relaxed finite-pool convex design solved by Frank--Wolfe,
#   followed by deterministic top-B rounding.
#
# Outcomes:
#   primary   = benchmark test classification error
#   secondary = balanced classification error
#
# Monte Carlo uncertainty:
#   variability across repeated pilot/acquisition replications,
#   conditional on the fixed supplied benchmark train/test split.
#
# BEFORE FINAL RUN:
#   1. Run with TEST_MODE <- TRUE.
#   2. If clean, set TEST_MODE <- FALSE.
############################################################

rm(list = ls())
gc()

############################################################
# 0. Packages
############################################################

required_packages <- c("ggplot2", "dplyr")

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if(length(missing_packages) > 0){
  install.packages(missing_packages)
}

library(ggplot2)
library(dplyr)

############################################################
# 1. User settings
############################################################



TRAIN_FILE <- file.path("data", "landsat", "sat.trn")
TEST_FILE  <- file.path("data", "landsat", "sat.tst")

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR)

SEED_MASTER <- 20260902L
TEST_MODE <- FALSE

SEED_MASTER <- 20260902L

TEST_MODE <- FALSE

N_REP_FINAL <- 50L
N_REP <- if(TEST_MODE) 2L else N_REP_FINAL

# Central-pixel spectral bands, using 1-based attribute numbers.
FEATURE_COLS <- 17:20

PILOT_FRAC <- 0.05
BUDGETS <- c(0.10, 0.20, 0.30)

METHODS <- c(
  "Random",
  "Entropy",
  "Margin",
  "Fisher",
  "AdaptiveRisk"
)

EM_MAXIT <- if(TEST_MODE) 150L else 400L
EM_TOL <- 1e-7

COV_FLOOR_REL <- 1e-6
INFO_EIG_FLOOR_REL <- 1e-8

# Boundary-kernel settings for H_R.
BOUNDARY_BW_MIN <- 0.05
BOUNDARY_BW_MAX <- 1.00
BOUNDARY_KERNEL_CUTOFF <- 4.0
MIN_BOUNDARY_POINTS <- 10L

FW_MAXIT <- if(TEST_MODE) 400L else 1000L
FW_TOL <- 1e-5

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
# 3. Gaussian density calculations
############################################################

log_dmvnorm_chol <- function(X, mu, Sigma){
  d <- ncol(X)
  
  Sigma <- regularize_cov(Sigma)
  
  R <- chol(Sigma)
  
  XC <- sweep(X, 2L, mu, "-")
  
  Z <- t(
    forwardsolve(
      t(R),
      t(XC)
    )
  )
  
  mahal <- rowSums(Z^2)
  logdet <- 2 * sum(log(diag(R))
                    
  )
  
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
# 4. Semi-supervised QDA / Gaussian mixture EM
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
        "Too few labelled observations for at least one class. Counts: ",
        paste(names(pilot_counts), pilot_counts, collapse = ", ")
      )
    )
  }
  
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
    
    resp[,] <- softmax_rows(log_joint)
    
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
      
      mu_new[k, ] <- colSums(X * wk) / Nk[k]
      
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
  
  list(
    pi = pi_k,
    mu = mu_k,
    Sigma = Sigma_k,
    classes = classes,
    logLik = ll,
    iterations = iter,
    labelled_counts = pilot_counts
  )
}

############################################################
# 5. Gaussian-QDA parameter indexing
############################################################

vech_pairs <- function(d){
  which(
    lower.tri(matrix(0, d, d), diag = TRUE),
    arr.ind = TRUE
  )
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
# 6. Class-score components
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
  
  bar <- matrix(
    0,
    nrow = n,
    ncol = index$p
  )
  
  bar[, index$alpha] <- tau %*% A
  
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
# 7. Empirical information objects
############################################################

weighted_E_tt <- function(score_obj, weights){
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
    
    M[index$alpha, index$alpha] <-
      M[index$alpha, index$alpha] +
      sw * tcrossprod(ak)
    
    Sk <- Slist[[k]]
    Bk <- index$class[[k]]
    
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
    
    M[Bk, Bk] <-
      M[Bk, Bk] +
      crossprod(
        Sk * sqrt(wk)
      )
  }
  
  symmetrize(M / n)
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
  
  IY <- crossprod(score_obj$bar) / n
  IY <- symmetrize(IY)
  
  Jw <- weighted_J_information(
    score_obj,
    w
  )
  
  symmetrize(
    IY + Jw
  )
}

############################################################
# 8. Empirical active-face curvature H_R
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
      
      kern <- dnorm(ell[idx] / h) / h
      
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
      
      D[, index$alpha] <-
        matrix(
          A[k, ] - A[l, ],
          nrow = length(idx),
          ncol = index$n_alpha,
          byrow = TRUE
        )
      
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
# 9. Efficient pointwise marginal value
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
  
  term1 <- numeric(nc)
  term2 <- numeric(nc)
  
  Sc <- lapply(
    Slist,
    function(S){
      S[candidate, , drop = FALSE]
    }
  )
  
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
    
    term1 <- term1 + tau_c[, k] * qkk
    
    term2 <- term2 +
      (tau_c[, k]^2) * qkk
  }
  
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
      
      v_l <- as.numeric(
        G[Bl, index$alpha, drop = FALSE] %*% ak
      )
      
      part_l <- as.numeric(
        Sl %*% v_l
      )
      
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
  
  psi[psi < 0 & psi > -1e-8] <- 0
  
  psi
}

############################################################
# 10. Simple acquisition scores
############################################################

score_entropy <- function(tau){
  -rowSums(
    tau * log(pmax(tau, 1e-15))
  )
}

score_margin <- function(tau){
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

############################################################
# 11. Design objective and Frank--Wolfe solver
############################################################

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
  
  sum(L * t(Iinv))
}

solve_finite_pool_design <- function(
    score_obj,
    pilot_labelled,
    candidate,
    L,
    B,
    maxit = FW_MAXIT,
    tol = FW_TOL,
    verbose = FALSE
){
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
  
  n_total <- nrow(score_obj$tau)
  
  a <- rep(B / N, N)
  
  Icur <- build_information_from_design(
    score_obj,
    pilot_labelled,
    candidate,
    a
  )
  
  Iinv <- safe_inverse(Icur)
  
  obj <- sum(L * t(Iinv))
  
  converged <- FALSE
  fw_gap <- Inf
  relative_fw_gap <- Inf
  last_step <- NA_real_
  
  for(iter in seq_len(maxit)){
    
    Iinv <- safe_inverse(Icur)
    
    G <-
      Iinv %*%
      L %*%
      Iinv
    
    G <- symmetrize(G)
    
    psi <- classification_value_from_G(
      score_obj,
      G,
      candidate
    )
    
    if(any(!is.finite(psi))){
      stop("Non-finite Frank-Wolfe gradient.")
    }
    
    ord <- order(
      psi,
      decreasing = TRUE
    )
    
    s <- numeric(N)
    s[ord[seq_len(B)]] <- 1
    
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
        "\n",
        sep = ""
      )
    }
    
    if(
      is.finite(relative_fw_gap) &&
      fw_gap >= 0 &&
      relative_fw_gap <= tol
    ){
      converged <- TRUE
      last_step <- 0
      break
    }
    
    Is <- build_information_from_design(
      score_obj,
      pilot_labelled,
      candidate,
      s
    )
    
    line_objective <- function(gamma){
      Ig <-
        (1 - gamma) * Icur +
        gamma * Is
      
      Ig_inv <- safe_inverse(Ig)
      
      sum(L * t(Ig_inv))
    }
    
    ls <- optimize(
      line_objective,
      interval = c(0, 1),
      tol = 1e-8
    )
    
    gamma <- ls$minimum
    new_obj <- ls$objective
    
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
      break
    }
    
    a <-
      (1 - gamma) * a +
      gamma * s
    
    a <- pmin(
      1,
      pmax(
        0,
        a
      )
    )
    
    Icur <-
      (1 - gamma) * Icur +
      gamma * Is
    
    Icur <- symmetrize(Icur)
    
    obj <- new_obj
  }
  
  Icur <- build_information_from_design(
    score_obj,
    pilot_labelled,
    candidate,
    a
  )
  
  Iinv <- safe_inverse(Icur)
  obj <- sum(L * t(Iinv))
  
  G <-
    Iinv %*%
    L %*%
    Iinv
  
  G <- symmetrize(G)
  
  psi <- classification_value_from_G(
    score_obj,
    G,
    candidate
  )
  
  ord <- order(
    psi,
    decreasing = TRUE
  )
  
  s <- numeric(N)
  s[ord[seq_len(B)]] <- 1
  
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
  
  a_round <- numeric(length(candidate))
  
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
# 12. Evaluation
############################################################

evaluate_fit <- function(
    fit,
    X_test,
    y_test,
    classes
){
  pred <- predict_ss_qda(
    fit,
    X_test
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
# 13. Load the official benchmark split
############################################################

if(!file.exists(TRAIN_FILE)){
  stop("Cannot find training file: ", TRAIN_FILE)
}

if(!file.exists(TEST_FILE)){
  stop("Cannot find test file: ", TEST_FILE)
}

train_raw <- as.matrix(
  read.table(
    TRAIN_FILE,
    header = FALSE
  )
)

test_raw <- as.matrix(
  read.table(
    TEST_FILE,
    header = FALSE
  )
)

if(ncol(train_raw) != 37L || ncol(test_raw) != 37L){
  stop(
    "Expected 37 columns in each file: 36 predictors plus class label."
  )
}

if(nrow(train_raw) != 4435L){
  warning("Training set has ", nrow(train_raw), " rows; expected 4435.")
}

if(nrow(test_raw) != 2000L){
  warning("Test set has ", nrow(test_raw), " rows; expected 2000.")
}

X_train_raw <- train_raw[, FEATURE_COLS, drop = FALSE]
X_test_raw  <- test_raw[, FEATURE_COLS, drop = FALSE]

y_train_raw <- as.integer(train_raw[, 37L])
y_test_raw  <- as.integer(test_raw[, 37L])

classes <- c("1", "2", "3", "4", "5", "7")

if(!setequal(as.character(sort(unique(y_train_raw))), classes)){
  stop(
    "Unexpected training classes: ",
    paste(sort(unique(y_train_raw)), collapse = ", ")
  )
}

if(!setequal(as.character(sort(unique(y_test_raw))), classes)){
  stop(
    "Unexpected test classes: ",
    paste(sort(unique(y_test_raw)), collapse = ", ")
  )
}

y_train <- factor(
  as.character(y_train_raw),
  levels = classes
)

y_test <- factor(
  as.character(y_test_raw),
  levels = classes
)

############################################################
# 14. Training-only standardization
#
# All training features are observed before label acquisition,
# so using the full benchmark training feature pool to estimate
# the unsupervised scaling transformation is legitimate.
############################################################

train_center <- colMeans(X_train_raw)
train_scale <- apply(X_train_raw, 2L, sd)

if(any(!is.finite(train_scale)) || any(train_scale <= 0)){
  stop("A selected spectral feature has zero or non-finite SD.")
}

X_train <- sweep(
  X_train_raw,
  2L,
  train_center,
  "-"
)

X_train <- sweep(
  X_train,
  2L,
  train_scale,
  "/"
)

X_test <- sweep(
  X_test_raw,
  2L,
  train_center,
  "-"
)

X_test <- sweep(
  X_test,
  2L,
  train_scale,
  "/"
)

colnames(X_train) <- paste0("Band", 1:4)
colnames(X_test) <- paste0("Band", 1:4)

cat("\n============================================================\n")
cat("LANDSAT BENCHMARK DATA\n")
cat("============================================================\n")
cat("Training rows:", nrow(X_train), "\n")
cat("Test rows:", nrow(X_test), "\n")
cat("Selected attributes:", paste(FEATURE_COLS, collapse = ", "), "\n")
cat("Training class counts:\n")
print(table(y_train))
cat("Test class counts:\n")
print(table(y_test))

############################################################
# 15. Repeated pilot/acquisition experiment
############################################################

n_train <- nrow(X_train)
d <- ncol(X_train)
g <- length(classes)

pilot_n <- ceiling(PILOT_FRAC * n_train)
budget_n <- ceiling(BUDGETS * n_train)

cat("\nPilot size:", pilot_n, "\n")
cat(
  "Total label counts:",
  paste(
    paste0(
      sprintf("%.0f%%", 100 * BUDGETS),
      "=",
      budget_n
    ),
    collapse = ", "
  ),
  "\n"
)

raw_results <- list()
class_results <- list()
boundary_results <- list()
design_results <- list()
failure_log <- list()

raw_rr <- 1L
class_rr <- 1L
bound_rr <- 1L
design_rr <- 1L
fail_rr <- 1L

for(rep_id in seq_len(N_REP)){
  
  if(VERBOSE){
    cat(
      "\n============================================================\n",
      "Replication ", rep_id, " / ", N_REP, "\n",
      "============================================================\n",
      sep = ""
    )
  }
  
  set.seed(
    SEED_MASTER + 10000L * rep_id
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
  
  pilot_counts <- table(
    factor(
      y_train[pilot_idx],
      levels = classes
    )
  )
  
  if(VERBOSE){
    cat("Pilot class counts:\n")
    print(pilot_counts)
  }
  
  if(any(pilot_counts <= d + 1L)){
    
    failure_log[[fail_rr]] <- data.frame(
      Replication = rep_id,
      Method = "ALL",
      Budget = NA_real_,
      Reason = paste0(
        "Insufficient pilot count: ",
        paste(
          names(pilot_counts),
          as.integer(pilot_counts),
          collapse = ", "
        )
      )
    )
    
    fail_rr <- fail_rr + 1L
    
    warning(
      "Skipping replication ",
      rep_id,
      " due to insufficient pilot class count."
    )
    
    next
  }
  
  ##########################################################
  # Shared pilot fit
  ##########################################################
  
  pilot_fit <- tryCatch(
    fit_ss_qda(
      X_train,
      y_train,
      pilot_labelled,
      classes
    ),
    error = function(e) e
  )
  
  if(inherits(pilot_fit, "error")){
    
    failure_log[[fail_rr]] <- data.frame(
      Replication = rep_id,
      Method = "ALL",
      Budget = NA_real_,
      Reason = paste0(
        "Pilot fit failed: ",
        conditionMessage(pilot_fit)
      )
    )
    
    fail_rr <- fail_rr + 1L
    next
  }
  
  pilot_pred <- predict_ss_qda(
    pilot_fit,
    X_train
  )
  
  entropy_all <- score_entropy(
    pilot_pred$posterior
  )
  
  margin_all <- score_margin(
    pilot_pred$posterior
  )
  
  index <- build_param_index(
    g = g,
    d = d
  )
  
  score_obj <- compute_score_components(
    X_train,
    pilot_fit,
    index
  )
  
  ICC_fit <- weighted_E_tt(
    score_obj,
    rep(1, n_train)
  )
  
  ##########################################################
  # Pilot-fitted classification-risk curvature
  ##########################################################
  
  Hout <- tryCatch(
    estimate_HR(
      score_obj
    ),
    error = function(e) e
  )
  
  if(inherits(Hout, "error")){
    
    failure_log[[fail_rr]] <- data.frame(
      Replication = rep_id,
      Method = "AdaptiveRisk",
      Budget = NA_real_,
      Reason = paste0(
        "H_R estimation failed: ",
        conditionMessage(Hout)
      )
    )
    
    fail_rr <- fail_rr + 1L
    next
  }
  
  Hfit <- Hout$H
  
  if(
    !all(is.finite(Hfit)) ||
    sum(abs(Hfit)) <= 0
  ){
    
    failure_log[[fail_rr]] <- data.frame(
      Replication = rep_id,
      Method = "AdaptiveRisk",
      Budget = NA_real_,
      Reason = "H_R is non-finite or numerically zero."
    )
    
    fail_rr <- fail_rr + 1L
    next
  }
  
  if(!is.null(Hout$diagnostics)){
    
    dd <- Hout$diagnostics
    
    dd$Replication <- rep_id
    dd$ClassK <- classes[dd$k]
    dd$ClassL <- classes[dd$l]
    
    boundary_results[[bound_rr]] <- dd
    bound_rr <- bound_rr + 1L
  }
  
  ##########################################################
  # Each budget independently from the same pilot
  ##########################################################
  
  for(bi in seq_along(BUDGETS)){
    
    budget <- BUDGETS[bi]
    total_n <- budget_n[bi]
    
    additional_B <- total_n - pilot_n
    
    if(additional_B <= 0L){
      stop(
        "PILOT_FRAC must be smaller than every total budget."
      )
    }
    
    additional_B <- min(
      additional_B,
      length(candidate)
    )
    
    if(VERBOSE){
      cat(
        "\nBudget ",
        sprintf("%.0f%%", 100 * budget),
        " (additional labels=",
        additional_B,
        ")\n",
        sep = ""
      )
    }
    
    ########################################################
    # Random
    ########################################################
    
    set.seed(
      SEED_MASTER +
        20000L * rep_id +
        bi
    )
    
    sel_random <- sample(
      candidate,
      size = additional_B,
      replace = FALSE
    )
    
    ########################################################
    # Entropy
    ########################################################
    
    sel_entropy <- candidate[
      order(
        entropy_all[candidate],
        decreasing = TRUE
      )[seq_len(additional_B)]
    ]
    
    ########################################################
    # Margin
    ########################################################
    
    sel_margin <- candidate[
      order(
        margin_all[candidate],
        decreasing = TRUE
      )[seq_len(additional_B)]
    ]
    
    ########################################################
    # Fisher
    ########################################################
    
    Fisher_design <- tryCatch(
      solve_finite_pool_design(
        score_obj = score_obj,
        pilot_labelled = pilot_labelled,
        candidate = candidate,
        L = ICC_fit,
        B = additional_B,
        maxit = FW_MAXIT,
        tol = FW_TOL,
        verbose = FALSE
      ),
      error = function(e) e
    )
    
    if(inherits(Fisher_design, "error")){
      failure_log[[fail_rr]] <- data.frame(
        Replication = rep_id,
        Method = "Fisher",
        Budget = budget,
        Reason = conditionMessage(Fisher_design)
      )
      fail_rr <- fail_rr + 1L
      next
    }
    
    Fisher_round <- round_design_topB(
      design = Fisher_design,
      score_obj = score_obj,
      pilot_labelled = pilot_labelled,
      candidate = candidate,
      L = ICC_fit,
      B = additional_B
    )
    
    sel_Fisher <- candidate[
      Fisher_round$selected
    ]
    
    ########################################################
    # Adaptive classification-risk design
    ########################################################
    
    Adaptive_design <- tryCatch(
      solve_finite_pool_design(
        score_obj = score_obj,
        pilot_labelled = pilot_labelled,
        candidate = candidate,
        L = Hfit,
        B = additional_B,
        maxit = FW_MAXIT,
        tol = FW_TOL,
        verbose = FALSE
      ),
      error = function(e) e
    )
    
    if(inherits(Adaptive_design, "error")){
      failure_log[[fail_rr]] <- data.frame(
        Replication = rep_id,
        Method = "AdaptiveRisk",
        Budget = budget,
        Reason = conditionMessage(Adaptive_design)
      )
      fail_rr <- fail_rr + 1L
      next
    }
    
    Adaptive_round <- round_design_topB(
      design = Adaptive_design,
      score_obj = score_obj,
      pilot_labelled = pilot_labelled,
      candidate = candidate,
      L = Hfit,
      B = additional_B
    )
    
    sel_Adaptive <- candidate[
      Adaptive_round$selected
    ]
    
    ########################################################
    # Design diagnostics
    ########################################################
    
    design_results[[design_rr]] <- data.frame(
      Replication = rep_id,
      Budget = budget,
      Method = c(
        "Fisher",
        "AdaptiveRisk"
      ),
      AdditionalBudget = additional_B,
      Iterations = c(
        Fisher_design$iterations,
        Adaptive_design$iterations
      ),
      Converged = c(
        Fisher_design$converged,
        Adaptive_design$converged
      ),
      FWGap = c(
        Fisher_design$fw_gap,
        Adaptive_design$fw_gap
      ),
      RelativeFWGap = c(
        Fisher_design$relative_fw_gap,
        Adaptive_design$relative_fw_gap
      ),
      LastFWStep = c(
        Fisher_design$last_step,
        Adaptive_design$last_step
      ),
      FractionalCount = c(
        Fisher_round$fractional_count,
        Adaptive_round$fractional_count
      ),
      FractionalProp = c(
        Fisher_round$fractional_prop,
        Adaptive_round$fractional_prop
      ),
      RelaxedObjective = c(
        Fisher_round$relaxed_objective,
        Adaptive_round$relaxed_objective
      ),
      RoundedObjective = c(
        Fisher_round$rounded_objective,
        Adaptive_round$rounded_objective
      ),
      RoundingLossPct = c(
        Fisher_round$rounding_loss_pct,
        Adaptive_round$rounding_loss_pct
      )
    )
    
    design_rr <- design_rr + 1L
    
    ########################################################
    # Refit and evaluate all methods
    ########################################################
    
    selections <- list(
      Random = sel_random,
      Entropy = sel_entropy,
      Margin = sel_margin,
      Fisher = sel_Fisher,
      AdaptiveRisk = sel_Adaptive
    )
    
    for(method in METHODS){
      
      labelled <- pilot_labelled
      labelled[selections[[method]]] <- TRUE
      
      fit_budget <- tryCatch(
        fit_ss_qda(
          X_train,
          y_train,
          labelled,
          classes
        ),
        error = function(e) e
      )
      
      if(inherits(fit_budget, "error")){
        
        failure_log[[fail_rr]] <- data.frame(
          Replication = rep_id,
          Method = method,
          Budget = budget,
          Reason = paste0(
            "Final fit failed: ",
            conditionMessage(fit_budget)
          )
        )
        
        fail_rr <- fail_rr + 1L
        next
      }
      
      ev <- evaluate_fit(
        fit_budget,
        X_test,
        y_test,
        classes
      )
      
      raw_results[[raw_rr]] <- data.frame(
        Replication = rep_id,
        Budget = budget,
        Method = method,
        NLabelled = sum(labelled),
        Error = ev$error,
        BalancedError = ev$balanced_error,
        LogLik = fit_budget$logLik,
        EMIterations = fit_budget$iterations
      )
      
      raw_rr <- raw_rr + 1L
      
      class_results[[class_rr]] <- data.frame(
        Replication = rep_id,
        Budget = budget,
        Method = method,
        Class = classes,
        ClassError = as.numeric(ev$class_error)
      )
      
      class_rr <- class_rr + 1L
    }
  }
  
  ##########################################################
  # Save progress after every completed replication
  ##########################################################
  
  if(length(raw_results) > 0){
    write.csv(
      bind_rows(raw_results),
      "Landsat_raw_PROGRESS.csv",
      row.names = FALSE
    )
  }
  
  if(length(design_results) > 0){
    write.csv(
      bind_rows(design_results),
      "Landsat_design_diagnostics_PROGRESS.csv",
      row.names = FALSE
    )
  }
}

############################################################
# 16. Final outputs
############################################################

raw <- bind_rows(raw_results)
class_raw <- bind_rows(class_results)
boundary_raw <- bind_rows(boundary_results)
design_raw <- bind_rows(design_results)
failures <- bind_rows(failure_log)

write.csv(
  raw,
  "Landsat_raw.csv",
  row.names = FALSE
)

write.csv(
  class_raw,
  "Landsat_class_errors.csv",
  row.names = FALSE
)

write.csv(
  boundary_raw,
  "Landsat_boundary_diagnostics.csv",
  row.names = FALSE
)

write.csv(
  design_raw,
  "Landsat_design_diagnostics.csv",
  row.names = FALSE
)

write.csv(
  failures,
  "Landsat_failures.csv",
  row.names = FALSE
)

############################################################
# 17. Mean errors and Monte Carlo standard errors
############################################################

summary_results <-
  raw %>%
  group_by(
    Budget,
    Method
  ) %>%
  summarise(
    R = n(),
    MeanError = mean(Error),
    SDError = sd(Error),
    MCSEError = SDError / sqrt(R),
    MeanBalancedError = mean(BalancedError),
    SDBalancedError = sd(BalancedError),
    MCSEBalancedError = SDBalancedError / sqrt(R),
    .groups = "drop"
  )

write.csv(
  summary_results,
  "Landsat_summary.csv",
  row.names = FALSE
)

############################################################
# 18. Paired differences against AdaptiveRisk
#
# Difference = competitor error - AdaptiveRisk error.
# Positive values favor AdaptiveRisk.
############################################################

paired_list <- list()
pp <- 1L

for(budget in BUDGETS){
  
  dat_b <- raw[
    raw$Budget == budget,
    ,
    drop = FALSE
  ]
  
  adapt <- dat_b[
    dat_b$Method == "AdaptiveRisk",
    c("Replication", "Error", "BalancedError"),
    drop = FALSE
  ]
  
  names(adapt)[2:3] <- c(
    "AdaptiveError",
    "AdaptiveBalancedError"
  )
  
  for(comp in setdiff(METHODS, "AdaptiveRisk")){
    
    dc <- dat_b[
      dat_b$Method == comp,
      c("Replication", "Error", "BalancedError"),
      drop = FALSE
    ]
    
    names(dc)[2:3] <- c(
      "CompError",
      "CompBalancedError"
    )
    
    mm <- merge(
      adapt,
      dc,
      by = "Replication"
    )
    
    d_err <- mm$CompError - mm$AdaptiveError
    d_bal <- mm$CompBalancedError - mm$AdaptiveBalancedError
    
    n_pair <- length(d_err)
    
    mean_d <- mean(d_err)
    se_d <- sd(d_err) / sqrt(n_pair)
    
    mean_db <- mean(d_bal)
    se_db <- sd(d_bal) / sqrt(n_pair)
    
    paired_list[[pp]] <- data.frame(
      Budget = budget,
      Comparator = comp,
      NPairs = n_pair,
      MeanDifference = mean_d,
      MCSEDifference = se_d,
      Lower95MC = mean_d - 1.96 * se_d,
      Upper95MC = mean_d + 1.96 * se_d,
      MeanBalancedDifference = mean_db,
      MCSEBalancedDifference = se_db,
      Lower95MCBalanced = mean_db - 1.96 * se_db,
      Upper95MCBalanced = mean_db + 1.96 * se_db
    )
    
    pp <- pp + 1L
  }
}

paired_results <- bind_rows(paired_list)

write.csv(
  paired_results,
  "Landsat_paired_vs_AdaptiveRisk.csv",
  row.names = FALSE
)

############################################################
# 19. Design diagnostic summary
############################################################

design_summary <-
  design_raw %>%
  group_by(
    Budget,
    Method
  ) %>%
  summarise(
    R = n(),
    ConvergenceRate = mean(Converged),
    MeanIterations = mean(Iterations),
    MaxRelativeFWGap = max(RelativeFWGap),
    MeanFractionalProp = mean(FractionalProp),
    MaxFractionalProp = max(FractionalProp),
    MeanRoundingLossPct = mean(RoundingLossPct),
    MaxRoundingLossPct = max(RoundingLossPct),
    .groups = "drop"
  )

write.csv(
  design_summary,
  "Landsat_design_diagnostics_summary.csv",
  row.names = FALSE
)

############################################################
# 20. Performance figure
############################################################

plot_dat <-
  summary_results %>%
  mutate(
    Method = factor(
      Method,
      levels = METHODS
    )
  )

fig <-
  ggplot(
    plot_dat,
    aes(
      x = Budget,
      y = MeanError,
      group = Method,
      linetype = Method,
      shape = Method
    )
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.4
  ) +
  geom_errorbar(
    aes(
      ymin = MeanError - 1.96 * MCSEError,
      ymax = MeanError + 1.96 * MCSEError
    ),
    width = 0.008,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = BUDGETS,
    labels = paste0(100 * BUDGETS, "%")
  ) +
  labs(
    x = "Total labeling budget",
    y = "Benchmark test classification error",
    linetype = "Acquisition rule",
    shape = "Acquisition rule"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  "Landsat_performance.png",
  fig,
  width = 9.5,
  height = 6.5,
  dpi = 400
)

############################################################
# 21. Console output
############################################################

cat("\n\n============================================================\n")
cat("FINAL LANDSAT SUMMARY\n")
cat("============================================================\n")
print(summary_results)

cat("\n\n============================================================\n")
cat("PAIRED DIFFERENCES: comparator - AdaptiveRisk\n")
cat("Positive values favor AdaptiveRisk\n")
cat("============================================================\n")
print(paired_results)

cat("\n\n============================================================\n")
cat("DESIGN DIAGNOSTICS\n")
cat("============================================================\n")
print(design_summary)

cat("\n\n============================================================\n")
cat("FAILURES\n")
cat("============================================================\n")

if(nrow(failures) == 0L){
  cat("No failures recorded.\n")
} else {
  print(failures)
}

cat("\n\nRun complete.\n")
