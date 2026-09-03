############################################################
# LANDSAT REAL-DATA GEOMETRY / ACQUISITION FIGURE
# Classification-Risk-Optimal Label Acquisition
#
# This script reproduces ONE objectively chosen replication from
# the frozen 50-replication Landsat study.
#
# IMPORTANT DESIGN CHOICES:
#
#   - The official UCI sat.trn / sat.tst split is kept fixed.
#   - The fitted classification model uses ONLY attributes 17--20
#     (the four central-pixel spectral bands), exactly as in the
#     final Landsat experiment.
#   - The model itself is 4-dimensional Gaussian QDA.
#   - PC1 and PC2 below are used ONLY for visualization.
#   - Entropy, H_R, classification value, and acquisition are all
#     computed from the full 4-dimensional QDA model.
#
# Representative replication:
#
#   We use the 20% total labeling budget and select the replication
#   whose paired AdaptiveRisk-versus-Entropy improvement is closest
#   to the median paired improvement over all 50 replications:
#
#       Entropy error - AdaptiveRisk error.
#
#   This is objective and does NOT choose a replication because the
#   proposed method happens to look unusually favorable.
#
# Why 20%?
#   In the final Landsat experiment, 20% is the most stable moderate
#   label-scarcity regime. AdaptiveRisk and entropy have very similar
#   mean error, while AdaptiveRisk clearly improves on the Fisher
#   information design. Hence this budget is appropriate for showing
#   the geometry rather than a cherry-picked performance extreme.
#
# Outputs:
#   Landsat_geometry_acquisition.pdf
#   Landsat_geometry_acquisition.png
#   Landsat_geometry_plot_data.csv
#   Landsat_geometry_diagnostics.csv
#   Landsat_figure_replication_selection.csv
############################################################

rm(list = ls())
gc()

############################################################
# 0. Packages
############################################################

required_packages <- c(
  "ggplot2",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if(length(missing_packages) > 0){
  install.packages(missing_packages)
}

library(ggplot2)
library(patchwork)
library(scales)

############################################################
# 1. User settings
############################################################
TRAIN_FILE <- file.path("data", "landsat", "sat.trn")
TEST_FILE  <- file.path("data", "landsat", "sat.tst")
FINAL_RAW_FILE <- file.path("results", "Landsat_raw.csv")

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR)

SEED_MASTER <- 20260902L

FEATURE_COLS <- 17:20

PILOT_FRAC <- 0.05

FIG_BUDGET <- 0.20

CLASSES <- c("1", "2", "3", "4", "5", "7")

EM_MAXIT <- 400L
EM_TOL <- 1e-7

COV_FLOOR_REL <- 1e-6
INFO_EIG_FLOOR_REL <- 1e-8

BOUNDARY_BW_MIN <- 0.05
BOUNDARY_BW_MAX <- 1.00
BOUNDARY_KERNEL_CUTOFF <- 4.0
MIN_BOUNDARY_POINTS <- 10L

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
# 4. Semi-supervised QDA
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
  
  global_S <- regularize_cov(
    cov(X),
    cov_floor_rel
  )
  
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
    
    mu_k[k, ] <- colMeans(
      X[idx, , drop = FALSE]
    )
    
    Sk <- cov(
      X[idx, , drop = FALSE]
    )
    
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
    
    log_joint <- compute_log_joint(
      X,
      fit_now
    )
    
    resp[,] <- softmax_rows(
      log_joint
    )
    
    resp[lab_idx, ] <- 0
    
    lab_class_num <- match(
      as.character(y[lab_idx]),
      classes
    )
    
    resp[
      cbind(
        lab_idx,
        lab_class_num
      )
    ] <- 1
    
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
    
    Sigma_new <- vector(
      "list",
      g
    )
    
    for(k in seq_len(g)){
      wk <- resp[, k]
      
      mu_new[k, ] <-
        colSums(X * wk) /
        Nk[k]
      
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
      sum(
        log_sum_exp_rows(
          log_joint[
            unl_idx,
            ,
            drop = FALSE
          ]
        )
      )
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
      rel_change <-
        abs(ll - old_ll) /
        (1 + abs(old_ll))
      
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
    pilot_counts = pilot_counts
  )
}

############################################################
# 5. Parameter indexing
############################################################

vech_pairs <- function(d){
  which(
    lower.tri(
      matrix(0, d, d),
      diag = TRUE
    ),
    arr.ind = TRUE
  )
}

build_param_index <- function(g, d){
  n_alpha <- g - 1L
  n_cov <- d * (d + 1L) / 2L
  
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
    class_blocks[[k]] <- c(
      mu_idx,
      cov_idx
    )
  }
  
  list(
    alpha = alpha_idx,
    mu = mu_blocks,
    cov = cov_blocks,
    class = class_blocks,
    n_alpha = n_alpha,
    n_cov = n_cov,
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

cov_score_matrix <- function(
    X,
    mu,
    Sigma,
    pairs
){
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

compute_score_components <- function(
    X,
    fit,
    index
){
  n <- nrow(X)
  g <- length(fit$pi)
  
  pred <- predict_ss_qda(
    fit,
    X
  )
  
  tau <- pred$posterior
  log_joint <- pred$log_joint
  
  A <- prior_score_vectors(
    fit$pi
  )
  
  class_scores <- vector(
    "list",
    g
  )
  
  for(k in seq_len(g)){
    invS <- solve(
      fit$Sigma[[k]]
    )
    
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
# 7. Information calculations
############################################################

weighted_E_tt <- function(
    score_obj,
    weights
){
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

############################################################
# 8. Active-face curvature H_R
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
  
  tau2[
    cbind(
      seq_len(n),
      top1
    )
  ] <- -Inf
  
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
      
      ell <-
        log_joint[, k] -
        log_joint[, l]
      
      ell_active <- ell[active]
      
      s_ell <- sd(ell_active)
      
      if(
        !is.finite(s_ell) ||
        s_ell <= 0
      ){
        s_ell <- mad(
          ell_active,
          constant = 1
        )
      }
      
      h <-
        1.06 *
        s_ell *
        n_active^(-1 / 5)
      
      h <- max(
        BOUNDARY_BW_MIN,
        min(
          BOUNDARY_BW_MAX,
          h
        )
      )
      
      keep <-
        active &
        abs(ell) <=
        BOUNDARY_KERNEL_CUTOFF * h
      
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
      
      w <-
        sqrt(
          pmax(tau[idx, k], 0) *
            pmax(tau[idx, l], 0)
        ) *
        kern
      
      good <-
        is.finite(w) &
        w > 0
      
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
        Slist[[k]][
          idx,
          ,
          drop = FALSE
        ]
      
      D[, index$class[[l]]] <-
        -Slist[[l]][
          idx,
          ,
          drop = FALSE
        ]
      
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
  
  list(
    H = symmetrize(H),
    diagnostics = do.call(
      rbind,
      pair_diag
    )
  )
}

############################################################
# 9. Pointwise classification value
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
    candidate <- seq_len(
      nrow(tau)
    )
  }
  
  tau_c <- tau[
    candidate,
    ,
    drop = FALSE
  ]
  
  g <- ncol(tau)
  nc <- length(candidate)
  
  GAA <- G[
    index$alpha,
    index$alpha,
    drop = FALSE
  ]
  
  term1 <- numeric(nc)
  term2 <- numeric(nc)
  
  Sc <- lapply(
    Slist,
    function(S){
      S[
        candidate,
        ,
        drop = FALSE
      ]
    }
  )
  
  for(k in seq_len(g)){
    ak <- A[k, ]
    Bk <- index$class[[k]]
    Sk <- Sc[[k]]
    
    const <- as.numeric(
      t(ak) %*%
        GAA %*%
        ak
    )
    
    v_left <- as.numeric(
      G[
        Bk,
        index$alpha,
        drop = FALSE
      ] %*% ak
    )
    
    lin <- 2 * as.numeric(
      Sk %*% v_left
    )
    
    bil <- row_quad_bilinear(
      Sk,
      G[
        Bk,
        Bk,
        drop = FALSE
      ],
      Sk
    )
    
    qkk <- const + lin + bil
    
    term1 <-
      term1 +
      tau_c[, k] * qkk
    
    term2 <-
      term2 +
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
        t(ak) %*%
          GAA %*%
          al
      )
      
      v_l <- as.numeric(
        G[
          Bl,
          index$alpha,
          drop = FALSE
        ] %*% ak
      )
      
      part_l <- as.numeric(
        Sl %*% v_l
      )
      
      v_k <- as.numeric(
        G[
          Bk,
          index$alpha,
          drop = FALSE
        ] %*% al
      )
      
      part_k <- as.numeric(
        Sk %*% v_k
      )
      
      bil <- row_quad_bilinear(
        Sk,
        G[
          Bk,
          Bl,
          drop = FALSE
        ],
        Sl
      )
      
      qkl <-
        const +
        part_l +
        part_k +
        bil
      
      term2 <-
        term2 +
        2 *
        tau_c[, k] *
        tau_c[, l] *
        qkl
    }
  }
  
  psi <- term1 - term2
  
  psi[
    psi < 0 &
      psi > -1e-8
  ] <- 0
  
  psi
}

score_entropy <- function(tau){
  -rowSums(
    tau *
      log(
        pmax(
          tau,
          1e-15
        )
      )
  )
}

############################################################
# 10. Finite-pool convex design solver
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
  
  Iinv <- safe_inverse(
    Icur
  )
  
  sum(
    L *
      t(Iinv)
  )
}

solve_finite_pool_design <- function(
    score_obj,
    pilot_labelled,
    candidate,
    L,
    B,
    maxit = 1000L,
    tol = 1e-5,
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
    
    return(
      list(
        a = a0,
        objective = design_objective(
          score_obj,
          pilot_labelled,
          candidate,
          a0,
          L
        ),
        iterations = 0L,
        converged = TRUE,
        fw_gap = 0,
        relative_fw_gap = 0
      )
    )
  }
  
  n_total <- nrow(
    score_obj$tau
  )
  
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
    L *
      t(Iinv)
  )
  
  converged <- FALSE
  fw_gap <- Inf
  relative_fw_gap <- Inf
  
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
    
    ord <- order(
      psi,
      decreasing = TRUE
    )
    
    s <- numeric(N)
    
    s[
      ord[seq_len(B)]
    ] <- 1
    
    fw_gap <-
      sum(
        (s - a) *
          psi
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
    
    if(
      is.finite(relative_fw_gap) &&
      fw_gap >= 0 &&
      relative_fw_gap <= tol
    ){
      converged <- TRUE
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
      
      Ig_inv <- safe_inverse(
        Ig
      )
      
      sum(
        L *
          t(Ig_inv)
      )
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
    
    if(
      obj0 <= new_obj &&
      obj0 <= obj1
    ){
      gamma <- 0
      new_obj <- obj0
    } else if(
      obj1 < new_obj &&
      obj1 < obj0
    ){
      gamma <- 1
      new_obj <- obj1
    }
    
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
    
    Icur <- symmetrize(
      Icur
    )
    
    obj <- new_obj
  }
  
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
    L *
      t(Iinv)
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
  
  s <- numeric(N)
  
  s[
    ord[seq_len(B)]
  ] <- 1
  
  fw_gap <-
    sum(
      (s - a) *
        psi
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
    relative_fw_gap = relative_fw_gap
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
  selected <- order(
    design$a,
    decreasing = TRUE
  )[seq_len(B)]
  
  a_round <- numeric(
    length(candidate)
  )
  
  a_round[selected] <- 1
  
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
# 11. Load official Landsat data
############################################################

if(!file.exists(TRAIN_FILE)){
  stop("Cannot find ", TRAIN_FILE)
}

if(!file.exists(TEST_FILE)){
  stop("Cannot find ", TEST_FILE)
}

if(!file.exists(FINAL_RAW_FILE)){
  stop(
    paste0(
      "Cannot find ",
      FINAL_RAW_FILE,
      ". Run the final 50-replication Landsat experiment first."
    )
  )
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

X_train_raw <- train_raw[
  ,
  FEATURE_COLS,
  drop = FALSE
]

X_test_raw <- test_raw[
  ,
  FEATURE_COLS,
  drop = FALSE
]

y_train <- factor(
  as.character(
    as.integer(
      train_raw[, 37L]
    )
  ),
  levels = CLASSES
)

y_test <- factor(
  as.character(
    as.integer(
      test_raw[, 37L]
    )
  ),
  levels = CLASSES
)

classes <- CLASSES

############################################################
# 12. Training-only standardization
############################################################

train_center <- colMeans(
  X_train_raw
)

train_scale <- apply(
  X_train_raw,
  2L,
  sd
)

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

colnames(X_train) <-
  paste0(
    "Band",
    1:4
  )

colnames(X_test) <-
  paste0(
    "Band",
    1:4
  )

############################################################
# 13. Choose representative replication objectively
############################################################

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
    "Landsat_raw.csv does not contain the required columns."
  )
}

ad <- final_raw[
  abs(
    final_raw$Budget -
      FIG_BUDGET
  ) < 1e-12 &
    final_raw$Method ==
    "AdaptiveRisk",
  c(
    "Replication",
    "Error"
  ),
  drop = FALSE
]

names(ad)[2] <-
  "AdaptiveError"

en <- final_raw[
  abs(
    final_raw$Budget -
      FIG_BUDGET
  ) < 1e-12 &
    final_raw$Method ==
    "Entropy",
  c(
    "Replication",
    "Error"
  ),
  drop = FALSE
]

names(en)[2] <-
  "EntropyError"

paired_rep <- merge(
  ad,
  en,
  by = "Replication"
)

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

TARGET_REP <-
  paired_rep$Replication[1L]

write.csv(
  paired_rep,
  "Landsat_figure_replication_selection.csv",
  row.names = FALSE
)

cat("\n============================================================\n")
cat("LANDSAT REPRESENTATIVE VISUALIZATION\n")
cat("============================================================\n")
cat(
  "Target budget: ",
  sprintf(
    "%.0f%%",
    100 * FIG_BUDGET
  ),
  "\n",
  sep = ""
)
cat(
  "Median paired improvement (Entropy - AdaptiveRisk): ",
  signif(
    median_improvement,
    6
  ),
  "\n",
  sep = ""
)
cat(
  "Selected replication: ",
  TARGET_REP,
  "\n",
  sep = ""
)
cat(
  "Selected-rep improvement: ",
  signif(
    paired_rep$Improvement[1L],
    6
  ),
  "\n",
  sep = ""
)
cat(
  "AdaptiveRisk test error: ",
  signif(
    paired_rep$AdaptiveError[1L],
    6
  ),
  "\n",
  sep = ""
)
cat(
  "Entropy test error: ",
  signif(
    paired_rep$EntropyError[1L],
    6
  ),
  "\n",
  sep = ""
)

############################################################
# 14. Reproduce the selected replication exactly
############################################################

n_train <- nrow(
  X_train
)

pilot_n <- ceiling(
  PILOT_FRAC *
    n_train
)

total_n <- ceiling(
  FIG_BUDGET *
    n_train
)

additional_B <-
  total_n -
  pilot_n

rep_id <- TARGET_REP

set.seed(
  SEED_MASTER +
    10000L *
    rep_id
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

pilot_labelled[
  pilot_idx
] <- TRUE

candidate <- which(
  !pilot_labelled
)

cat(
  "Training n: ",
  n_train,
  "\n",
  sep = ""
)

cat(
  "Pilot n: ",
  pilot_n,
  "\n",
  sep = ""
)

cat(
  "Additional labels at ",
  sprintf(
    "%.0f%%",
    100 * FIG_BUDGET
  ),
  " budget: ",
  additional_B,
  "\n",
  sep = ""
)

############################################################
# 15. Pilot fit and full 4D acquisition quantities
############################################################

pilot_fit <- fit_ss_qda(
  X_train,
  y_train,
  pilot_labelled,
  classes,
  maxit = EM_MAXIT,
  tol = EM_TOL
)

pilot_pred <- predict_ss_qda(
  pilot_fit,
  X_train
)

entropy_all <- score_entropy(
  pilot_pred$posterior
)

index <- build_param_index(
  g = length(classes),
  d = ncol(X_train)
)

score_obj <- compute_score_components(
  X_train,
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
  stop(
    "Estimated H_R is non-finite or numerically zero."
  )
}

############################################################
# 16. Entropy and AdaptiveRisk selections
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
  maxit = 1000L,
  tol = 1e-5,
  verbose = FALSE
)

if(!isTRUE(Adaptive_design$converged)){
  stop(
    paste0(
      "Adaptive design did not converge. Relative FW gap = ",
      signif(
        Adaptive_design$relative_fw_gap,
        6
      )
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

cat(
  "Adaptive design converged: ",
  Adaptive_design$converged,
  "\n",
  sep = ""
)

cat(
  "Relative FW gap: ",
  signif(
    Adaptive_design$relative_fw_gap,
    6
  ),
  "\n",
  sep = ""
)

cat(
  "Rounding loss (%): ",
  signif(
    Adaptive_round$rounding_loss_pct,
    6
  ),
  "\n",
  sep = ""
)

############################################################
# 17. Classification value at fitted adaptive design
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

psi_all[
  psi_all < 0 &
    psi_all > -1e-8
] <- 0

psi_nonneg <- pmax(
  psi_all,
  0
)

# Robust colour scaling for visualization only.
# This prevents a few extreme values from flattening the colour scale.
psi_cap <- as.numeric(
  quantile(
    psi_nonneg,
    0.99,
    na.rm = TRUE
  )
)

if(
  !is.finite(psi_cap) ||
  psi_cap <= 0
){
  psi_cap <- max(
    psi_nonneg
  )
}

relative_psi <- pmin(
  psi_nonneg / psi_cap,
  1
)

############################################################
# 18. Two-dimensional visualization coordinates
#
# PCA is used ONLY here for plotting.
# It is fitted to the full standardized benchmark training features.
############################################################

pca_vis <- prcomp(
  X_train,
  center = FALSE,
  scale. = FALSE
)

X_vis <- pca_vis$x[
  ,
  1:2,
  drop = FALSE
]

############################################################
# 19. Plot data
############################################################

class_names <- c(
  "1" = "Red soil",
  "2" = "Cotton crop",
  "3" = "Grey soil",
  "4" = "Damp grey soil",
  "5" = "Vegetation stubble",
  "7" = "Very damp grey soil"
)

plot_df <- data.frame(
  Row = seq_len(n_train),
  PC1 = X_vis[, 1L],
  PC2 = X_vis[, 2L],
  ClassCode = y_train,
  Class = factor(
    unname(
      class_names[
        as.character(y_train)
      ]
    ),
    levels = unname(
      class_names
    )
  ),
  Entropy = entropy_all,
  RelativeValue = relative_psi,
  Pilot = pilot_labelled,
  Candidate = !pilot_labelled,
  EntropySelected =
    seq_len(n_train) %in%
    sel_entropy,
  AdaptiveSelected =
    seq_len(n_train) %in%
    sel_adaptive
)

write.csv(
  plot_df,
  "Landsat_geometry_plot_data.csv",
  row.names = FALSE
)

selection_summary <- data.frame(
  Replication = TARGET_REP,
  Budget = FIG_BUDGET,
  PilotN = pilot_n,
  AdditionalBudget = additional_B,
  AdaptiveConverged =
    Adaptive_design$converged,
  AdaptiveIterations =
    Adaptive_design$iterations,
  AdaptiveFWGap =
    Adaptive_design$fw_gap,
  AdaptiveRelativeFWGap =
    Adaptive_design$relative_fw_gap,
  AdaptiveRoundingLossPct =
    Adaptive_round$rounding_loss_pct,
  MedianPairedImprovement =
    median_improvement,
  SelectedRepPairedImprovement =
    paired_rep$Improvement[1L]
)

write.csv(
  selection_summary,
  "Landsat_geometry_diagnostics.csv",
  row.names = FALSE
)

############################################################
# 20. Publication-quality four-panel figure
############################################################

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
    alpha = 0.46,
    size = 0.70
  ) +
  labs(
    title = "(a) Landsat feature geometry",
    x = expression(PC[1]),
    y = expression(PC[2]),
    colour = "Land-cover class"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        size = 2.2
      )
    )
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 8.7
    ),
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
    name = expression(
      "Relative " * psi[a](y)
    )
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
    Selected =
      plot_df$AdaptiveSelected,
    Rule =
      "Classification-risk optimal"
  ),
  data.frame(
    PC1 = plot_df$PC1,
    PC2 = plot_df$PC2,
    Selected =
      plot_df$EntropySelected,
    Rule =
      "Entropy"
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
    alpha = 0.28,
    size = 0.52
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
    title = paste0(
      "(d) Additional labels selected under a ",
      sprintf(
        "%.0f%%",
        100 * FIG_BUDGET
      ),
      " total budget"
    ),
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
      "Statlog Landsat Satellite data; representative replication ",
      TARGET_REP,
      " selected by the median paired ",
      sprintf(
        "%.0f%%",
        100 * FIG_BUDGET
      ),
      " AdaptiveRisk-versus-Entropy difference"
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
  "Landsat_geometry_acquisition.pdf",
  final_figure,
  width = 13.5,
  height = 9.5,
  device = cairo_pdf
)

ggsave(
  "Landsat_geometry_acquisition.png",
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
      "Landsat_geometry_acquisition.pdf",
      "Landsat_geometry_acquisition.png",
      "Landsat_geometry_plot_data.csv",
      "Landsat_geometry_diagnostics.csv",
      "Landsat_figure_replication_selection.csv"
    ),
    collapse = "\n"
  ),
  "\n"
)

cat("\nDONE.\n")
