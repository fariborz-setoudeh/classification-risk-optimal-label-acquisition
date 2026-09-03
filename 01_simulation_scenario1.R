# ============================================================
# STUDY 1
# Classification-risk-optimal label acquisition
#
# Three-class QDA with curved Bayes boundaries
#
# Outputs:
#   Study1_raw.csv
#   Study1_summary.csv
#   Figure1_geometry.png
#   Figure2_performance.png
#
# Methods:
#   Random
#   Entropy
#   Margin
#   A-optimal
#   Risk-optimal
# ============================================================

rm(list = ls())

set.seed(20260901)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# ============================================================
# 1. TRUE THREE-CLASS QDA MODEL
# ============================================================

g <- 3
p <- 2

pi0 <- c(
  0.36,
  0.34,
  0.30
)

mu0 <- rbind(
  c(-1.60, 0.00),
  c( 1.50, 0.20),
  c( 0.00, 1.90)
)

Sigma0 <- array(NA_real_, dim = c(p, p, g))

Sigma0[, , 1] <- matrix(
  c(
    1.00,  0.45,
    0.45,  0.70
  ),
  2, 2, byrow = TRUE
)

Sigma0[, , 2] <- matrix(
  c(
    0.80, -0.35,
    -0.35,  1.10
  ),
  2, 2, byrow = TRUE
)

Sigma0[, , 3] <- matrix(
  c(
    1.25, 0.15,
    0.15, 0.55
  ),
  2, 2, byrow = TRUE
)

true_par <- list(
  pi = pi0,
  mu = mu0,
  Sigma = Sigma0
)


# ============================================================
# 2. BASIC GAUSSIAN FUNCTIONS
# ============================================================

log_dmvnorm <- function(x, mu, Sigma) {
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  mu <- as.numeric(mu)
  
  Sigma <- as.matrix(Sigma)
  storage.mode(Sigma) <- "double"
  
  q <- ncol(x)
  
  if (length(mu) != q) {
    stop(
      "Dimension mismatch in log_dmvnorm(): length(mu) = ",
      length(mu),
      ", but ncol(x) = ",
      q
    )
  }
  
  R <- chol(Sigma)
  
  z <- x - matrix(
    mu,
    nrow = nrow(x),
    ncol = q,
    byrow = TRUE
  )
  
  z <- backsolve(
    R,
    t(z),
    transpose = TRUE
  )
  
  quad <- colSums(z^2)
  
  logdet <- 2 * sum(
    log(diag(R))
  )
  
  -0.5 * (
    q * log(2 * pi) +
      logdet +
      quad
  )
}


rmvnorm_fast <- function(n, mu, Sigma) {
  
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  
  q <- length(mu)
  
  Z <- matrix(
    rnorm(n * q),
    nrow = n,
    ncol = q
  )
  
  X <- Z %*% chol(Sigma)
  
  X + matrix(
    mu,
    nrow = n,
    ncol = q,
    byrow = TRUE
  )
}


simulate_qda <- function(n, pars) {
  
  z <- sample(
    seq_along(pars$pi),
    n,
    replace = TRUE,
    prob = pars$pi
  )
  
  Y <- matrix(
    NA_real_,
    nrow = n,
    ncol = ncol(pars$mu)
  )
  
  for (k in seq_along(pars$pi)) {
    
    ind <- which(z == k)
    
    if (length(ind) > 0) {
      
      Y[ind, ] <- rmvnorm_fast(
        length(ind),
        pars$mu[k, ],
        pars$Sigma[, , k]
      )
    }
  }
  
  list(
    Y = Y,
    z = z
  )
}


log_weighted_density <- function(Y, pars) {
  
  Y <- as.matrix(Y)
  
  out <- matrix(
    NA_real_,
    nrow(Y),
    length(pars$pi)
  )
  
  for (k in seq_along(pars$pi)) {
    
    out[, k] <-
      log(pars$pi[k]) +
      log_dmvnorm(
        Y,
        pars$mu[k, ],
        pars$Sigma[, , k]
      )
  }
  
  out
}


posterior_prob <- function(Y, pars) {
  
  L <- log_weighted_density(Y, pars)
  
  M <- apply(L, 1, max)
  
  E <- exp(L - M)
  
  E / rowSums(E)
}


predict_qda <- function(Y, pars) {
  
  max.col(
    log_weighted_density(Y, pars),
    ties.method = "first"
  )
}


# ============================================================
# 3. PARAMETER SCORE t_k(y)
#
# theta =
#   alpha_1,...,alpha_{g-1},
#   mu_1,...,mu_g,
#   vech(Sigma_1),...,vech(Sigma_g)
#
# For p = 2:
# vech(Sigma) = (sigma11, sigma21, sigma22)
# ============================================================

qdim <- (g - 1) + g * p + g * 3


class_score <- function(y, k, pars) {
  
  out <- rep(0, qdim)
  
  # ----------------------------------------------------------
  # mixing proportions: baseline logit parameterization
  # ----------------------------------------------------------
  
  alpha_score <- -pars$pi[1:(g - 1)]
  
  if (k < g) {
    alpha_score[k] <- alpha_score[k] + 1
  }
  
  out[1:(g - 1)] <- alpha_score
  
  
  # ----------------------------------------------------------
  # mean block
  # ----------------------------------------------------------
  
  mean_start <-
    (g - 1) +
    (k - 1) * p +
    1
  
  mean_ind <- mean_start:(mean_start + p - 1)
  
  Sinv <- solve(pars$Sigma[, , k])
  
  d <- y - pars$mu[k, ]
  
  out[mean_ind] <- as.vector(
    Sinv %*% d
  )
  
  
  # ----------------------------------------------------------
  # covariance block
  # ----------------------------------------------------------
  
  cov_base <-
    (g - 1) +
    g * p
  
  cov_start <-
    cov_base +
    (k - 1) * 3 +
    1
  
  cov_ind <- cov_start:(cov_start + 2)
  
  G <- 0.5 *
    Sinv %*%
    (
      tcrossprod(d) -
        pars$Sigma[, , k]
    ) %*%
    Sinv
  
  # derivative with respect to
  # (sigma11, sigma21, sigma22)
  
  out[cov_ind] <- c(
    G[1, 1],
    2 * G[2, 1],
    G[2, 2]
  )
  
  out
}


# ============================================================
# 4. J(y), FEATURE SCORE, AND INFORMATION
# ============================================================

point_information <- function(y, pars) {
  
  tau <- as.vector(
    posterior_prob(
      matrix(y, nrow = 1),
      pars
    )
  )
  
  Tmat <- sapply(
    1:g,
    function(k)
      class_score(y, k, pars)
  )
  
  tbar <- as.vector(
    Tmat %*% tau
  )
  
  J <- matrix(
    0,
    qdim,
    qdim
  )
  
  for (k in 1:g) {
    
    u <- Tmat[, k] - tbar
    
    J <- J +
      tau[k] * tcrossprod(u)
  }
  
  list(
    J = J,
    tbar = tbar,
    tau = tau,
    T = Tmat
  )
}


information_objects <- function(Y, pars) {
  
  n <- nrow(Y)
  
  Jlist <- vector(
    "list",
    n
  )
  
  feature_info <- matrix(
    0,
    qdim,
    qdim
  )
  
  for (i in 1:n) {
    
    obj <- point_information(
      Y[i, ],
      pars
    )
    
    Jlist[[i]] <- obj$J
    
    feature_info <-
      feature_info +
      tcrossprod(obj$tbar)
  }
  
  feature_info <- feature_info / n
  
  list(
    J = Jlist,
    IY = feature_info
  )
}


# ============================================================
# 5. ACTIVE QDA BAYES FACES AND H_R
# ============================================================

compute_H <- function(
    pars,
    ngrid = 220,
    return_segments = FALSE) {
  
  sds <- matrix(
    NA_real_,
    g,
    2
  )
  
  for (k in 1:g) {
    sds[k, ] <- sqrt(
      diag(pars$Sigma[, , k])
    )
  }
  
  xlim <- c(
    min(pars$mu[, 1] - 4.5 * sds[, 1]),
    max(pars$mu[, 1] + 4.5 * sds[, 1])
  )
  
  ylim <- c(
    min(pars$mu[, 2] - 4.5 * sds[, 2]),
    max(pars$mu[, 2] + 4.5 * sds[, 2])
  )
  
  xx <- seq(
    xlim[1],
    xlim[2],
    length.out = ngrid
  )
  
  yy <- seq(
    ylim[1],
    ylim[2],
    length.out = ngrid
  )
  
  grid <- expand.grid(
    x = xx,
    y = yy
  )
  
  L <- log_weighted_density(
    as.matrix(
      grid[, c("x", "y"), drop = FALSE]
    ),
    pars
  )
  
  H <- matrix(
    0,
    qdim,
    qdim
  )
  
  seg_all <- list()
  
  seg_counter <- 0
  
  for (k in 1:(g - 1)) {
    
    for (l in (k + 1):g) {
      
      zz <- matrix(
        L[, k] - L[, l],
        nrow = length(xx),
        ncol = length(yy)
      )
      
      CL <- contourLines(
        x = xx,
        y = yy,
        z = zz,
        levels = 0
      )
      
      if (length(CL) == 0) next
      
      for (cl in CL) {
        
        if (length(cl$x) < 2) next
        
        for (s in 1:(length(cl$x) - 1)) {
          
          y1 <- c(
            cl$x[s],
            cl$y[s]
          )
          
          y2 <- c(
            cl$x[s + 1],
            cl$y[s + 1]
          )
          
          mid <- 0.5 * (
            y1 + y2
          )
          
          ds <- sqrt(
            sum(
              (y2 - y1)^2
            )
          )
          
          lm <- as.vector(
            log_weighted_density(
              matrix(mid, nrow = 1),
              pars
            )
          )
          
          other <- setdiff(
            1:g,
            c(k, l)
          )
          
          pair_level <- mean(
            lm[c(k, l)]
          )
          
          # active face only
          active <-
            length(other) == 0 ||
            pair_level >= max(lm[other])
          
          if (!active) next
          
          Sinv_k <- solve(
            pars$Sigma[, , k]
          )
          
          Sinv_l <- solve(
            pars$Sigma[, , l]
          )
          
          dkl <-
            Sinv_l %*%
            (mid - pars$mu[l, ]) -
            Sinv_k %*%
            (mid - pars$mu[k, ])
          
          norm_d <- sqrt(
            sum(dkl^2)
          )
          
          if (norm_d < 1e-10) next
          
          tk <- class_score(
            mid,
            k,
            pars
          )
          
          tl <- class_score(
            mid,
            l,
            pars
          )
          
          dt <- tk - tl
          
          rkl <- exp(
            pair_level
          )
          
          H <- H +
            rkl *
            tcrossprod(dt) /
            norm_d *
            ds
          
          if (return_segments) {
            
            seg_counter <-
              seg_counter + 1
            
            seg_all[[seg_counter]] <-
              data.frame(
                x = y1[1],
                y = y1[2],
                xend = y2[1],
                yend = y2[2],
                pair = paste0(k, "-", l)
              )
          }
        }
      }
    }
  }
  
  H <- 0.5 * (
    H + t(H)
  )
  
  out <- list(
    H = H,
    xlim = xlim,
    ylim = ylim
  )
  
  if (return_segments) {
    
    out$segments <-
      if (length(seg_all) > 0)
        bind_rows(seg_all)
    else
      data.frame()
  }
  
  out
}


# ============================================================
# 6. CLASSIFICATION VALUE psi_a(y)
# ============================================================

psi_from_G <- function(
    Y,
    pars,
    Gmat) {
  
  Y <- as.matrix(Y)
  
  out <- numeric(
    nrow(Y)
  )
  
  for (i in 1:nrow(Y)) {
    
    obj <- point_information(
      Y[i, ],
      pars
    )
    
    val <- 0
    
    for (k in 1:g) {
      
      u <-
        obj$T[, k] -
        obj$tbar
      
      val <- val +
        obj$tau[k] *
        as.numeric(
          crossprod(
            u,
            Gmat %*% u
          )
        )
    }
    
    out[i] <- val
  }
  
  out
}


# ============================================================
# 7. PROJECTION ONTO
#       0 <= a_i <= 1,
#       sum a_i = B
# ============================================================

project_capped_simplex <- function(
    v,
    B,
    tol = 1e-10) {
  
  n <- length(v)
  
  B <- max(
    0,
    min(B, n)
  )
  
  if (B <= tol)
    return(rep(0, n))
  
  if (B >= n - tol)
    return(rep(1, n))
  
  lo <- min(v) - 1
  hi <- max(v) + 1
  
  for (iter in 1:100) {
    
    lambda <- 0.5 * (
      lo + hi
    )
    
    a <- pmin(
      1,
      pmax(
        0,
        v - lambda
      )
    )
    
    if (sum(a) > B) {
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


# ============================================================
# 8. CONVEX DESIGN SOLVER
#
# Minimize
#
# tr{ H [Ibase + scale sum_i a_i J_i]^{-1} }
#
# subject to
#
# 0 <= a_i <= 1
# sum a_i = B
# ============================================================

solve_design <- function(
    Jlist,
    Ibase,
    H,
    B,
    scale,
    maxit = 500,
    tol = 1e-8,
    verbose = FALSE) {
  
  N <- length(Jlist)
  
  B <- min(
    max(B, 0),
    N
  )
  
  a <- rep(
    B / N,
    N
  )
  
  build_I <- function(a) {
    
    M <- Ibase
    
    nz <- which(
      a > 1e-12
    )
    
    for (i in nz) {
      
      M <- M +
        scale *
        a[i] *
        Jlist[[i]]
    }
    
    M
  }
  
  
  objective <- function(a) {
    
    M <- build_I(a)
    
    invM <- solve(M)
    
    sum(
      H * t(invM)
    )
  }
  
  
  old_obj <- objective(a)
  
  for (iter in 1:maxit) {
    
    M <- build_I(a)
    
    invM <- solve(M)
    
    Gmat <-
      invM %*%
      H %*%
      invM
    
    psi <- vapply(
      Jlist,
      function(J)
        sum(
          Gmat * t(J)
        ),
      numeric(1)
    )
    
    # gradient is proportional to -psi
    grad <- -psi
    
    step <- 1
    
    accepted <- FALSE
    
    for (bt in 1:30) {
      
      anew <-
        project_capped_simplex(
          a - step * grad,
          B
        )
      
      new_obj <- objective(anew)
      
      if (
        is.finite(new_obj) &&
        new_obj <= old_obj + 1e-12
      ) {
        
        accepted <- TRUE
        
        break
      }
      
      step <- step / 2
    }
    
    if (!accepted)
      break
    
    rel <-
      abs(
        old_obj - new_obj
      ) /
      max(
        1,
        abs(old_obj)
      )
    
    a <- anew
    
    old_obj <- new_obj
    
    if (
      verbose &&
      iter %% 25 == 0
    ) {
      
      cat(
        "iter =",
        iter,
        " objective =",
        old_obj,
        " rel =",
        rel,
        "\n"
      )
    }
    
    if (rel < tol)
      break
  }
  
  M <- build_I(a)
  
  invM <- solve(M)
  
  Gmat <-
    invM %*%
    H %*%
    invM
  
  psi <- vapply(
    Jlist,
    function(J)
      sum(
        Gmat * t(J)
      ),
    numeric(1)
  )
  
  list(
    a = a,
    psi = psi,
    objective = old_obj,
    iterations = iter
  )
}


# ============================================================
# 9. SEMI-SUPERVISED QDA EM
# ============================================================

initial_from_labels <- function(
    Y,
    labels,
    g,
    ridge = 1e-5) {
  
  n <- nrow(Y)
  
  global_mu <- colMeans(Y)
  
  global_S <-
    cov(Y) +
    diag(ridge, ncol(Y))
  
  pi_hat <- numeric(g)
  
  mu_hat <- matrix(
    NA_real_,
    g,
    ncol(Y)
  )
  
  S_hat <- array(
    NA_real_,
    dim = c(
      ncol(Y),
      ncol(Y),
      g
    )
  )
  
  for (k in 1:g) {
    
    ind <- which(
      labels == k
    )
    
    pi_hat[k] <-
      (length(ind) + 1) /
      (
        sum(!is.na(labels)) +
          g
      )
    
    if (length(ind) >= 3) {
      
      mu_hat[k, ] <-
        colMeans(
          Y[ind, , drop = FALSE]
        )
      
      Sk <- cov(
        Y[ind, , drop = FALSE]
      )
      
      S_hat[, , k] <-
        Sk +
        diag(
          ridge,
          ncol(Y)
        )
      
    } else {
      
      # fallback only if a pilot class is very sparse
      mu_hat[k, ] <-
        global_mu +
        rnorm(
          ncol(Y),
          sd = 0.10
        )
      
      S_hat[, , k] <-
        global_S
    }
  }
  
  pi_hat <-
    pi_hat /
    sum(pi_hat)
  
  list(
    pi = pi_hat,
    mu = mu_hat,
    Sigma = S_hat
  )
}


observed_loglik <- function(
    Y,
    labels,
    pars) {
  
  L <- log_weighted_density(
    Y,
    pars
  )
  
  ans <- 0
  
  lab <- which(
    !is.na(labels)
  )
  
  unlab <- which(
    is.na(labels)
  )
  
  if (length(lab) > 0) {
    
    ans <- ans +
      sum(
        L[
          cbind(
            lab,
            labels[lab]
          )
        ]
      )
  }
  
  if (length(unlab) > 0) {
    
    Lu <- L[
      unlab,
      ,
      drop = FALSE
    ]
    
    m <- apply(
      Lu,
      1,
      max
    )
    
    ans <- ans +
      sum(
        m +
          log(
            rowSums(
              exp(
                Lu - m
              )
            )
          )
      )
  }
  
  ans
}


fit_ss_qda <- function(
    Y,
    labels,
    init,
    maxit = 500,
    tol = 1e-7,
    ridge = 1e-5) {
  
  pars <- init
  
  n <- nrow(Y)
  
  old_ll <- -Inf
  
  for (iter in 1:maxit) {
    
    tau <- posterior_prob(
      Y,
      pars
    )
    
    lab <- which(
      !is.na(labels)
    )
    
    if (length(lab) > 0) {
      
      tau[lab, ] <- 0
      
      tau[
        cbind(
          lab,
          labels[lab]
        )
      ] <- 1
    }
    
    nk <- colSums(tau)
    
    pi_new <- nk / n
    
    mu_new <- matrix(
      NA_real_,
      g,
      p
    )
    
    S_new <- array(
      NA_real_,
      dim = c(p, p, g)
    )
    
    for (k in 1:g) {
      
      wk <- tau[, k]
      
      mu_new[k, ] <-
        colSums(
          Y * wk
        ) /
        nk[k]
      
      D <-
        sweep(
          Y,
          2,
          mu_new[k, ],
          "-"
        )
      
      S_new[, , k] <-
        crossprod(
          D,
          D * wk
        ) /
        nk[k]
      
      S_new[, , k] <-
        0.5 *
        (
          S_new[, , k] +
            t(S_new[, , k])
        ) +
        diag(ridge, p)
    }
    
    new_pars <- list(
      pi = pi_new,
      mu = mu_new,
      Sigma = S_new
    )
    
    ll <- observed_loglik(
      Y,
      labels,
      new_pars
    )
    
    if (
      is.finite(old_ll) &&
      abs(ll - old_ll) /
      max(1, abs(old_ll)) <
      tol
    ) {
      
      pars <- new_pars
      
      break
    }
    
    pars <- new_pars
    
    old_ll <- ll
  }
  
  list(
    pars = pars,
    logLik = ll,
    iterations = iter
  )
}


# ============================================================
# 10. ACQUISITION RULES AND DESIGN DIAGNOSTICS
# ============================================================

entropy_score <- function(tau) {
  
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


margin_score <- function(tau) {
  
  apply(
    tau,
    1,
    function(x) {
      
      sx <- sort(
        x,
        decreasing = TRUE
      )
      
      # larger score = more uncertain
      -(sx[1] - sx[2])
    }
  )
}


top_B <- function(
    score,
    B) {
  
  if (B <= 0) {
    return(integer(0))
  }
  
  order(
    score,
    decreasing = TRUE
  )[seq_len(B)]
}


# ------------------------------------------------------------
# Evaluate a binary finite-pool design
#
# This is used to compare the relaxed optimum returned by
# solve_design() with the actual top-B binary acquisition rule.
# ------------------------------------------------------------

binary_design_objective <- function(
    selected,
    Jlist,
    Ibase,
    H,
    scale) {
  
  M <- Ibase
  
  if (length(selected) > 0) {
    
    for (i in selected) {
      
      M <-
        M +
        scale *
        Jlist[[i]]
    }
  }
  
  invM <- solve(M)
  
  sum(
    H * t(invM)
  )
}


# ------------------------------------------------------------
# Round the relaxed design to exactly B labels and quantify
# the induced criterion loss.
# ------------------------------------------------------------

round_design_topB <- function(
    design,
    Jlist,
    Ibase,
    H,
    B,
    scale,
    fractional_tol = 1e-6) {
  
  N <- length(design$a)
  
  B <- min(
    max(as.integer(B), 0L),
    N
  )
  
  if (B == 0L) {
    
    selected <- integer(0)
    
  } else {
    
    selected <-
      order(
        design$a,
        decreasing = TRUE
      )[seq_len(B)]
  }
  
  rounded_objective <-
    binary_design_objective(
      selected = selected,
      Jlist = Jlist,
      Ibase = Ibase,
      H = H,
      scale = scale
    )
  
  fractional <-
    design$a > fractional_tol &
    design$a < 1 - fractional_tol
  
  fractional_count <-
    sum(fractional)
  
  fractional_prop <-
    fractional_count / N
  
  denom <-
    max(
      abs(design$objective),
      1e-15
    )
  
  rounding_loss_pct <-
    100 *
    (
      rounded_objective -
        design$objective
    ) /
    denom
  
  list(
    selected = selected,
    relaxed_objective = design$objective,
    rounded_objective = rounded_objective,
    rounding_loss_pct = rounding_loss_pct,
    fractional_count = fractional_count,
    fractional_prop = fractional_prop
  )
}


# ============================================================
# 11. POPULATION GEOMETRY FIGURE
# ============================================================

make_geometry_figure <- function(
    pars,
    rho = 0.20,
    Npool = 5000,
    ngrid = 180) {
  
  cat(
    "Creating geometry figure...\n"
  )
  
  Hobj <- compute_H(
    pars,
    ngrid = 240,
    return_segments = TRUE
  )
  
  H <- Hobj$H
  
  xx <- seq(
    Hobj$xlim[1],
    Hobj$xlim[2],
    length.out = ngrid
  )
  
  yy <- seq(
    Hobj$ylim[1],
    Hobj$ylim[2],
    length.out = ngrid
  )
  
  grid <- expand.grid(
    x = xx,
    y = yy
  )
  
  grid_xy <-
    as.matrix(
      grid[, c("x", "y"), drop = FALSE]
    )
  
  tau_grid <-
    posterior_prob(
      grid_xy,
      pars
    )
  
  grid$class <-
    factor(
      max.col(
        tau_grid,
        ties.method = "first"
      )
    )
  
  grid$entropy <-
    entropy_score(
      tau_grid
    ) /
    log(g)
  
  
  # ----------------------------------------------------------
  # Approximate population information under uniform design
  # ----------------------------------------------------------
  
  pool <-
    simulate_qda(
      Npool,
      pars
    )$Y
  
  info <-
    information_objects(
      pool,
      pars
    )
  
  I_uniform <- info$IY
  
  for (i in seq_len(Npool)) {
    
    I_uniform <-
      I_uniform +
      rho / Npool *
      info$J[[i]]
  }
  
  invI <- solve(
    I_uniform
  )
  
  G_uniform <-
    invI %*%
    H %*%
    invI
  
  psi_grid <-
    psi_from_G(
      grid_xy,
      pars,
      G_uniform
    )
  
  # Visualization only:
  # cap the extreme 1% so that the interior structure remains
  # visible without changing the acquisition calculations.
  
  cap <-
    as.numeric(
      quantile(
        psi_grid,
        0.99,
        na.rm = TRUE
      )
    )
  
  if (
    !is.finite(cap) ||
    cap <= 0
  ) {
    cap <- max(
      psi_grid,
      na.rm = TRUE
    )
  }
  
  grid$risk_value <-
    pmin(
      psi_grid / cap,
      1
    )
  
  
  # ----------------------------------------------------------
  # Compare population acquisition locations
  # ----------------------------------------------------------
  
  tau_pool <-
    posterior_prob(
      pool,
      pars
    )
  
  entropy_pool <-
    entropy_score(
      tau_pool
    )
  
  B <-
    round(
      rho * Npool
    )
  
  entropy_ind <-
    top_B(
      entropy_pool,
      B
    )
  
  
  # ----------------------------------------------------------
  # Classification-risk-optimal population design
  # ----------------------------------------------------------
  
  risk_design <-
    solve_design(
      Jlist = info$J,
      Ibase = info$IY,
      H = H,
      B = B,
      scale = 1 / Npool,
      maxit = 400,
      tol = 1e-9
    )
  
  risk_round <-
    round_design_topB(
      design = risk_design,
      Jlist = info$J,
      Ibase = info$IY,
      H = H,
      B = B,
      scale = 1 / Npool
    )
  
  risk_ind <-
    risk_round$selected
  
  
  cat(
    "Population risk-design rounding loss =",
    round(
      risk_round$rounding_loss_pct,
      4
    ),
    "%\n"
  )
  
  cat(
    "Population fractional weights =",
    risk_round$fractional_count,
    "of",
    Npool,
    "\n"
  )
  
  
  sel_entropy <-
    data.frame(
      x = pool[entropy_ind, 1],
      y = pool[entropy_ind, 2],
      method = "Entropy"
    )
  
  sel_risk <-
    data.frame(
      x = pool[risk_ind, 1],
      y = pool[risk_ind, 2],
      method =
        "Classification-risk optimal"
    )
  
  selected <-
    bind_rows(
      sel_entropy,
      sel_risk
    )
  
  background <-
    data.frame(
      x = pool[, 1],
      y = pool[, 2]
    )
  
  boundary <- Hobj$segments
  
  
  # ----------------------------------------------------------
  # Panel A
  # ----------------------------------------------------------
  
  p1 <-
    ggplot(
      grid,
      aes(x, y)
    ) +
    geom_raster(
      aes(fill = class),
      alpha = 0.58
    ) +
    geom_segment(
      data = boundary,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    scale_fill_manual(
      values = c(
        "#8ecae6",
        "#ffb703",
        "#b8de8f"
      )
    ) +
    coord_equal() +
    labs(
      title =
        "(a) Bayes classification geometry",
      x = expression(y[1]),
      y = expression(y[2]),
      fill = "Class"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  
  # ----------------------------------------------------------
  # Panel B
  # ----------------------------------------------------------
  
  p2 <-
    ggplot(
      grid,
      aes(x, y)
    ) +
    geom_raster(
      aes(fill = entropy)
    ) +
    geom_segment(
      data = boundary,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      linewidth = 0.45
    ) +
    scale_fill_viridis_c(
      option = "magma"
    ) +
    coord_equal() +
    labs(
      title =
        "(b) Posterior uncertainty",
      x = expression(y[1]),
      y = expression(y[2]),
      fill = "Entropy"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  
  # ----------------------------------------------------------
  # Panel C
  # ----------------------------------------------------------
  
  p3 <-
    ggplot(
      grid,
      aes(x, y)
    ) +
    geom_raster(
      aes(fill = risk_value)
    ) +
    geom_segment(
      data = boundary,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      linewidth = 0.45
    ) +
    scale_fill_viridis_c(
      option = "inferno"
    ) +
    coord_equal() +
    labs(
      title =
        "(c) Classification value",
      subtitle =
        expression(psi[a](y)),
      x = expression(y[1]),
      y = expression(y[2]),
      fill =
        "Relative value"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  
  # ----------------------------------------------------------
  # Panel D
  # ----------------------------------------------------------
  
  p4 <-
    ggplot() +
    geom_point(
      data = background,
      aes(x, y),
      alpha = 0.075,
      size = 0.55
    ) +
    geom_point(
      data = selected,
      aes(x, y),
      alpha = 0.48,
      size = 0.9
    ) +
    geom_segment(
      data = boundary,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linewidth = 0.45
    ) +
    facet_wrap(
      ~ method,
      nrow = 1
    ) +
    coord_equal() +
    labs(
      title =
        paste0(
          "(d) Labels selected under a ",
          round(100 * rho),
          "% budget"
        ),
      x = expression(y[1]),
      y = expression(y[2])
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      panel.grid = element_blank()
    )
  
  
  fig <-
    (p1 | p2) /
    (p3 | p4) +
    plot_annotation(
      title =
        "Classification uncertainty and classification value are different"
    )
  
  
  ggsave(
    "Figure1_geometry.png",
    fig,
    width = 12.5,
    height = 10,
    dpi = 400
  )
  
  fig
}


# ============================================================
# 12. MONTE CARLO STUDY
# ============================================================

run_simulation <- function(
    R = 100,
    n = 1000,
    pilot_n = 60,
    budgets = c(
      0.10,
      0.20,
      0.30
    ),
    test_n = 100000,
    H_ngrid = 150) {
  
  methods <- c(
    "Random",
    "Entropy",
    "Margin",
    "Fisher",
    "Adaptive risk-optimal",
    "Oracle risk-optimal"
  )
  
  
  # ----------------------------------------------------------
  # True classification-risk curvature:
  # computed once because it does not change across replications
  # ----------------------------------------------------------
  
  cat(
    "Computing true classification-risk curvature...\n"
  )
  
  Htrue <-
    compute_H(
      true_par,
      ngrid = H_ngrid
    )$H
  
  
  # ----------------------------------------------------------
  # Common large test set
  # ----------------------------------------------------------
  
  test <-
    simulate_qda(
      test_n,
      true_par
    )
  
  bayes_pred <-
    predict_qda(
      test$Y,
      true_par
    )
  
  bayes_error <-
    mean(
      bayes_pred != test$z
    )
  
  cat(
    "Bayes error =",
    round(
      bayes_error,
      5
    ),
    "\n"
  )
  
  
  ans <- list()
  diagnostic_ans <- list()
  
  row_id <- 0
  diagnostic_id <- 0
  
  
  for (r in seq_len(R)) {
    
    cat(
      "\nReplication",
      r,
      "of",
      R,
      "\n"
    )
    
    
    # --------------------------------------------------------
    # Generate training sample
    # --------------------------------------------------------
    
    dat <-
      simulate_qda(
        n,
        true_par
      )
    
    Y <- dat$Y
    z <- dat$z
    
    
    # --------------------------------------------------------
    # Common random pilot for every acquisition method
    # --------------------------------------------------------
    
    pilot <-
      sample(
        seq_len(n),
        pilot_n,
        replace = FALSE
      )
    
    labels_pilot <-
      rep(
        NA_integer_,
        n
      )
    
    labels_pilot[pilot] <-
      z[pilot]
    
    
    init <-
      initial_from_labels(
        Y,
        labels_pilot,
        g
      )
    
    pilot_fit <-
      tryCatch(
        fit_ss_qda(
          Y,
          labels_pilot,
          init
        ),
        error = function(e)
          NULL
      )
    
    if (is.null(pilot_fit)) {
      
      cat(
        "Pilot fit failed; replication skipped.\n"
      )
      
      next
    }
    
    pfit <-
      pilot_fit$pars
    
    
    # --------------------------------------------------------
    # Plug-in quantities
    # --------------------------------------------------------
    
    tau <-
      posterior_prob(
        Y,
        pfit
      )
    
    entropy <-
      entropy_score(
        tau
      )
    
    margin <-
      margin_score(
        tau
      )
    
    
    info <-
      information_objects(
        Y,
        pfit
      )
    
    
    # --------------------------------------------------------
    # Estimated I_Y plus information from pilot labels
    # --------------------------------------------------------
    
    Ibase <- info$IY
    
    for (i in pilot) {
      
      Ibase <-
        Ibase +
        info$J[[i]] / n
    }
    
    
    candidate <-
      setdiff(
        seq_len(n),
        pilot
      )
    
    Jcand <-
      info$J[candidate]
    
    
    # --------------------------------------------------------
    # Estimated complete-classification information
    #
    # I_CC = I_Y + E{J(Y)}
    #
    # This defines the Fisher-information comparator:
    #
    # tr{ I_CC I(a)^(-1) }
    #
    # Unlike ordinary A-optimality with H = I, this criterion
    # is invariant to smooth one-to-one reparameterization.
    # --------------------------------------------------------
    
    ICC_fit <- info$IY
    
    for (i in seq_len(n)) {
      
      ICC_fit <-
        ICC_fit +
        info$J[[i]] / n
    }
    
    
    # --------------------------------------------------------
    # Estimated classification-risk curvature
    # --------------------------------------------------------
    
    Hfit <-
      tryCatch(
        compute_H(
          pfit,
          ngrid = H_ngrid
        )$H,
        error = function(e)
          NULL
      )
    
    if (
      is.null(Hfit) ||
      !all(is.finite(Hfit)) ||
      sum(abs(Hfit)) < 1e-12
    ) {
      
      cat(
        "Estimated H computation failed; replication skipped.\n"
      )
      
      next
    }
    
    
    # --------------------------------------------------------
    # Oracle quantities
    #
    # Oracle means that the acquisition criterion knows the
    # true generating parameter. The classifier fitted after
    # acquisition is still estimated from the observed sample.
    # --------------------------------------------------------
    
    info_oracle <-
      information_objects(
        Y,
        true_par
      )
    
    Ibase_oracle <-
      info_oracle$IY
    
    for (i in pilot) {
      
      Ibase_oracle <-
        Ibase_oracle +
        info_oracle$J[[i]] / n
    }
    
    Jcand_oracle <-
      info_oracle$J[candidate]
    
    
    # --------------------------------------------------------
    # Each labeling budget
    # --------------------------------------------------------
    
    for (rho in budgets) {
      
      total_B <-
        round(
          rho * n
        )
      
      additional_B <-
        total_B -
        pilot_n
      
      if (additional_B <= 0) {
        
        stop(
          paste0(
            "pilot_n = ",
            pilot_n,
            " must be smaller than every total label budget."
          )
        )
      }
      
      additional_B <-
        min(
          additional_B,
          length(candidate)
        )
      
      
      # ======================================================
      # RANDOM
      # ======================================================
      
      sel_random <-
        sample(
          candidate,
          additional_B,
          replace = FALSE
        )
      
      
      # ======================================================
      # ENTROPY
      # ======================================================
      
      sel_entropy <-
        candidate[
          top_B(
            entropy[candidate],
            additional_B
          )
        ]
      
      
      # ======================================================
      # MARGIN
      # ======================================================
      
      sel_margin <-
        candidate[
          top_B(
            margin[candidate],
            additional_B
          )
        ]
      
      
      # ======================================================
      # FISHER-INFORMATION DESIGN
      #
      # min tr{ I_CC I(a)^(-1) }
      # ======================================================
      
      Fisher_design <-
        solve_design(
          Jlist = Jcand,
          Ibase = Ibase,
          H = ICC_fit,
          B = additional_B,
          scale = 1 / n,
          maxit = 350,
          tol = 1e-8
        )
      
      
      Fisher_round <-
        round_design_topB(
          design = Fisher_design,
          Jlist = Jcand,
          Ibase = Ibase,
          H = ICC_fit,
          B = additional_B,
          scale = 1 / n
        )
      
      
      sel_Fisher <-
        candidate[
          Fisher_round$selected
        ]
      
      
      # ======================================================
      # ADAPTIVE CLASSIFICATION-RISK-OPTIMAL DESIGN
      #
      # Uses pilot-fitted J(y), I_Y and H_R.
      # ======================================================
      
      Adaptive_design <-
        solve_design(
          Jlist = Jcand,
          Ibase = Ibase,
          H = Hfit,
          B = additional_B,
          scale = 1 / n,
          maxit = 350,
          tol = 1e-8
        )
      
      
      Adaptive_round <-
        round_design_topB(
          design = Adaptive_design,
          Jlist = Jcand,
          Ibase = Ibase,
          H = Hfit,
          B = additional_B,
          scale = 1 / n
        )
      
      
      sel_Adaptive <-
        candidate[
          Adaptive_round$selected
        ]
      
      
      # ======================================================
      # ORACLE CLASSIFICATION-RISK-OPTIMAL DESIGN
      #
      # Uses the true QDA parameter and true H_R.
      # This is a theoretical benchmark, not an implementable
      # acquisition strategy.
      # ======================================================
      
      Oracle_design <-
        solve_design(
          Jlist = Jcand_oracle,
          Ibase = Ibase_oracle,
          H = Htrue,
          B = additional_B,
          scale = 1 / n,
          maxit = 350,
          tol = 1e-8
        )
      
      
      Oracle_round <-
        round_design_topB(
          design = Oracle_design,
          Jlist = Jcand_oracle,
          Ibase = Ibase_oracle,
          H = Htrue,
          B = additional_B,
          scale = 1 / n
        )
      
      
      sel_Oracle <-
        candidate[
          Oracle_round$selected
        ]
      
      
      # ======================================================
      # SAVE DESIGN-ROUNDING DIAGNOSTICS
      # ======================================================
      
      design_diagnostics <-
        list(
          Fisher = Fisher_round,
          `Adaptive risk-optimal` =
            Adaptive_round,
          `Oracle risk-optimal` =
            Oracle_round
        )
      
      
      design_iterations <-
        c(
          Fisher =
            Fisher_design$iterations,
          `Adaptive risk-optimal` =
            Adaptive_design$iterations,
          `Oracle risk-optimal` =
            Oracle_design$iterations
        )
      
      
      for (
        design_name in
        names(design_diagnostics)
      ) {
        
        DD <-
          design_diagnostics[[design_name]]
        
        diagnostic_id <-
          diagnostic_id + 1
        
        
        diagnostic_ans[[diagnostic_id]] <-
          data.frame(
            replication = r,
            n = n,
            pilot_n = pilot_n,
            rho = rho,
            total_budget = total_B,
            additional_budget =
              additional_B,
            design = design_name,
            iterations =
              as.numeric(
                design_iterations[
                  design_name
                ]
              ),
            fractional_count =
              DD$fractional_count,
            fractional_prop =
              DD$fractional_prop,
            relaxed_objective =
              DD$relaxed_objective,
            rounded_objective =
              DD$rounded_objective,
            rounding_loss_pct =
              DD$rounding_loss_pct
          )
      }
      
      
      # ======================================================
      # SELECTION SETS
      # ======================================================
      
      selections <-
        list(
          Random =
            sel_random,
          Entropy =
            sel_entropy,
          Margin =
            sel_margin,
          Fisher =
            sel_Fisher,
          `Adaptive risk-optimal` =
            sel_Adaptive,
          `Oracle risk-optimal` =
            sel_Oracle
        )
      
      
      # ======================================================
      # REFIT CLASSIFIER FOR EACH ACQUISITION STRATEGY
      # ======================================================
      
      for (method in methods) {
        
        labels <-
          rep(
            NA_integer_,
            n
          )
        
        labelled <-
          c(
            pilot,
            selections[[method]]
          )
        
        labelled <-
          unique(
            labelled
          )
        
        labels[labelled] <-
          z[labelled]
        
        
        fit <-
          tryCatch(
            fit_ss_qda(
              Y,
              labels,
              pfit
            ),
            error = function(e)
              NULL
          )
        
        if (is.null(fit)) {
          
          cat(
            "Final fit failed:",
            method,
            "rho =",
            rho,
            "\n"
          )
          
          next
        }
        
        
        pred <-
          predict_qda(
            test$Y,
            fit$pars
          )
        
        
        err <-
          mean(
            pred != test$z
          )
        
        
        row_id <-
          row_id + 1
        
        
        ans[[row_id]] <-
          data.frame(
            replication = r,
            n = n,
            pilot_n = pilot_n,
            rho = rho,
            n_labelled =
              length(labelled),
            method = method,
            error = err,
            bayes_error =
              bayes_error,
            excess_error =
              err - bayes_error,
            logLik =
              fit$logLik,
            EM_iterations =
              fit$iterations
          )
      }
    }
    
    
    # --------------------------------------------------------
    # Save progress continuously
    # --------------------------------------------------------
    
    if (
      length(ans) > 0 &&
      r %% 5 == 0
    ) {
      
      write.csv(
        bind_rows(ans),
        "Study1_raw_PROGRESS.csv",
        row.names = FALSE
      )
    }
    
    
    if (
      length(diagnostic_ans) > 0 &&
      r %% 5 == 0
    ) {
      
      write.csv(
        bind_rows(
          diagnostic_ans
        ),
        "Study1_design_diagnostics_PROGRESS.csv",
        row.names = FALSE
      )
    }
  }
  
  
  # ==========================================================
  # RAW CLASSIFICATION RESULTS
  # ==========================================================
  
  raw <-
    bind_rows(
      ans
    )
  
  
  summary <-
    raw %>%
    group_by(
      rho,
      method
    ) %>%
    summarise(
      R = n(),
      mean_error =
        mean(
          error,
          na.rm = TRUE
        ),
      sd_error =
        sd(
          error,
          na.rm = TRUE
        ),
      se_error =
        sd_error /
        sqrt(R),
      mean_excess =
        mean(
          excess_error,
          na.rm = TRUE
        ),
      sd_excess =
        sd(
          excess_error,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
  
  
  write.csv(
    raw,
    "Study1_raw.csv",
    row.names = FALSE
  )
  
  
  write.csv(
    summary,
    "Study1_summary.csv",
    row.names = FALSE
  )
  
  
  # ==========================================================
  # DESIGN-ROUNDING DIAGNOSTICS
  # ==========================================================
  
  design_diagnostics <-
    bind_rows(
      diagnostic_ans
    )
  
  
  design_diagnostics_summary <-
    design_diagnostics %>%
    group_by(
      rho,
      design
    ) %>%
    summarise(
      R = n(),
      mean_fractional_count =
        mean(
          fractional_count
        ),
      mean_fractional_prop =
        mean(
          fractional_prop
        ),
      max_fractional_prop =
        max(
          fractional_prop
        ),
      mean_rounding_loss_pct =
        mean(
          rounding_loss_pct
        ),
      max_rounding_loss_pct =
        max(
          rounding_loss_pct
        ),
      mean_iterations =
        mean(
          iterations
        ),
      .groups = "drop"
    )
  
  
  write.csv(
    design_diagnostics,
    "Study1_design_diagnostics.csv",
    row.names = FALSE
  )
  
  
  write.csv(
    design_diagnostics_summary,
    "Study1_design_diagnostics_summary.csv",
    row.names = FALSE
  )
  
  
  # ==========================================================
  # PERFORMANCE FIGURE
  # ==========================================================
  
  summary_plot <-
    summary %>%
    mutate(
      method =
        factor(
          method,
          levels = methods
        )
    )
  
  
  fig2 <-
    ggplot(
      summary_plot,
      aes(
        x = rho,
        y = mean_error,
        group = method,
        linetype = method,
        shape = method
      )
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    geom_point(
      size = 2.5
    ) +
    geom_errorbar(
      aes(
        ymin =
          mean_error -
          1.96 * se_error,
        ymax =
          mean_error +
          1.96 * se_error
      ),
      width = 0.008,
      linewidth = 0.55
    ) +
    geom_hline(
      yintercept =
        bayes_error,
      linetype = "dotted",
      linewidth = 0.55
    ) +
    scale_x_continuous(
      breaks = budgets,
      labels =
        paste0(
          100 * budgets,
          "%"
        )
    ) +
    labs(
      x =
        "Proportion of acquired labels",
      y =
        "Test classification error",
      linetype =
        "Acquisition rule",
      shape =
        "Acquisition rule",
      title =
        "Classification error under a fixed labeling budget",
      subtitle =
        paste0(
          "Dotted horizontal line: Bayes error = ",
          sprintf(
            "%.4f",
            bayes_error
          )
        )
    ) +
    theme_bw(
      base_size = 13
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor =
        element_blank()
    )
  
  
  ggsave(
    "Figure2_performance.png",
    fig2,
    width = 9.5,
    height = 6.5,
    dpi = 400
  )
  
  
  list(
    raw =
      raw,
    summary =
      summary,
    design_diagnostics =
      design_diagnostics,
    design_diagnostics_summary =
      design_diagnostics_summary,
    bayes_error =
      bayes_error,
    figure =
      fig2
  )
}


# ============================================================
# 13. RUN
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "STEP 1: population geometry\n"
)

cat(
  "============================================\n"
)


fig1 <-
  make_geometry_figure(
    true_par,
    rho = 0.20,
    Npool = 5000,
    ngrid = 180
  )


print(
  fig1
)


cat(
  "\n============================================\n"
)

cat(
  "STEP 2: Monte Carlo diagnostic\n"
)

cat(
  "============================================\n"
)


# ------------------------------------------------------------
# IMPORTANT:
#
# Keep R = 10 for this diagnostic.
#
# We first inspect:
#
#   1. classification performance,
#   2. adaptive versus oracle risk-optimal,
#   3. Fisher comparator,
#   4. number of fractional design weights,
#   5. loss caused by top-B rounding.
#
# Only after these are satisfactory should R be increased.
# ------------------------------------------------------------

RES <-
  run_simulation(
    R = 10,
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


cat(
  "\n\n============================================\n"
)

cat(
  "CLASSIFICATION RESULTS\n"
)

cat(
  "============================================\n"
)


print(
  RES$summary
)


cat(
  "\n\n============================================\n"
)

cat(
  "DESIGN ROUNDING DIAGNOSTICS\n"
)

cat(
  "============================================\n"
)


print(
  RES$design_diagnostics_summary
)


# ============================================================
# DO NOT CHANGE R = 10 TO R = 100 YET.
#
# After this run, inspect:
#
#   Study1_summary.csv
#   Study1_design_diagnostics_summary.csv
#   Figure1_geometry.png
#   Figure2_performance.png
#
# Then decide whether the discrete rounding is sufficiently
# faithful and whether the adaptive/oracle gap is reasonable.
# ============================================================