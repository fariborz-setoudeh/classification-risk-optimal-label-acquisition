# Simulation convergence and precision validation for the reviewed acquisition paper.
# Source this file and call validation_09_simulation(...).
# Start with n_rep=2. The same frozen protocol can be extended to 100 via resume_dir.
# Original scripts, historical results, and manuscript files are not changed.

validation_09_simulation <- function(project_root, n_rep = 2L, resume_dir = NULL) {
  stopifnot(dir.exists(project_root), length(n_rep) == 1L, is.finite(n_rep),
    n_rep >= 1L, n_rep <= 100L, n_rep == as.integer(n_rep))
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  shadows <- c("pi", "sample", "sample.int", "set.seed", "cov", "sd", "eigen", "solve",
    "chol", "dnorm", "pnorm", "order", "max.col", "pmax", "rnorm", "integrate")
  shadows <- shadows[vapply(shadows, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(shadows)) stop("Global objects shadow numerical functions/constants: ", paste(shadows, collapse = ", "))
  files <- list.files(project_root, recursive = TRUE, full.names = TRUE)
  find_one <- function(name) {
    hits <- files[basename(files) == name]
    if (length(hits) != 1L) stop("Expected exactly one ", name, "; found ", length(hits), ".")
    normalizePath(hits, winslash = "/", mustWork = TRUE)
  }
  paths <- list(sim1 = find_one("01_simulation_scenario1.R"),
    sim2 = find_one("02_simulation_scenario2.R"), core = find_one("03_landsat_analysis.R"))
  hashes <- unname(tools::md5sum(unlist(paths)))
  if (!identical(hashes, c("fd4eb80059402b1984b69cac17d5050f", "7b9b34b51a9664005000e4cc78d78f97",
      "67b711404f884efec1a3f1d51340cbf5"))) stop("Use the reviewed, unedited original scripts.")
  hash_raw <- function(bytes) {
    path <- tempfile("prediction_digest_")
    on.exit(unlink(path), add = TRUE)
    writeBin(bytes, path)
    unname(tools::md5sum(path))
  }
  protocol <- list(version = "2026-10-04-simulation-joint-v1", hashes = hashes,
    implementation_md5 = hash_raw(charToRaw(paste(deparse(sys.function(), width.cutoff = 500L), collapse = "\n"))),
    RNG = c("Mersenne-Twister", "Inversion", "Rejection"), seed_base = 2026104000L,
    n = 1000L, pilot_n = 60L, test_n = 100000L, budgets = c(.1, .2, .3),
    EM_MAXIT = 20000L, EM_TOL = 1e-11, EM_PARAMETER_TOLS = c(1e-8, 1e-9, 1e-10, 1e-11),
    EM_STABLE_ITERATIONS = 5L, EM_NEWTON_TOL = 1e-6, COV_FLOOR_REL = 1e-6,
    FW_TOL = 1e-5, FW_MAXITS = c(2000L, 5000L),
    mc_reps = 8L, true_mc_n = 524288L, pilot_mc_sizes = c(131072L, 262144L, 524288L),
    max_H_MCSE = .01, max_white_H_MCSE = .01, max_objective_MCSE = .002,
    R_version = R.version.string)
  out <- if (is.null(resume_dir)) tempfile("simulation_validated_", tmpdir = project_root) else
    normalizePath(resume_dir, winslash = "/", mustWork = TRUE)
  if (is.null(resume_dir)) {
    if (!dir.create(out)) stop("Cannot create output directory.")
    dir.create(file.path(out, "fits")); saveRDS(protocol, file.path(out, "protocol.rds"))
  } else if (!identical(readRDS(file.path(out, "protocol.rds")), protocol)) {
    stop("Resume requires the same implementation, source files, protocol, and R version.")
  }
  sink(file.path(out, "simulation_09_log.txt"), append = !is.null(resume_dir), split = TRUE)
  initial_rng <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  initial_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(initial_rng)))
    if (had_seed) assign(".Random.seed", initial_seed, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    sink()
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  cat("Simulation joint-convergence diagnostic; ", format(Sys.time()), "\nOutput: ", out, "\n", sep = "")
  cat("Requested replications per scenario: ", n_rep, "; the manuscript requires 100.\n", sep = "")
  cat("Both scenarios, all three budgets, and all six methods are retained.\n")
  cat("New explicit training, pilot, random-acquisition, evaluation, and integration seeds are separated.\n")
  cat("The evaluation sample is fixed within each scenario and is not used to choose acquisition or precision.\n")
  cat("The M-step uses an eigenvalue floor, rather than adding a fixed ridge to every covariance.\n")
  cat("Likelihood/parameter stopping, observed-Hessian/Newton checks, and FW gap certificates are required.\n")
  cat("Curvature MCSEs and quadrature errors are diagnostics, not rigorous error bounds.\n")
  print(data.frame(path = unlist(paths), md5 = hashes)); print(protocol); print(sessionInfo())
  METHODS <- c("Random", "Entropy", "Margin", "Fisher", "AdaptiveRisk", "OracleRisk")
  records <- list(); current_checks <- list(); failures <- list(); scenario_records <- list()
  add <- function(ds, rep_id, stage, metric, value, method = "", budget = NA_real_, detail = "") {
    current_checks[[length(current_checks) + 1L]] <<- data.frame(Dataset = ds, Replication = rep_id,
      Stage = stage, Method = method, Budget = budget, Metric = metric, Value = as.character(value),
      Detail = detail, stringsAsFactors = FALSE)
  }
  atomic_save <- function(object, path) {
    temporary <- paste0(path, ".temporary")
    saveRDS(object, temporary)
    backup <- paste0(path, ".previous")
    existed <- file.exists(path)
    if (existed && !file.rename(path, backup)) stop("Cannot preserve old checkpoint: ", path)
    if (!file.rename(temporary, path)) {
      if (existed) file.rename(backup, path)
      stop("Cannot save checkpoint: ", path)
    }
    if (existed) unlink(backup)
  }

  definitions <- function(path, constants) {
    e <- new.env(parent = globalenv()); list2env(constants, e)
    for (expr in parse(file = path, keep.source = FALSE)) {
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
      rhs <- expr[[3L]]
      if (is.symbol(expr[[2L]]) && is.call(rhs) && identical(rhs[[1L]], as.name("function"))) eval(expr, e)
    }
    e
  }

  sym <- function(x) (x + t(x)) / 2

  fnorm <- function(x) sqrt(sum(x^2))

  inverse <- function(x) chol2inv(chol(sym(x)))

  relative <- function(x, y) fnorm(x - y) / max(fnorm(y), 1e-12)

  row_max <- function(M) do.call(pmax, lapply(seq_len(ncol(M)), function(k) M[, k]))

  replace_once <- function(text, old, new) {
    positions <- gregexpr(old, text, fixed = TRUE)[[1L]]
    if (positions[1L] < 0L || length(positions) != 1L) {
      stop("Expected one source marker: ", old)
    }
    sub(old, new, text, fixed = TRUE)
  }

  repaired_EM <- function(original) {
    src <- paste(deparse(original, width.cutoff = 500L), collapse = "\n")
    src <- replace_once(src, "unl_idx <- which(!labelled)", paste(
      "fit_new <- list(pi = pi_new, mu = mu_new, Sigma = Sigma_new, classes = classes)",
      "log_joint <- compute_log_joint(X, fit_new)",
      "unl_idx <- which(!labelled)", sep = "\n"))
    src <- replace_once(src, "ll <- ll_unl + ll_lab", paste(
      "ll <- ll_unl + ll_lab",
      "ll_history <- c(ll_history, ll)",
      "old_parameters <- c(pi_k, mu_k, unlist(Sigma_k))",
      "new_parameters <- c(pi_new, mu_new, unlist(Sigma_new))",
      "parameter_change <- max(abs(new_parameters - old_parameters) / (1 + abs(old_parameters)))",
      sep = "\n"))
    src <- replace_once(src, "old_ll <- -Inf", paste(
      "if (!is.null(start_fit)) {",
      "  pi_k <- start_fit$pi; mu_k <- start_fit$mu; Sigma_k <- start_fit$Sigma",
      "}", "old_ll <- -Inf", sep = "\n"))
    src <- replace_once(src, "if (rel_change < tol)", paste(
      "stable_joint_iterations <- if (rel_change < tol && parameter_change < parameter_tol) stable_joint_iterations + 1L else 0L",
      "if (stable_joint_iterations >= stable_required)", sep = "\n"))
    src <- replace_once(src, "break", "{ converged <- TRUE; break }")
    f <- eval(parse(text = src), envir = environment(original))
    formals(f) <- as.pairlist(c(as.list(formals(f)),
      list(start_fit = NULL, parameter_tol = 1e-8, stable_required = 5L)))
    b <- as.list(body(f))
    prefix <- quote({
      converged <- FALSE
      ll_history <- if (is.null(start_fit)) numeric(0) else start_fit$logLik_history
      stable_joint_iterations <- 0L
      parameter_change <- Inf
      rel_change <- Inf
    })
    last <- b[[length(b)]]
    b[[length(b)]] <- substitute({
      answer <- LAST
      answer$converged <- converged
      answer$stopping_reason <- if (converged) "joint_likelihood_parameter_tolerance" else "iteration_limit"
      answer$additional_iterations <- answer$iterations
      answer$iterations <- answer$iterations + if (is.null(start_fit)) 0L else start_fit$iterations
      answer$stable_joint_iterations <- stable_joint_iterations
      answer$parameter_tolerance <- parameter_tol
      answer$relative_logLik_change <- rel_change
      answer$max_scaled_parameter_change <- parameter_change
      answer$logLik_history <- ll_history
      answer
    }, list(LAST = last))
    body(f) <- as.call(c(list(as.name("{"), prefix), b[-1L]))
    f
  }

  make_model <- function(par, simulation_order = FALSE) {
    g <- length(par$pi); d <- ncol(par$mu)
    S <- if (is.list(par$Sigma)) par$Sigma else lapply(seq_len(g), function(k) par$Sigma[, , k])
    stopifnot(g >= 2L, d >= 2L, nrow(par$mu) == g, length(S) == g,
              all(par$pi > 0), abs(sum(par$pi) - 1) < 1e-10)
    P <- lapply(S, inverse); b <- t(vapply(seq_len(g), function(k) as.vector(P[[k]] %*% par$mu[k, ]), numeric(d)))
    cc <- vapply(seq_len(g), function(k) log(par$pi[k]) - d * log(2 * pi) / 2 -
      sum(log(diag(chol(S[[k]])))) - sum(par$mu[k, ] * b[k, ]) / 2, numeric(1))
    index <- cores$Landsat$build_param_index(g, d)
    if (simulation_order) {
      nc <- nrow(index$vech_pairs)
      for (k in seq_len(g)) {
        index$mu[[k]] <- (g - 1L) + d * (k - 1L) + seq_len(d)
        index$cov[[k]] <- (g - 1L) + g * d + nc * (k - 1L) + seq_len(nc)
        index$class[[k]] <- c(index$mu[[k]], index$cov[[k]])
      }
    }
    marginal <- lapply(seq_len(d), function(j) {
      other <- setdiff(seq_len(d), j)
      Sm <- lapply(S, function(s) s[other, other, drop = FALSE])
      list(other = other, P = lapply(Sm, inverse), R = lapply(Sm, chol),
        cc = vapply(seq_len(g), function(k) log(par$pi[k]) - (d - 1L) * log(2 * pi) / 2 -
          sum(log(diag(chol(Sm[[k]])))), numeric(1)))
    })
    list(pi = par$pi, mu = par$mu, S = S, P = P, b = b, cc = cc,
         g = g, d = d, p = index$p, index = index, marginal = marginal)
  }

  log_joint <- function(X, m) {
    X <- as.matrix(X)
    z <- vapply(seq_len(m$g), function(k) -rowSums((X %*% m$P[[k]]) * X) / 2 +
      as.vector(X %*% m$b[k, ]) + m$cc[k], numeric(nrow(X)))
    matrix(z, nrow(X), m$g)
  }

  block_score <- function(X, k, m) {
    u <- sweep(as.matrix(X), 2L, m$mu[k, ], "-") %*% m$P[[k]]
    pairs <- m$index$vech_pairs
    cs <- vapply(seq_len(nrow(pairs)), function(j) {
      a <- pairs[j, 1L]; b <- pairs[j, 2L]
      (u[, a] * u[, b] - m$P[[k]][a, b]) * if (a == b) 0.5 else 1
    }, numeric(nrow(u)))
    cbind(u, matrix(cs, nrow(u), nrow(pairs)))
  }

  score_object <- function(X, m) {
    lj <- log_joint(X, m); z <- sweep(lj, 1L, row_max(lj), "-")
    tau <- exp(z); tau <- tau / rowSums(tau)
    A <- matrix(-m$pi[seq_len(m$g - 1L)], m$g, m$g - 1L, byrow = TRUE)
    for (k in seq_len(m$g - 1L)) A[k, k] <- A[k, k] + 1
    ss <- lapply(seq_len(m$g), function(k) block_score(X, k, m))
    bar <- matrix(0, nrow(X), m$p); bar[, m$index$alpha] <- tau %*% A
    for (k in seq_len(m$g)) bar[, m$index$class[[k]]] <- ss[[k]] * tau[, k]
    list(tau = tau, log_joint = lj, A = A, class_scores = ss, bar = bar, index = m$index)
  }

  sample_marginal <- function(n, j, m, fixture = FALSE) {
    mm <- m$marginal[[j]]; half <- n %/% 2L
    if (fixture) {
      comp <- (seq_len(half) - 1L) %% m$g + 1L
      z <- vapply(seq_len(m$d - 1L), function(col) {
        qnorm(((seq_len(half) * (2L * col + 1L) + 7L * j) %% 67L + 0.5) / 67)
      }, numeric(half))
      z <- matrix(z, half, m$d - 1L)
    } else {
      comp <- sample.int(m$g, half, replace = TRUE, prob = m$pi)
      z <- matrix(rnorm(half * (m$d - 1L)), half, m$d - 1L)
    }
    X <- matrix(0, n, m$d - 1L)
    for (k in seq_len(m$g)) {
      ii <- which(comp == k)
      if (!length(ii)) next
      zz <- z[ii, , drop = FALSE] %*% mm$R[[k]]
      X[ii, ] <- sweep(zz, 2L, m$mu[k, mm$other], "+")
      X[ii + half, ] <- sweep(-zz, 2L, m$mu[k, mm$other], "+")
    }
    X
  }

  marginal_log_density <- function(X, j, m) {
    mm <- m$marginal[[j]]
    lj <- vapply(seq_len(m$g), function(k) {
      xc <- sweep(X, 2L, m$mu[k, mm$other], "-")
      mm$cc[k] - rowSums((xc %*% mm$P[[k]]) * xc) / 2
    }, numeric(nrow(X)))
    lj <- matrix(lj, nrow(X), m$g); mx <- row_max(lj)
    mx + log(rowSums(exp(sweep(lj, 1L, mx, "-"))))
  }

  atlas <- function(m, n, seed, chart_weights = rep(1, m$d), fixture = FALSE) {
    stopifnot(n %% 2L == 0L, length(chart_weights) == m$d, all(chart_weights > 0))
    set.seed(seed); H <- matrix(0, m$p, m$p); dd <- list(); rr <- 1L
    max_residual <- 0; max_scaled_residual <- 0
    for (j in seq_len(m$d)) {
      other <- m$marginal[[j]]$other; Xm <- sample_marginal(n, j, m, fixture)
      lp <- marginal_log_density(Xm, j, m)
      X <- matrix(0, n, m$d); X[, other] <- Xm
      for (k in seq_len(m$g - 1L)) for (l in (k + 1L):m$g) {
        A <- (m$P[[l]] - m$P[[k]]) / 2; b <- m$b[k, ] - m$b[l, ]; cc <- m$cc[k] - m$cc[l]
        aa <- A[j, j]; bb <- as.vector(2 * Xm %*% A[other, j, drop = FALSE]) + b[j]
        cv <- rowSums((Xm %*% A[other, other, drop = FALSE]) * Xm) + as.vector(Xm %*% b[other]) + cc
        roots <- list()
        if (aa == 0) {
          ii <- which(bb != 0)
          if (length(ii)) roots[[1L]] <- list(i = ii, z = -cv[ii] / bb[ii])
        } else {
          disc <- bb^2 - 4 * aa * cv; ii <- which(disc > 0)
          if (length(ii)) {
            zz <- -0.5 * (bb[ii] + ifelse(bb[ii] >= 0, 1, -1) * sqrt(disc[ii]))
            roots <- list(list(i = ii, z = zz / aa), list(i = ii, z = cv[ii] / zz))
          }
        }
        pair_H <- matrix(0, m$p, m$p); nr <- 0L; wm <- 0
        for (root in roots) {
          XX <- X[root$i, , drop = FALSE]; XX[, j] <- root$z
          lj <- log_joint(XX, m)
          active <- (lj[, k] + lj[, l]) / 2 >= row_max(lj) - 1e-10
          if (!any(active)) next
          XX <- XX[active, , drop = FALSE]; ii <- root$i[active]; lj <- lj[active, , drop = FALSE]
          grad <- sweep(2 * XX %*% A, 2L, b, "+")
          denom <- as.vector((grad^2) %*% chart_weights)
          if (any(denom <= 0)) stop("A sampled boundary is not regular.")
          w <- exp((lj[, k] + lj[, l]) / 2 - lp[ii]) * chart_weights[j] * abs(grad[, j]) / denom / n
          if (any(!is.finite(w))) stop("Nonfinite boundary importance weight.")
          prior <- matrix(0, nrow(XX), m$g - 1L)
          if (k < m$g) prior[, k] <- 1
          if (l < m$g) prior[, l] <- -1
          D <- cbind(prior, block_score(XX, k, m), -block_score(XX, l, m))
          ix <- c(m$index$alpha, m$index$class[[k]], m$index$class[[l]])
          update <- crossprod(D * sqrt(w))
          pair_H[ix, ix] <- pair_H[ix, ix, drop = FALSE] + update
          nr <- nr + nrow(XX); wm <- wm + sum(w)
          res <- abs(lj[, k] - lj[, l])
          max_residual <- max(max_residual, res)
          max_scaled_residual <- max(max_scaled_residual, res / (1 + abs(lj[, k]) + abs(lj[, l])))
        }
        H <- H + pair_H
        dd[[rr]] <- data.frame(k = k, l = l, axis = j, active_roots = nr,
          importance_weight_sum = wm, H_trace = sum(diag(pair_H)))
        rr <- rr + 1L
      }
    }
    if (!all(is.finite(H)) || max_scaled_residual > 1e-9) stop("Boundary integration failed a numerical integrity check.")
    list(H = sym(H), diagnostics = do.call(rbind, dd), max_root_residual = max_residual,
         max_scaled_root_residual = max_scaled_residual)
  }

  aggregate_H <- function(hs) {
    avg <- Reduce(`+`, hs) / length(hs)
    se <- sqrt(sum(vapply(hs, function(x) sum((x - avg)^2), numeric(1))) /
      (length(hs) * (length(hs) - 1L)))
    list(H = sym(avg), absolute_MCSE = se, relative_MCSE = se / fnorm(avg), estimates = hs)
  }

  install_core <- function(e) {
    original_cov <- e$regularize_cov
    monitor <- new.env(parent = emptyenv()); monitor$active_calls <- 0L; monitor$clipped_eigenvalues <- 0L
    e$regularize_cov <- function(S, rel_floor = e$COV_FLOOR_REL) {
      vals <- eigen(sym(S), symmetric = TRUE, only.values = TRUE)$values
      floor_value <- rel_floor * max(mean(diag(S)), 1e-10)
      count <- sum(vals < floor_value)
      if (count > 0L) monitor$active_calls <- monitor$active_calls + 1L
      monitor$clipped_eigenvalues <- monitor$clipped_eigenvalues + count
      original_cov(S, rel_floor)
    }
    e$fit_ss_qda <- repaired_EM(e$fit_ss_qda)
    # Compiled row maxima replace per-row R callbacks; arithmetic and target
    # likelihoods are unchanged. Independent fixtures verify the score arithmetic.
    e$log_sum_exp_rows <- function(M) {
      mx <- row_max(M)
      mx + log(rowSums(exp(M - mx)))
    }
    e$safe_inverse <- function(M, rel_floor = e$INFO_EIG_FLOOR_REL) inverse(M)
    e$.__original_score_components <- e$compute_score_components
    e$compute_score_components <- function(X, fit, index) score_object(X, make_model(fit))
    e$.__cov_monitor <- monitor
    e
  }

  observed_stationarity <- function(X, y, labelled, classes, fit) {
    m <- make_model(fit); scores <- score_object(X, m)
    n <- nrow(X); p <- m$p; index <- m$index
    resp <- scores$tau; lab <- which(labelled)
    resp[lab, ] <- 0
    resp[cbind(lab, match(as.character(y[lab]), classes))] <- 1
    unl <- as.numeric(!labelled)
    missing <- cores$Landsat$weighted_E_tt(scores, unl) -
      crossprod(scores$bar * sqrt(unl)) / n
    complete <- matrix(0, p, p); normalizer <- matrix(0, p, p)
    alpha_block <- n * (diag(m$pi[seq_len(m$g - 1L)], nrow = m$g - 1L) -
      tcrossprod(m$pi[seq_len(m$g - 1L)]))
    complete[index$alpha, index$alpha] <- alpha_block
    normalizer[index$alpha, index$alpha] <- alpha_block
    gradient <- numeric(p)
    gradient[index$alpha] <- colSums(resp)[seq_len(m$g - 1L)] -
      n * m$pi[seq_len(m$g - 1L)]
    bases <- lapply(seq_len(nrow(index$vech_pairs)), function(a) {
      E <- matrix(0, m$d, m$d); i <- index$vech_pairs[a, 1L]; j <- index$vech_pairs[a, 2L]
      E[i, j] <- 1; E[j, i] <- 1; E
    })
    trace_matrix <- function(M) sum(diag(M))
    for (k in seq_len(m$g)) {
      P <- m$P[[k]]; w <- resp[, k]; nk <- sum(w)
      u <- scores$class_scores[[k]][, seq_len(m$d), drop = FALSE]
      U1 <- as.vector(crossprod(u, w)); U2 <- crossprod(u * sqrt(w))
      mi <- index$mu[[k]]; ci <- index$cov[[k]]
      complete[mi, mi] <- nk * P; normalizer[mi, mi] <- nk * P
      cross <- vapply(bases, function(E) as.vector(P %*% E %*% U1), numeric(m$d))
      complete[mi, ci] <- cross; complete[ci, mi] <- t(cross)
      for (a in seq_along(bases)) for (b in seq_along(bases)) {
        value <- trace_matrix(P %*% bases[[b]] %*% P %*% bases[[a]])
        complete[ci[a], ci[b]] <- trace_matrix(bases[[a]] %*% P %*% bases[[b]] %*% U2) - 0.5 * nk * value
        normalizer[ci[a], ci[b]] <- 0.5 * nk * value
      }
      gradient[index$class[[k]]] <- as.vector(crossprod(scores$class_scores[[k]], w))
    }
    gradient <- gradient / n
    observed <- sym(complete / n - missing); normalizer <- sym(normalizer / n)
    C <- backsolve(chol(normalizer), diag(p))
    minimum <- min(eigen(sym(t(C) %*% observed %*% C), symmetric = TRUE, only.values = TRUE)$values)
    theta <- c(log(m$pi[seq_len(m$g - 1L)] / m$pi[m$g]),
      unlist(lapply(seq_len(m$g), function(k) c(m$mu[k, ],
        m$S[[k]][cbind(index$vech_pairs[, 1L], index$vech_pairs[, 2L])]))))
    maximum_step <- if (minimum > 1e-8) {
      step <- as.vector(inverse(observed) %*% gradient)
      max(abs(step) / (1 + abs(theta)))
    } else Inf
    list(minimum_scaled_observed_Hessian_eigenvalue = minimum,
      maximum_scaled_Newton_step = maximum_step, maximum_absolute_score_per_observation = max(abs(gradient)),
      observed_information = observed, score_per_observation = gradient)
  }

  certify_design <- function(e, score, labelled, candidate, target, B, solution) {
    M <- e$build_information_from_design(score, labelled, candidate, solution$a)
    inv <- inverse(M); objective <- sum(target * t(inv))
    psi <- e$classification_value_from_G(score, sym(inv %*% target %*% inv), candidate)
    vertex <- numeric(length(candidate)); vertex[order(psi, decreasing = TRUE)[seq_len(B)]] <- 1
    gap <- sum((vertex - solution$a) * psi) / nrow(score$tau)
    noise <- 64 * .Machine$double.eps * max(1, abs(objective))
    if (!isTRUE(solution$converged) || gap < -noise || gap / objective > 1.001e-5 ||
        abs(sum(solution$a) - B) > 1e-7 || min(solution$a) < -1e-10 || max(solution$a) > 1 + 1e-10)
      stop("A design failed its independent optimality or feasibility check.")
    selected <- candidate[order(solution$a, decreasing = TRUE)[seq_len(B)]]
    rounded <- numeric(length(candidate)); rounded[order(solution$a, decreasing = TRUE)[seq_len(B)]] <- 1
    Mr <- e$build_information_from_design(score, labelled, candidate, rounded)
    rounded_objective <- sum(target * t(inverse(Mr)))
    if (rounded_objective < objective - max(0, gap) - 4 * noise) stop("Rounding contradicts the certified lower bound.")
    list(a = solution$a, selected = selected, M = M, inverse = inv, objective = objective,
      rounded_objective = rounded_objective, gap = max(0, gap), relative_gap = max(0, gap) / objective,
      psi = psi, iterations = solution$iterations, minimum_information_eigenvalue = min(eigen(M, symmetric = TRUE, only.values = TRUE)$values))
  }

  precision_at_design <- function(agg, design) {
    C <- backsolve(chol(design$M), diag(nrow(design$M)))
    whitened <- lapply(agg$estimates, function(H) sym(t(C) %*% H %*% C))
    summary <- aggregate_H(whitened)
    objective_replicates <- vapply(agg$estimates, function(H) sum(H * t(design$inverse)), numeric(1))
    list(relative_white_H_MCSE = summary$relative_MCSE,
      relative_objective_MCSE = sd(objective_replicates) / sqrt(length(objective_replicates)) / mean(objective_replicates),
      white_norm_sq = vapply(whitened, function(H) sum(H^2), numeric(1)), objective_replicates = objective_replicates)
  }

  constants <- list(EM_MAXIT = protocol$EM_MAXIT, EM_TOL = protocol$EM_TOL,
    COV_FLOOR_REL = protocol$COV_FLOOR_REL, INFO_EIG_FLOOR_REL = 1e-8,
    FW_MAXIT = max(protocol$FW_MAXITS), FW_TOL = protocol$FW_TOL)
  cores <- list(Landsat = definitions(paths$core, constants))
  cores <- lapply(cores, install_core)
  e <- cores$Landsat
  simulator <- definitions(paths$sim1, list(g = 3L, p = 2L, qdim = 17L))
  CLASSES <- as.character(1:3)
  DATASETS <- c("Scenario1", "Scenario2")

  # Separate scenario offsets and replication offsets make all stream seeds
  # unique through replication 100. The largest seed remains below 2^31 - 1.
  stream_seed <- function(ds, rep_id, event) {
    seed <- as.double(protocol$seed_base) + match(ds, DATASETS) * 10000000 + rep_id * 10000 + event
    if (!is.finite(seed) || seed > .Machine$integer.max) stop("Seed outside integer range.")
    as.integer(seed)
  }
  seed_grid <- expand.grid(Dataset = DATASETS, Replication = 0:100,
    Event = c(11L, 12L, 21L, 101:103, 201:218, 301:308, 501L, 511L))
  frozen_seeds <- mapply(stream_seed, seed_grid$Dataset, seed_grid$Replication, seed_grid$Event)
  if (anyDuplicated(frozen_seeds)) stop("The frozen seed streams overlap.")

  extract_truth <- function(path, scenario) {
    local <- new.env(parent = baseenv())
    if (scenario == "Scenario1") {
      allowed <- c("g", "p", "pi0", "mu0", "Sigma0", "true_par")
      covariance <- "Sigma0"; result <- "true_par"
    } else {
      allowed <- c("g", "p", "pi0_s2", "mu0_s2", "Sigma0_s2", "true_par_s2")
      covariance <- "Sigma0_s2"; result <- "true_par_s2"
    }
    for (expr in parse(file = path, keep.source = FALSE)) {
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
      lhs <- expr[[2L]]
      permitted <- is.symbol(lhs) && as.character(lhs) %in% allowed
      permitted <- permitted || (is.call(lhs) && identical(lhs[[1L]], as.name("[")) &&
        is.symbol(lhs[[2L]]) && identical(as.character(lhs[[2L]]), covariance))
      if (permitted) eval(expr, local)
    }
    get(result, envir = local, inherits = FALSE)
  }
  truth <- list(Scenario1 = extract_truth(paths$sim1, "Scenario1"),
    Scenario2 = extract_truth(paths$sim2, "Scenario2"))
  expected_S <- list(
    list(matrix(c(1,.45,.45,.7),2,2), matrix(c(.8,-.35,-.35,1.1),2,2), matrix(c(1.25,.15,.15,.55),2,2)),
    list(matrix(c(1,.15,.15,.9),2,2), matrix(c(.95,.12,.12,.95),2,2), matrix(c(1.05,.18,.18,.85),2,2)))
  for (i in seq_along(DATASETS)) {
    par <- truth[[i]]
    if (!isTRUE(all.equal(as.numeric(par$pi), c(.36,.34,.30), tolerance = 1e-14)) ||
        max(abs(par$mu - rbind(c(-1.6,0),c(1.5,.2),c(0,1.9)))) > 1e-14 ||
        max(vapply(1:3, function(k) max(abs(par$Sigma[,,k] - expected_S[[i]][[k]])), numeric(1))) > 1e-14)
      stop("The extracted simulation truth does not match the reviewed model.")
  }

  collect <- function(name) {
    xs <- lapply(records, function(x) x[[name]])
    if (length(xs)) do.call(rbind, xs) else data.frame()
  }
  checkpoint <- function() {
    for (name in c("raw", "class_errors", "designs", "curvature", "checks"))
      write.csv(collect(name), file.path(out, paste0("simulation_09_", name, ".csv")), row.names = FALSE)
    if (length(current_checks)) write.csv(do.call(rbind, current_checks),
      file.path(out, "simulation_09_current_checks.csv"), row.names = FALSE)
    if (length(failures)) write.csv(do.call(rbind, failures),
      file.path(out, "simulation_09_failures.csv"), row.names = FALSE)
    if (length(scenario_records)) {
      reference <- do.call(rbind, lapply(scenario_records, function(x) x$reference))
      write.csv(reference, file.path(out, "simulation_09_bayes_reference.csv"), row.names = FALSE)
    }
  }
  on.exit(checkpoint(), add = TRUE)

  # Exact conditional normal interval probabilities reduce 2D risk to 1D
  # adaptive quadrature. All pairwise roots are included; intervals are
  # classified using the maximum of all three joint densities.
  risk_quadrature <- function(classifier, population, window = 8) {
    stopifnot(classifier$d == 2L, population$d == 2L)
    lower <- min(population$mu[,1] - window * sqrt(vapply(population$S, function(S) S[1,1], numeric(1))))
    upper <- max(population$mu[,1] + window * sqrt(vapply(population$S, function(S) S[1,1], numeric(1))))
    scalar_integrand <- function(x) {
      roots <- numeric(0)
      for (k in 1:(classifier$g - 1L)) for (l in (k + 1L):classifier$g) {
        A <- (classifier$P[[l]] - classifier$P[[k]]) / 2
        b <- classifier$b[k,] - classifier$b[l,]; cc <- classifier$cc[k] - classifier$cc[l]
        aa <- A[2,2]; bb <- 2*x*A[1,2] + b[2]; cv <- A[1,1]*x*x + b[1]*x + cc
        if (aa == 0) {
          if (bb != 0) roots <- c(roots, -cv/bb)
        } else {
          disc <- bb*bb - 4*aa*cv
          if (disc > 0) {
            q <- -0.5*(bb + if (bb >= 0) sqrt(disc) else -sqrt(disc))
            roots <- c(roots, q/aa, cv/q)
          }
        }
      }
      roots <- sort(unique(roots[is.finite(roots)]))
      cuts <- c(-Inf, roots, Inf)
      a <- head(cuts, -1L); b <- tail(cuts, -1L)
      probe <- numeric(length(a))
      for (i in seq_along(a)) {
        probe[i] <- if (is.finite(a[i]) && is.finite(b[i])) a[i]/2 + b[i]/2 else
          if (is.finite(a[i])) a[i] + max(1, abs(a[i])*.1) else
          if (is.finite(b[i])) b[i] - max(1, abs(b[i])*.1) else 0
      }
      winner <- max.col(log_joint(cbind(rep(x,length(probe)),probe), classifier), ties.method = "first")
      error_density <- 0
      for (k in seq_len(population$g)) {
        S <- population$S[[k]]
        mean_cond <- population$mu[k,2] + S[2,1]/S[1,1]*(x-population$mu[k,1])
        sd_cond <- sqrt(S[2,2] - S[2,1]^2/S[1,1])
        above <- a > mean_cond
        masses <- pnorm(b,mean_cond,sd_cond) - pnorm(a,mean_cond,sd_cond)
        masses[above] <- pnorm(a[above],mean_cond,sd_cond,lower.tail=FALSE) -
          pnorm(b[above],mean_cond,sd_cond,lower.tail=FALSE)
        error_density <- error_density + population$pi[k]*dnorm(x,population$mu[k,1],sqrt(S[1,1]))*
          sum(masses[winner != k])
      }
      error_density
    }
    result <- integrate(function(x) vapply(x,scalar_integrand,numeric(1)), lower, upper,
      rel.tol = 1e-9, abs.tol = 1e-11, subdivisions = 2000L, stop.on.error = TRUE)
    list(risk = result$value, quadrature_error = result$abs.error, message = result$message,
      feature_tail_bound = 2*pnorm(-window), window = window, limits = c(lower,upper))
  }
  perturb_model <- function(par, coordinate, delta) {
    m <- make_model(par); index <- m$index
    theta <- c(log(m$pi[1:(m$g-1L)]/m$pi[m$g]),
      unlist(lapply(seq_len(m$g), function(k) c(m$mu[k,],
        m$S[[k]][cbind(index$vech_pairs[,1],index$vech_pairs[,2])]))))
    theta[coordinate] <- theta[coordinate] + delta
    logits <- c(theta[index$alpha],0); weights <- exp(logits-max(logits))
    updated <- par; updated$pi <- weights/sum(weights)
    for (k in seq_len(m$g)) {
      updated$mu[k,] <- theta[index$mu[[k]]]
      S <- matrix(0,m$d,m$d)
      S[cbind(index$vech_pairs[,1],index$vech_pairs[,2])] <- theta[index$cov[[k]]]
      S <- S + t(S) - diag(diag(S))
      updated$Sigma[,,k] <- S
    }
    make_model(updated)
  }

  # These fixtures verify the parameter-order conversion and the six-method
  # engine against the original simulation functions and the retained core.
  binary <- make_model(list(pi=c(.5,.5),mu=rbind(c(-1,0),c(1,0)),
    Sigma=list(diag(2),diag(2))))
  binary_risk <- risk_quadrature(binary,binary,8)
  binary_H <- atlas(binary,4096L,stream_seed("Scenario1",0L,511L))
  binary_risk_error <- abs(binary_risk$risk-pnorm(-1))
  binary_curvature_error <- abs(binary_H$H[1,1]-dnorm(1)/4)
  if (binary_risk_error>1e-9 || binary_curvature_error>1e-10)
    stop("Analytic equal-covariance Gaussian risk/curvature normalization fixture failed.")
  preflight <- list(binary=list(risk=binary_risk,prior_curvature=binary_H$H[1,1],
    analytic_risk=pnorm(-1),analytic_prior_curvature=dnorm(1)/4,
    risk_discrepancy=binary_risk_error,curvature_discrepancy=binary_curvature_error))
  cat("Analytic binary Gaussian risk and curvature fixtures: passed.\n")
  fixtures <- rbind(c(0,0),c(-1.1,.7),c(1.9,-.4),c(.3,2.2))
  for (ds in DATASETS) {
    m <- make_model(truth[[ds]]); score <- score_object(fixtures,m)
    permutation <- c(m$index$alpha, unlist(m$index$mu), unlist(m$index$cov))
    score_error <- 0; density_error <- max(abs(log_joint(fixtures,m) - simulator$log_weighted_density(fixtures,truth[[ds]])))
    for (i in seq_len(nrow(fixtures))) for (k in 1:3) {
      full <- numeric(m$p); full[m$index$alpha] <- score$A[k,]
      full[m$index$class[[k]]] <- score$class_scores[[k]][i,]
      score_error <- max(score_error,max(abs(full[permutation] - simulator$class_score(fixtures[i,],k,truth[[ds]]))))
    }
    native <- e$.__original_score_components(fixtures,
      list(pi=m$pi,mu=m$mu,Sigma=m$S,classes=CLASSES),m$index)
    native_error <- max(abs(score$bar-native$bar), abs(score$tau-native$tau),
      abs(e$weighted_E_tt(score,rep(1,nrow(fixtures))) - e$weighted_E_tt(native,rep(1,nrow(fixtures)))))
    if (score_error > 1e-10 || density_error > 1e-10 || native_error > 1e-10)
      stop("Simulation score/density/information fixture failed.")
    preflight[[ds]] <- list(original_score_max_error=score_error,
      original_density_max_error=density_error, original_core_max_error=native_error)
    cat(ds," original-score/density/information fixtures: passed.\n")
  }
  saveRDS(preflight,file.path(out,"simulation_09_preflight.rds"))

  fit_once <- function(X, z, labelled, ds, rep_id, method, budget = NA_real_, init_event = 21L) {
    cat("  ",ds," replication ",rep_id,"; fitting ",method,
      if (is.finite(budget)) paste0("; budget=",budget) else "","\n",sep="")
    flush.console()
    monitor <- e$.__cov_monitor; monitor$active_calls <- 0L; monitor$clipped_eigenvalues <- 0L
    labels <- rep(NA_integer_,nrow(X)); labels[labelled] <- z[labelled]
    counts <- tabulate(labels[labelled],nbins=3L)
    if (any(counts <= ncol(X)+1L)) stop("A labeled class is too small for this preflight; do not resample it. Counts: ",paste(counts,collapse=","))
    set.seed(stream_seed(ds,rep_id,init_event))
    init <- simulator$initial_from_labels(X,labels,3L)
    warm <- list(pi=init$pi,mu=init$mu,Sigma=lapply(1:3,function(k) init$Sigma[,,k]),
      classes=CLASSES,logLik_history=numeric(0),iterations=0L)
    used <- 0L
    for (level in seq_along(protocol$EM_PARAMETER_TOLS)) {
      remaining <- protocol$EM_MAXIT - used
      if (remaining < protocol$EM_STABLE_ITERATIONS+1L) stop("EM exhausted the fixed iteration limit.")
      ptol <- protocol$EM_PARAMETER_TOLS[level]
      fit <- e$fit_ss_qda(X,z,labelled,CLASSES,maxit=remaining,tol=protocol$EM_TOL,
        start_fit=warm,parameter_tol=ptol,stable_required=protocol$EM_STABLE_ITERATIONS)
      used <- used + fit$additional_iterations
      lj <- log_joint(X,make_model(fit)); lab <- which(labelled); unl <- which(!labelled)
      fresh <- sum(lj[cbind(lab,z[lab])]) + if (length(unl)) sum(e$log_sum_exp_rows(lj[unl,,drop=FALSE])) else 0
      mismatch <- abs(fresh-fit$logLik); increment <- min(diff(fit$logLik_history))
      if (!isTRUE(fit$converged) || fit$stable_joint_iterations < protocol$EM_STABLE_ITERATIONS ||
          fit$relative_logLik_change >= protocol$EM_TOL || fit$max_scaled_parameter_change >= ptol ||
          !is.finite(fresh) || mismatch > 1e-9 + 64*.Machine$double.eps*(1+abs(fresh)) ||
          increment < -1e-10*(1+max(abs(fit$logLik_history)))) stop("Joint EM integrity check failed.")
      for (S in fit$Sigma) if (min(eigen(sym(S),symmetric=TRUE,only.values=TRUE)$values) <=
          1.000001*protocol$COV_FLOOR_REL*max(mean(diag(S)),1e-10))
        stop("A final covariance is on the floor; unconstrained stationarity needs review.")
      diagnostic <- observed_stationarity(X,z,labelled,CLASSES,fit)
      add(ds,rep_id,"fit","converged",fit$converged,method,budget,paste0("level=",level))
      add(ds,rep_id,"fit","fresh_logLik_discrepancy",mismatch,method,budget,paste0("level=",level))
      add(ds,rep_id,"fit","minimum_logLik_increment",increment,method,budget,paste0("level=",level))
      add(ds,rep_id,"fit","relative_logLik_change",fit$relative_logLik_change,method,budget,paste0("level=",level))
      add(ds,rep_id,"fit","maximum_scaled_parameter_change",fit$max_scaled_parameter_change,method,budget,paste0("level=",level))
      add(ds,rep_id,"stationarity","minimum_scaled_observed_Hessian_eigenvalue",
        diagnostic$minimum_scaled_observed_Hessian_eigenvalue,method,budget,paste0("level=",level))
      add(ds,rep_id,"stationarity","maximum_scaled_Newton_step",diagnostic$maximum_scaled_Newton_step,
        method,budget,paste0("level=",level))
      add(ds,rep_id,"stationarity","maximum_absolute_score_per_observation",
        diagnostic$maximum_absolute_score_per_observation,method,budget,paste0("level=",level))
      if (diagnostic$minimum_scaled_observed_Hessian_eigenvalue <= 1e-8)
        stop("Observed likelihood curvature is not reliably positive definite; this fit needs review.")
      if (diagnostic$maximum_scaled_Newton_step <= protocol$EM_NEWTON_TOL) break
      warm <- fit
    }
    if (diagnostic$maximum_scaled_Newton_step > protocol$EM_NEWTON_TOL)
      stop("Observed-likelihood stationarity failed after all fixed parameter tolerances.")
    fit$stationarity <- diagnostic[c("minimum_scaled_observed_Hessian_eigenvalue",
      "maximum_scaled_Newton_step","maximum_absolute_score_per_observation")]
    fit$initialization <- init; fit$initialization_seed <- stream_seed(ds,rep_id,init_event)
    fit$covariance_floor_active_calls <- monitor$active_calls
    fit$covariance_clipped_eigenvalues <- monitor$clipped_eigenvalues
    add(ds,rep_id,"fit","EM_iterations",fit$iterations,method,budget)
    add(ds,rep_id,"fit","covariance_floor_active_calls",monitor$active_calls,method,budget)
    add(ds,rep_id,"fit","covariance_clipped_eigenvalues",monitor$clipped_eigenvalues,method,budget)
    fit
  }
  solve_one <- function(score,labelled,candidate,target,B,ds,rep_id,method,budget) {
    for (limit in protocol$FW_MAXITS) {
      solution <- e$solve_finite_pool_design(score,labelled,candidate,target,B,maxit=limit,tol=protocol$FW_TOL)
      add(ds,rep_id,"design","solver_relative_FW_gap",solution$relative_fw_gap,method,budget,paste0("limit=",limit))
      if (isTRUE(solution$converged)) break
    }
    design <- certify_design(e,score,labelled,candidate,target,B,solution)
    add(ds,rep_id,"design","independent_relative_FW_gap",design$relative_gap,method,budget)
    design
  }
  integrate_H <- function(model,n,seeds,ds,label) {
    estimates <- vector("list",length(seeds)); diagnostics <- vector("list",length(seeds))
    for (i in seq_along(seeds)) {
      cat("  ",ds," ",label,"; integration ",i,"/",length(seeds),"; n/chart=",n,"\n",sep="")
      z <- atlas(model,n,seeds[i])
      estimates[[i]] <- z$H
      diagnostics[[i]] <- list(seed=seeds[i],table=z$diagnostics,
        maximum_root_residual=z$max_root_residual,maximum_scaled_root_residual=z$max_scaled_root_residual)
    }
    result <- aggregate_H(estimates); result$diagnostics <- diagnostics
    result
  }

  build_scenario <- function(ds) {
    cat("\nSTART scenario reference: ",ds,"\n",sep="")
    model <- make_model(truth[[ds]])
    mc_seeds <- vapply(301:308,function(event) stream_seed(ds,0L,event),integer(1))
    agg <- integrate_H(model,protocol$true_mc_n,mc_seeds,ds,"truth")
    if (!is.finite(agg$relative_MCSE) || agg$relative_MCSE > protocol$max_H_MCSE)
      stop("True-model curvature precision failed at the frozen integration effort.")
    q8 <- risk_quadrature(model,model,8); q10 <- risk_quadrature(model,model,10)
    if (abs(q8$risk-q10$risk) > max(1e-8,5*(q8$quadrature_error+q10$quadrature_error)))
      stop("Bayes risk quadrature did not stabilize with the feature window.")
    # Independent directional risk differences also check normalization.
    directions <- c(Prior=model$index$alpha[1],Mean=model$index$mu[[1]][2],
      Covariance=model$index$cov[[1]][2])
    curvature_rows <- list(); cr <- 1L
    for (direction in names(directions)) {
      coordinate <- directions[[direction]]
      values <- vapply(agg$estimates,function(H) H[coordinate,coordinate],numeric(1))
      prediction <- mean(values); mcse <- sd(values)/sqrt(length(values))
      for (step in c(.01,.005)) {
        plus <- risk_quadrature(perturb_model(truth[[ds]],coordinate,step),model,8)
        minus <- risk_quadrature(perturb_model(truth[[ds]],coordinate,-step),model,8)
        fd <- (plus$risk+minus$risk-2*q8$risk)/(step*step)
        numerical <- (plus$quadrature_error+minus$quadrature_error+2*q8$quadrature_error)/(step*step)
        tolerance <- 6*mcse + .01*abs(prediction) + 10*numerical
        passed <- abs(fd-prediction) <= tolerance
        curvature_rows[[cr]] <- data.frame(Dataset=ds,Direction=direction,Coordinate=coordinate,
          Step=step,BoundaryCurvature=prediction,BoundaryMCSE=mcse,DirectRiskCurvature=fd,
          PropagatedQuadratureDiagnostic=numerical,CompatibilityTolerance=tolerance,Passed=passed)
        cr <- cr+1L
      }
    }
    direct <- do.call(rbind,curvature_rows)
    if (!all(direct$Passed)) {
      write.csv(direct,file.path(out,paste0(ds,"_direct_risk_FAILED.csv")),row.names=FALSE)
      stop("An independent direct-risk curvature check failed; review before any replication.")
    }
    # The fixed evaluation sample is generated only after precision checks.
    test_seed <- stream_seed(ds,0L,501L); set.seed(test_seed)
    test <- simulator$simulate_qda(protocol$test_n,truth[[ds]])
    pred <- max.col(log_joint(test$Y,model),ties.method="first")
    count <- sum(pred != test$z); error <- count/protocol$test_n
    reference <- data.frame(Dataset=ds,NTest=protocol$test_n,TestSeed=test_seed,
      BayesErrorCount=count,BayesMonteCarloError=error,
      BayesMonteCarloMCSE=sqrt(error*(1-error)/protocol$test_n),
      BayesQuadratureRisk=q8$risk,BayesQuadratureDiagnostic=q8$quadrature_error,
      BayesQuadratureWindowChange=abs(q8$risk-q10$risk),
      TrueCurvatureRelativeMCSE=agg$relative_MCSE)
    result <- list(Dataset=ds,protocol_version=protocol$version,truth=truth[[ds]],
      H=agg$H,mc_estimates=agg$estimates,mc_seeds=mc_seeds,mc_n=protocol$true_mc_n,
      mc_relative_se=agg$relative_MCSE,mc_diagnostics=agg$diagnostics,
      test=test,test_seed=test_seed,bayes_prediction_md5=hash_raw(as.raw(pred)),
      quadrature=list(window8=q8,window10=q10),direct_risk=direct,reference=reference,
      reference_checks=if (length(current_checks)) do.call(rbind,current_checks) else data.frame())
    atomic_save(result,file.path(out,paste0(ds,"_reference.rds")))
    cat("END scenario reference: ",ds,"; Bayes quadrature=",format(q8$risk,digits=12),"\n",sep="")
    result
  }

  run_replication <- function(ds,rep_id) {
    cat("\nSTART ",ds," replication ",rep_id,"/",n_rep,"\n",sep="")
    current_checks <<- list()
    reference <- scenario_records[[ds]]
    training_seed <- stream_seed(ds,rep_id,11L); set.seed(training_seed)
    training <- simulator$simulate_qda(protocol$n,truth[[ds]])
    X <- training$Y; z <- training$z
    pilot_seed <- stream_seed(ds,rep_id,12L); set.seed(pilot_seed)
    pilot <- sort(sample.int(protocol$n,protocol$pilot_n,replace=FALSE))
    labelled <- seq_len(protocol$n) %in% pilot; candidate <- which(!labelled)
    pilot_counts <- tabulate(z[pilot],nbins=3L)
    add(ds,rep_id,"pilot","class_counts",paste(pilot_counts,collapse=","))
    cat("  Pilot class counts: ",paste(pilot_counts,collapse=" / "),"\n",sep="")
    fit <- fit_once(X,z,labelled,ds,rep_id,"Pilot")
    model <- make_model(fit); score <- score_object(X,model)
    truth_score <- score_object(X,make_model(truth[[ds]]))
    total_n <- as.integer(round(protocol$n*protocol$budgets))
    B <- total_n-protocol$pilot_n
    fisher_target <- e$weighted_E_tt(score,rep(1,protocol$n))
    fisher <- lapply(seq_along(B),function(bi) solve_one(score,labelled,candidate,
      fisher_target,B[bi],ds,rep_id,"Fisher",protocol$budgets[bi]))
    oracle <- lapply(seq_along(B),function(bi) solve_one(truth_score,labelled,candidate,
      reference$H,B[bi],ds,rep_id,"OracleRisk",protocol$budgets[bi]))
    truth_agg <- list(estimates=reference$mc_estimates)
    oracle_precision <- lapply(oracle,function(design) precision_at_design(truth_agg,design))
    if (!all(vapply(oracle_precision,function(p) is.finite(p$relative_white_H_MCSE) &&
        is.finite(p$relative_objective_MCSE) && p$relative_white_H_MCSE <= protocol$max_white_H_MCSE &&
        p$relative_objective_MCSE <= protocol$max_objective_MCSE,logical(1))))
      stop("Oracle-design curvature precision failed; the fixed truth integration needs review.")
    mc_seeds <- vapply(301:308,function(event) stream_seed(ds,rep_id,event),integer(1))
    accepted_level <- NA_integer_; old_H <- NULL
    for (level in seq_along(protocol$pilot_mc_sizes)) {
      agg <- integrate_H(model,protocol$pilot_mc_sizes[level],mc_seeds,ds,paste0("pilot rep ",rep_id))
      add(ds,rep_id,"boundary_precision","relative_H_MCSE",agg$relative_MCSE,detail=paste0("level=",level))
      if (!is.null(old_H)) add(ds,rep_id,"boundary_precision","relative_H_change_from_previous_effort",
        relative(agg$H,old_H),detail=paste0("level=",level))
      old_H <- agg$H
      if (!is.finite(agg$relative_MCSE) || agg$relative_MCSE > protocol$max_H_MCSE) next
      adaptive <- lapply(seq_along(B),function(bi) solve_one(score,labelled,candidate,
        agg$H,B[bi],ds,rep_id,"AdaptiveRisk",protocol$budgets[bi]))
      precision <- lapply(adaptive,function(design) precision_at_design(agg,design))
      for (bi in seq_along(B)) {
        add(ds,rep_id,"boundary_precision","relative_information_weighted_H_MCSE",
          precision[[bi]]$relative_white_H_MCSE,"AdaptiveRisk",protocol$budgets[bi],paste0("level=",level))
        add(ds,rep_id,"boundary_precision","relative_objective_MCSE",
          precision[[bi]]$relative_objective_MCSE,"AdaptiveRisk",protocol$budgets[bi],paste0("level=",level))
      }
      passed <- all(vapply(precision,function(p) is.finite(p$relative_white_H_MCSE) &&
        is.finite(p$relative_objective_MCSE) && p$relative_white_H_MCSE <= protocol$max_white_H_MCSE &&
        p$relative_objective_MCSE <= protocol$max_objective_MCSE,logical(1)))
      if (passed) {accepted_level <- level; break}
    }
    if (is.na(accepted_level)) stop("Pilot curvature precision failed at maximum effort; review this replication without skipping it.")
    entropy <- e$score_entropy(score$tau); margin <- e$score_margin(score$tau)
    selections <- vector("list",length(B)); selection_seeds <- integer(length(B))
    for (bi in seq_along(B)) {
      selection_seeds[bi] <- stream_seed(ds,rep_id,100L+bi); set.seed(selection_seeds[bi])
      selections[[bi]] <- list(Random=sample(candidate,B[bi],replace=FALSE),
        Entropy=candidate[order(entropy[candidate],decreasing=TRUE)[seq_len(B[bi])]],
        Margin=candidate[order(margin[candidate],decreasing=TRUE)[seq_len(B[bi])]],
        Fisher=fisher[[bi]]$selected,AdaptiveRisk=adaptive[[bi]]$selected,OracleRisk=oracle[[bi]]$selected)
    }
    # All designs and curvature decisions are now frozen before test use.
    final <- list(); raws <- list(); class_rows <- list(); design_rows <- list(); curvature_rows <- list(); rr <- 1L
    for (bi in seq_along(B)) {
      budget <- protocol$budgets[bi]
      for (mi in seq_along(METHODS)) {
        method <- METHODS[mi]; selected <- selections[[bi]][[method]]
        if (length(selected) != B[bi] || anyDuplicated(selected) || length(intersect(selected,pilot)))
          stop("Invalid exact acquisition set.")
        lab <- labelled; lab[selected] <- TRUE
        final_fit <- fit_once(X,z,lab,ds,rep_id,method,budget,init_event=200L+6L*(bi-1L)+mi)
        prediction <- max.col(log_joint(reference$test$Y,make_model(final_fit)),ties.method="first")
        error_count <- sum(prediction != reference$test$z)
        class_n <- tabulate(reference$test$z,nbins=3L)
        class_err <- vapply(1:3,function(k) sum(prediction[reference$test$z==k] != k),integer(1))
        error <- error_count/protocol$test_n; balanced <- mean(class_err/class_n)
        key <- paste0(method,"_",budget)
        final[[key]] <- list(selected=selected,fit=final_fit,error_count=error_count,
          class_error_count=class_err,class_n=class_n,error=error,balanced_error=balanced,
          prediction_md5=hash_raw(as.raw(prediction)))
        raws[[rr]] <- data.frame(Dataset=ds,Replication=rep_id,Budget=budget,Method=method,
          NLabelled=sum(lab),NTest=protocol$test_n,ErrorCount=error_count,Error=error,
          BalancedError=balanced,LogLik=final_fit$logLik,EMIterations=final_fit$iterations,
          EMConverged=final_fit$converged,
          MaximumScaledNewtonStep=final_fit$stationarity$maximum_scaled_Newton_step)
        class_rows[[rr]] <- data.frame(Dataset=ds,Replication=rep_id,Budget=budget,Method=method,
          Class=CLASSES,NTest=class_n,ErrorCount=class_err,ClassError=class_err/class_n)
        cat("  ",method,"; budget=",budget,"; error=",sprintf("%.6f",error),
          "; labels=",sum(lab),"; EM=",final_fit$iterations,"\n",sep="")
        rr <- rr+1L
      }
      for (method in c("Fisher","AdaptiveRisk","OracleRisk")) {
        design <- switch(method,Fisher=fisher[[bi]],AdaptiveRisk=adaptive[[bi]],OracleRisk=oracle[[bi]])
        design_rows[[length(design_rows)+1L]] <- data.frame(Dataset=ds,Replication=rep_id,Budget=budget,Method=method,
          AdditionalBudget=B[bi],Iterations=design$iterations,RelativeFWGap=design$relative_gap,
          RelaxedObjective=design$objective,RoundedObjective=design$rounded_objective,
          RoundingLossPct=100*(design$rounded_objective-design$objective)/design$objective,
          FractionalCount=sum(design$a>1e-8 & design$a<1-1e-8),
          MinimumInformationEigenvalue=design$minimum_information_eigenvalue)
      }
      for (method in c("AdaptiveRisk","OracleRisk")) {
        p <- if (method=="AdaptiveRisk") precision[[bi]] else oracle_precision[[bi]]
        curvature_rows[[length(curvature_rows)+1L]] <- data.frame(Dataset=ds,Replication=rep_id,Budget=budget,Method=method,
          NPerChart=if (method=="AdaptiveRisk") protocol$pilot_mc_sizes[accepted_level] else protocol$true_mc_n,
          IntegrationReplicates=protocol$mc_reps,
          RelativeHMCSE=if (method=="AdaptiveRisk") agg$relative_MCSE else reference$mc_relative_se,
          RelativeInformationWeightedHMCSE=p$relative_white_H_MCSE,
          RelativeObjectiveMCSE=p$relative_objective_MCSE,NumericalPrecisionPassed=TRUE)
      }
    }
    record <- list(Dataset=ds,Replication=rep_id,protocol_version=protocol$version,
      training=training,training_seed=training_seed,pilot=pilot,pilot_seed=pilot_seed,
      pilot_fit=fit,random_selection_seeds=selection_seeds,H=agg$H,mc_estimates=agg$estimates,
      mc_seeds=mc_seeds,mc_n=protocol$pilot_mc_sizes[accepted_level],mc_relative_se=agg$relative_MCSE,
      mc_diagnostics=agg$diagnostics,precision=precision,oracle_precision=oracle_precision,
      fisher=fisher,fisher_target=fisher_target,adaptive=adaptive,oracle=oracle,
      entropy=entropy,margin=margin,final=final,raw=do.call(rbind,raws),
      class_errors=do.call(rbind,class_rows),designs=do.call(rbind,design_rows),
      curvature=do.call(rbind,curvature_rows),checks=do.call(rbind,current_checks))
    atomic_save(record,file.path(out,"fits",sprintf("%s_replication_%03d.rds",ds,rep_id)))
    cat("END ",ds," replication ",rep_id," -- numerical checks passed\n",sep="")
    record
  }

  guarded <- function(task,ds,rep_id) {
    tryCatch(withCallingHandlers(task(),warning=function(w) {
      add(ds,rep_id,"execution","warning",conditionMessage(w))
      cat("WARNING: ",conditionMessage(w),"\n",sep=""); invokeRestart("muffleWarning")
    }),error=function(err) {
      failures[[length(failures)+1L]] <<- data.frame(Dataset=ds,Replication=rep_id,Reason=conditionMessage(err))
      checkpoint(); cat("\nRUN STOPPED: ",conditionMessage(err),"\nOutput: ",out,"\n",sep="")
      stop(err)
    })
  }
  for (ds in DATASETS) {
    reference_path <- file.path(out,paste0(ds,"_reference.rds"))
    reference <- if (file.exists(reference_path)) readRDS(reference_path) else
      guarded(function() build_scenario(ds),ds,0L)
    if (!identical(reference$protocol_version,protocol$version) || reference$Dataset!=ds)
      stop("Cached scenario reference does not match the protocol.")
    scenario_records[[ds]] <- reference; checkpoint()
    write.csv(reference$direct_risk,file.path(out,paste0(ds,"_direct_risk.csv")),row.names=FALSE)
    for (rep_id in seq_len(n_rep)) {
      fit_path <- file.path(out,"fits",sprintf("%s_replication_%03d.rds",ds,rep_id))
      record <- if (file.exists(fit_path)) readRDS(fit_path) else
        guarded(function() run_replication(ds,rep_id),ds,rep_id)
      if (!identical(record$protocol_version,protocol$version) || record$Dataset!=ds || record$Replication!=rep_id)
        stop("Cached replication does not match the protocol.")
      records[[paste0(ds,"_",rep_id)]] <- record; checkpoint(); gc(verbose=FALSE)
    }
  }
  raw <- collect("raw")
  expected <- length(DATASETS)*n_rep*length(protocol$budgets)*length(METHODS)
  if (nrow(raw)!=expected || anyDuplicated(raw[,c("Dataset","Replication","Budget","Method")]))
    stop("The requested experiment is incomplete or has duplicate records.")
  summaries <- list(); paired_rows <- list(); sr <- 1L; pr <- 1L
  for (ds in DATASETS) for (budget in protocol$budgets) for (outcome in c("Error","BalancedError")) {
    subset <- raw[raw$Dataset==ds & raw$Budget==budget,]
    for (method in METHODS) {
      values <- subset[subset$Method==method,outcome]
      summaries[[sr]] <- data.frame(Dataset=ds,Budget=budget,Outcome=outcome,Method=method,N=length(values),
        Mean=mean(values),MCSE=if (length(values)>1L) sd(values)/sqrt(length(values)) else NA_real_)
      sr <- sr+1L
    }
    adaptive <- subset[subset$Method=="AdaptiveRisk",c("Replication",outcome)]
    for (method in setdiff(METHODS,"AdaptiveRisk")) {
      comparator <- subset[subset$Method==method,c("Replication",outcome)]
      comparator <- comparator[match(adaptive$Replication,comparator$Replication),]
      differences <- adaptive[[outcome]]-comparator[[outcome]]
      nn <- length(differences); se <- if (nn>1L) sd(differences)/sqrt(nn) else NA_real_
      critical <- if (nn>1L) qt(.975,nn-1L) else NA_real_
      p <- if (is.finite(se) && se>0) 2*pt(-abs(mean(differences)/se),nn-1L) else
        if (nn>1L && all(differences==0)) 1 else NA_real_
      paired_rows[[pr]] <- data.frame(Dataset=ds,Budget=budget,Outcome=outcome,
        Contrast=paste0("AdaptiveRisk-minus-",method),N=nn,MeanDifference=mean(differences),MCSE=se,
        Lower95=mean(differences)-critical*se,Upper95=mean(differences)+critical*se,PValue=p)
      pr <- pr+1L
    }
  }
  summary <- do.call(rbind,summaries); paired <- do.call(rbind,paired_rows)
  paired$HolmPValue <- NA_real_
  for (ds in DATASETS) for (outcome in c("Error","BalancedError")) {
    ii <- which(paired$Dataset==ds & paired$Outcome==outcome)
    paired$HolmPValue[ii] <- p.adjust(paired$PValue[ii],method="holm")
  }
  status <- data.frame(Dataset=DATASETS,RequestedReplications=n_rep,
    CompleteReplications=vapply(DATASETS,function(ds) length(unique(raw$Replication[raw$Dataset==ds])),integer(1)),
    RawRows=vapply(DATASETS,function(ds) sum(raw$Dataset==ds),integer(1)),
    ExecutionComplete=TRUE,ManuscriptReplicationTargetMet=n_rep==100L)
  write.csv(summary,file.path(out,"simulation_09_summary.csv"),row.names=FALSE)
  write.csv(paired,file.path(out,"simulation_09_paired.csv"),row.names=FALSE)
  write.csv(status,file.path(out,"simulation_09_status.csv"),row.names=FALSE)
  review <- list(protocol=protocol,status=status,preflight=preflight,scenarios=scenario_records,
    records=records,summary=summary,paired=paired)
  atomic_save(review,file.path(out,"simulation_09_review.rds")); checkpoint()
  cat("\nCompleted requested diagnostic:\n"); print(status)
  cat("\nFor the first two-replication run, ZIP the ENTIRE output folder and return that ZIP.\n")
  cat("Include the review RDS, both scenario references, every fit RDS, all CSV files, and the log.\n")
  cat("For an interrupted run, retain the folder and use its exact printed path as resume_dir.\n")
  cat("Do not use two-replication results to replace manuscript means or select methods/settings.\n")
  cat("All 100 replications per scenario, independent result checks, and regenerated figures are still required.\n")
  cat("Negative paired differences favor AdaptiveRisk; intervals are pointwise.\n")
  cat("Holm adjustment covers 15 comparisons per scenario/outcome, including the oracle benchmark.\n")
  cat("Convergence checks establish local stationarity, not a global maximum or submission readiness.\n")
  invisible(out)
}
