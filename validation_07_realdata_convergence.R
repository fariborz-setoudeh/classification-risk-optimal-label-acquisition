# Joint-convergence real-data experiments for Classification-risk-optimal label acquisition.
# Source this file, then call validation_07_realdata(...).
# Every completed replication is saved. Resume an interrupted run by passing
# its printed output directory as resume_dir. Original files are not edited.

validation_07_realdata <- function(project_root, validation_06_dir,
                                   n_rep = 50L, datasets = c("Landsat", "DryBean"),
                                   resume_dir = NULL) {
  stopifnot(dir.exists(project_root), length(n_rep) == 1L, n_rep >= 1L,
    n_rep <= 50L, n_rep == as.integer(n_rep), length(datasets) >= 1L,
    !anyDuplicated(datasets), all(datasets %in% c("Landsat", "DryBean")))
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  validation_06_dir <- normalizePath(validation_06_dir, winslash = "/", mustWork = TRUE)
  validation_06_rds <- file.path(validation_06_dir, "realdata_final_review.rds")
  if (!file.exists(validation_06_rds)) stop("Missing validation 06 review cache.")
  files <- list.files(project_root, recursive = TRUE, full.names = TRUE)
  find_one <- function(names) {
    hits <- files[basename(files) %in% names]
    if (length(hits) != 1L) stop("Expected one of ", paste(names, collapse = "/"), "; found ", length(hits), ".")
    normalizePath(hits, winslash = "/", mustWork = TRUE)
  }
  paths <- list(landsat = find_one("03_landsat_analysis.R"),
    bean = find_one(c("03_drybean_analysis.R.R", "03_drybean_analysis.R")),
    trn = find_one("sat.trn"), tst = find_one("sat.tst"), xlsx = find_one("Dry_Bean_Dataset.xlsx"))
  expected <- c("67b711404f884efec1a3f1d51340cbf5", "d92057b86d189070ed12cc3caff12621",
    "2c5ba2900da0183cab2c41fdb279fa5b", "02c995991fecc864e809b2c4c42cd983", "a024b18787414f9881d21faae54ce9c7")
  hashes <- unname(tools::md5sum(unlist(paths)))
  if (!identical(hashes, expected)) stop("Use the reviewed, unedited original scripts and datasets.")
  if (unname(tools::md5sum(validation_06_rds)) != "9bf0e016c869f3d82916c759830ceeb5") stop("Use the fully audited validation 06 cache.")
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Your existing readxl installation is required.")
  shadows <- c("pi", "sample", "set.seed", "cov", "sd", "prcomp", "eigen", "solve", "chol", "dnorm", "order", "max.col", "pmax")
  shadows <- shadows[vapply(shadows, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(shadows)) stop("Global objects shadow numerical functions/constants: ", paste(shadows, collapse = ", "))
  previous <- readRDS(validation_06_rds)
  if (length(previous$records) != 100L || any(previous$status$CompleteReplications != 50L)) stop("The validation 06 cache is incomplete.")
  protocol <- list(version = "2026-10-03-realdata-joint-v1", hashes = hashes,
    validation_06_review_md5 = "9bf0e016c869f3d82916c759830ceeb5",
    RNG = c("Mersenne-Twister", "Inversion", "Rejection"), budgets = c(0.10, 0.20, 0.30), pilot_fraction = 0.05,
    EM_MAXIT = 20000L, EM_TOL = 1e-11, EM_PARAMETER_TOLS = c(1e-8, 1e-9, 1e-10, 1e-11),
    EM_STABLE_ITERATIONS = 5L, EM_NEWTON_TOL = 1e-6, COV_FLOOR_REL = 1e-6, FW_TOL = 1e-5,
    seed_master = c(Landsat = 20260902L, DryBean = 20260901L), mc_seed_base = 2026105000L,
    mc_case_id = c(Landsat = 5L, DryBean = 6L), mc_reps = 8L,
    mc_sizes = c(131072L, 262144L, 524288L), max_H_MCSE = 0.01,
    max_white_H_MCSE = 0.01, max_objective_MCSE = 0.002,
    R_version = R.version.string, readxl_version = as.character(utils::packageVersion("readxl")))
  out <- if (is.null(resume_dir)) tempfile("realdata_converged_", tmpdir = project_root) else
    normalizePath(resume_dir, winslash = "/", mustWork = TRUE)
  if (is.null(resume_dir)) {
    if (!dir.create(out)) stop("Cannot create output directory.")
    dir.create(file.path(out, "fits")); saveRDS(protocol, file.path(out, "protocol.rds"))
  } else {
    if (!identical(readRDS(file.path(out, "protocol.rds")), protocol)) stop("Resume requires the same protocol, inputs, and R environment.")
  }
  logfile <- file.path(out, "realdata_final_log.txt")
  sink(logfile, append = !is.null(resume_dir), split = TRUE)
  initial_rng <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  initial_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(initial_rng)))
    if (had_seed) assign(".Random.seed", initial_seed, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    sink()
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  cat("Joint-convergence real-data rerun; ", format(Sys.time()), "\nOutput: ", out, "\n", sep = "")
  cat("Target replications per requested dataset: ", n_rep, "\n", sep = "")
  cat("Primary curvature: fitted-Gaussian boundary integral without a bandwidth.\n")
  cat("Integration precision is selected using pilot fits and training features, before test evaluation.\n")
  cat("All three budgets and all five methods are retained. Numerical failures stop the run.\n")
  cat("EM requires likelihood and parameter stability for five updates, plus observed-likelihood stationarity.\n")
  cat("Cached fits are continued only when their labeled sets are identical; changed sets use the original initialization.\n")
  cat("The original raw results, scripts, datasets, and manuscript are unchanged.\n")
  print(data.frame(path = unlist(paths), md5 = hashes)); print(protocol); print(sessionInfo())

  METHODS <- c("Random", "Entropy", "Margin", "Fisher", "AdaptiveRisk")
  tables <- list(raw = list(), class_errors = list(), designs = list(), curvature = list(), checks = list())
  current_checks <- list(); failures <- list()
  add <- function(ds, rep_id, stage, metric, value, method = "", budget = NA_real_, detail = "") {
    current_checks[[length(current_checks) + 1L]] <<- data.frame(Dataset = ds, Replication = rep_id,
      Stage = stage, Method = method, Budget = budget, Metric = metric, Value = as.character(value), Detail = detail,
      stringsAsFactors = FALSE)
  }
  append_record <- function(record) {
    for (name in names(tables)) tables[[name]][[length(tables[[name]]) + 1L]] <<- record[[name]]
  }
  merge_tables <- function() lapply(tables, function(xs) if (length(xs)) do.call(rbind, xs) else data.frame())
  checkpoint <- function() {
    merged <- merge_tables()
    for (name in names(merged)) write.csv(merged[[name]], file.path(out, paste0("realdata_final_", name, ".csv")), row.names = FALSE)
    if (length(failures)) write.csv(do.call(rbind, failures), file.path(out, "realdata_final_failures.csv"), row.names = FALSE)
    if (length(current_checks)) write.csv(do.call(rbind, current_checks), file.path(out, "current_replication_checks.csv"), row.names = FALSE)
  }
  on.exit(checkpoint(), add = TRUE)
  definitions <- function(path, constants) {
    e <- new.env(parent = globalenv()); list2env(constants, e)
    for (expr in parse(file = path, keep.source = FALSE)) {
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
      rhs <- expr[[3L]]
      if (is.symbol(expr[[2L]]) && is.call(rhs) && identical(rhs[[1L]], as.name("function"))) eval(expr, e)
    }
    e
  }
  constants <- list(EM_MAXIT = 20000L, EM_TOL = 1e-11, COV_FLOOR_REL = 1e-6, INFO_EIG_FLOOR_REL = 1e-8,
    BOUNDARY_BW_MIN = 0.05, BOUNDARY_BW_MAX = 1, BOUNDARY_KERNEL_CUTOFF = 4, MIN_BOUNDARY_POINTS = 10L,
    FW_MAXIT = 1200L, FW_TOL = 1e-5, N_PC = 5L)
  cores <- list(Landsat = definitions(paths$landsat, constants), DryBean = definitions(paths$bean, constants))
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

  # For ell_kl=0, partition each regular face by w_j=lambda_j*ell_j^2 /
  # sum_i(lambda_i*ell_i^2). Coarea on coordinate j then gives the integrand
  # r_kl * D D' * lambda_j*abs(ell_j) / sum_i(lambda_i*ell_i^2).
  # Importance sampling the (d-1)-dimensional fitted mixture marginal
  # divides this by its marginal density. Summing ALL coordinate charts
  # reproduces the original boundary integral. There is no bandwidth and
  # no fixed integration window. The positive outer products preserve PSD.
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
  cores <- lapply(cores, install_core)
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

  fit_once <- function(e, X, y, labelled, classes, ds, rep_id, method, budget = NA_real_) {
    monitor <- e$.__cov_monitor; monitor$active_calls <- 0L; monitor$clipped_eigenvalues <- 0L
    base <- previous$records[[paste0(ds, "_", rep_id)]]
    if (method == "Pilot") {
      warm <- base$pilot_fit
      base_labels <- base$pilot
    } else {
      old <- base$final[[paste0(method, "_", budget)]]
      base_labels <- c(base$pilot, old$selected)
      warm <- if (identical(sort(which(labelled)), sort(base_labels))) old$fit else NULL
    }
    if (!identical(sort(which(labelled)), sort(base_labels)) && method == "Pilot") stop("Pilot labels disagree with the audited run.")
    initial_warm <- !is.null(warm); initial_fit <- warm
    initial_iterations <- if (initial_warm) as.integer(warm$iterations) else 0L
    used_iterations <- initial_iterations
    for (level in seq_along(protocol$EM_PARAMETER_TOLS)) {
      remaining <- protocol$EM_MAXIT - used_iterations
      if (remaining < protocol$EM_STABLE_ITERATIONS + 1L) stop("EM exhausted its iteration limit before stationarity.")
      parameter_tol <- protocol$EM_PARAMETER_TOLS[level]
      fit <- e$fit_ss_qda(X, y, labelled, classes, maxit = remaining, tol = protocol$EM_TOL,
        start_fit = warm, parameter_tol = parameter_tol, stable_required = protocol$EM_STABLE_ITERATIONS)
      used_iterations <- used_iterations + fit$additional_iterations
      lj <- e$compute_log_joint(X, fit); lab <- which(labelled); unl <- which(!labelled)
      fresh_ll <- sum(lj[cbind(lab, match(as.character(y[lab]), classes))]) +
        if (length(unl)) sum(e$log_sum_exp_rows(lj[unl, , drop = FALSE])) else 0
      mismatch <- abs(fresh_ll - fit$logLik)
      min_increment <- min(diff(fit$logLik_history))
      add(ds, rep_id, "fit", "converged", fit$converged, method, budget, paste0("level=", level))
      add(ds, rep_id, "fit", "fresh_logLik_discrepancy", mismatch, method, budget, paste0("level=", level))
      add(ds, rep_id, "fit", "minimum_logLik_increment", min_increment, method, budget, paste0("level=", level))
      add(ds, rep_id, "fit", "maximum_scaled_parameter_change", fit$max_scaled_parameter_change, method, budget, paste0("level=", level))
      add(ds, rep_id, "fit", "stable_joint_iterations", fit$stable_joint_iterations, method, budget, paste0("level=", level))
      add(ds, rep_id, "fit", "parameter_tolerance", parameter_tol, method, budget, paste0("level=", level))
      if (!isTRUE(fit$converged) || fit$stable_joint_iterations < protocol$EM_STABLE_ITERATIONS ||
          fit$max_scaled_parameter_change >= parameter_tol || !is.finite(fresh_ll) ||
          mismatch > 1e-9 + 64 * .Machine$double.eps * (1 + abs(fresh_ll)) ||
          min_increment < -1e-10 * (1 + max(abs(fit$logLik_history)))) stop("Joint EM integrity check failed for ", method, ".")
      for (S in fit$Sigma) {
        if (min(eigen(sym(S), symmetric = TRUE, only.values = TRUE)$values) <=
            1.000001 * 1e-6 * max(mean(diag(S)), 1e-10)) stop("A final covariance is on its numerical floor; unconstrained stationarity requires review.")
      }
      diagnostic <- observed_stationarity(X, y, labelled, classes, fit)
      add(ds, rep_id, "stationarity", "minimum_scaled_observed_Hessian_eigenvalue",
        diagnostic$minimum_scaled_observed_Hessian_eigenvalue, method, budget, paste0("level=", level))
      add(ds, rep_id, "stationarity", "maximum_scaled_Newton_step", diagnostic$maximum_scaled_Newton_step, method, budget, paste0("level=", level))
      add(ds, rep_id, "stationarity", "maximum_absolute_score_per_observation", diagnostic$maximum_absolute_score_per_observation,
        method, budget, paste0("level=", level))
      if (diagnostic$minimum_scaled_observed_Hessian_eigenvalue <= 1e-8) stop("The observed likelihood Hessian is not reliably positive definite; this fit needs review.")
      if (diagnostic$maximum_scaled_Newton_step <= protocol$EM_NEWTON_TOL) break
      warm <- fit
    }
    if (diagnostic$maximum_scaled_Newton_step > protocol$EM_NEWTON_TOL) stop("Observed-likelihood stationarity failed after all parameter tolerances.")
    fit$stationarity <- diagnostic[c("minimum_scaled_observed_Hessian_eigenvalue", "maximum_scaled_Newton_step", "maximum_absolute_score_per_observation")]
    fit$continuation_initially_used <- initial_warm
    fit$validation_07_iterations <- used_iterations - initial_iterations
    add(ds, rep_id, "fit", "continued_identical_label_fit", initial_warm, method, budget)
    add(ds, rep_id, "fit", "validation_07_EM_iterations", fit$validation_07_iterations, method, budget)
    if (initial_warm) {
      add(ds, rep_id, "fit_sensitivity", "scaled_parameter_change_from_validation_06",
        compare_fit(fit, initial_fit), method, budget, "identical labeled set")
      add(ds, rep_id, "fit_sensitivity", "logLik_gain_from_validation_06",
        fit$logLik - initial_fit$logLik, method, budget, "identical labeled set")
    }
    add(ds, rep_id, "fit", "covariance_floor_active_calls", monitor$active_calls, method, budget)
    add(ds, rep_id, "fit", "covariance_clipped_eigenvalues", monitor$clipped_eigenvalues, method, budget)
    fit
  }

  compare_fit <- function(fit, ref) max(abs(c(fit$pi, fit$mu, unlist(fit$Sigma)) -
    c(ref$pi, ref$mu, unlist(ref$Sigma))) / (1 + abs(c(ref$pi, ref$mu, unlist(ref$Sigma)))))
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
  solve_one <- function(e, score, labelled, candidate, target, B, ds) {
    sol <- e$solve_finite_pool_design(score, labelled, candidate, target, B,
      maxit = if (ds == "Landsat") 1000L else 1200L, tol = 1e-5)
    certify_design(e, score, labelled, candidate, target, B, sol)
  }

  train_raw <- as.matrix(read.table(paths$trn, header = FALSE)); test_raw <- as.matrix(read.table(paths$tst, header = FALSE))
  stopifnot(identical(dim(train_raw), c(4435L, 37L)), identical(dim(test_raw), c(2000L, 37L)))
  center <- colMeans(train_raw[, 17:20, drop = FALSE]); scale_values <- apply(train_raw[, 17:20, drop = FALSE], 2L, sd)
  landsat_data <- list(X_train = sweep(sweep(train_raw[, 17:20, drop = FALSE], 2L, center, "-"), 2L, scale_values, "/"),
    X_test = sweep(sweep(test_raw[, 17:20, drop = FALSE], 2L, center, "-"), 2L, scale_values, "/"),
    y_train = factor(as.character(train_raw[, 37L]), levels = c("1", "2", "3", "4", "5", "7")),
    y_test = factor(as.character(test_raw[, 37L]), levels = c("1", "2", "3", "4", "5", "7")),
    classes = c("1", "2", "3", "4", "5", "7"), preprocessing = list(center = center, scale = scale_values),
    train_indices = seq_len(4435L), test_indices = seq_len(2000L))
  colnames(landsat_data$X_train) <- colnames(landsat_data$X_test) <- paste0("Band", 1:4)
  bean <- as.data.frame(readxl::read_excel(paths$xlsx)); bean$Class <- factor(bean$Class)
  bean_source_rows <- which(!duplicated(bean)); bean <- bean[bean_source_rows, , drop = FALSE]
  stopifnot(nrow(bean) == 13543L, length(bean_source_rows) == 13543L)
  bean_X <- as.matrix(bean[, setdiff(names(bean), "Class"), drop = FALSE])
  data_for <- function(ds, rep_id) {
    if (ds == "Landsat") return(landsat_data)
    e <- cores$DryBean; sp <- e$stratified_split(bean$Class, 0.70, 20260901L + 10000L * rep_id)
    prep <- e$fit_preprocess(bean_X[sp$train, , drop = FALSE], 5L)
    list(X_train = e$apply_preprocess(bean_X[sp$train, , drop = FALSE], prep),
      X_test = e$apply_preprocess(bean_X[sp$test, , drop = FALSE], prep),
      y_train = factor(bean$Class[sp$train], levels = levels(bean$Class)),
      y_test = factor(bean$Class[sp$test], levels = levels(bean$Class)), classes = levels(bean$Class),
      preprocessing = prep, train_indices = sp$train, test_indices = sp$test)
  }
  fixture_cases <- list(
    list(ds = "Landsat", rep = 1L, fit = "Pilot", expected = c(0.093225679319841631, 0.00028522543851156799, 7.3935125488331995e-05)),
    list(ds = "DryBean", rep = 1L, fit = "Pilot", expected = c(0.10476111955580666, 0.00040892122671347798, 3.6549232707527447e-05)),
    list(ds = "Landsat", rep = 29L, fit = "Pilot", expected = c(0.0024867889601560529, 0.0089103715307649635, 6.9414686785716306e-05)),
    list(ds = "Landsat", rep = 18L, fit = "Entropy_0.1", expected = c(0.013362431697222165, 0.0017827209690248983, 0.00012090468656252847)),
    list(ds = "Landsat", rep = 47L, fit = "Entropy_0.1", expected = c(0.0056422785559519613, 0.0024206208938609671, 7.5395213934383547e-05)),
    list(ds = "DryBean", rep = 40L, fit = "Margin_0.3", expected = c(0.13564001131601078, 0.00089621487666573241, 2.0556403618977163e-05)))
  preflight <- list()
  for (ci in seq_along(fixture_cases)) {
    f <- fixture_cases[[ci]]; dat <- data_for(f$ds, f$rep)
    old <- previous$records[[paste0(f$ds, "_", f$rep)]]
    labels <- old$pilot
    if (f$fit == "Pilot") pilot_fit <- old$pilot_fit else {
      pilot_fit <- old$final[[f$fit]]$fit; labels <- c(labels, old$final[[f$fit]]$selected)
    }
    labelled <- rep(FALSE, nrow(dat$X_train)); labelled[labels] <- TRUE
    z <- observed_stationarity(dat$X_train, dat$y_train, labelled, dat$classes, pilot_fit)
    values <- c(z$minimum_scaled_observed_Hessian_eigenvalue, z$maximum_scaled_Newton_step,
      z$maximum_absolute_score_per_observation)
    discrepancy <- max(abs(values - f$expected) / (1 + abs(f$expected)))
    if (discrepancy > 2e-7) stop("Observed-likelihood preflight disagrees with independent Python calculations.")
    score_discrepancy <- NA_real_
    if (f$rep == 1L && f$fit == "Pilot") {
      e <- cores[[f$ds]]; index <- e$build_param_index(length(dat$classes), ncol(dat$X_train))
      original_score <- e$.__original_score_components(dat$X_train, pilot_fit, index)
      new_score <- score_object(dat$X_train, make_model(pilot_fit))
      original_values <- c(original_score$tau, original_score$bar, unlist(original_score$class_scores))
      new_values <- c(new_score$tau, new_score$bar, unlist(new_score$class_scores))
      score_discrepancy <- max(abs(new_values - original_values) / (1 + abs(original_values)))
      if (score_discrepancy > 1e-8) stop("Vectorized scores disagree with the reviewed score implementation.")
    }
    preflight[[ci]] <- data.frame(Dataset = f$ds, Replication = f$rep, Fit = f$fit,
      IndependentStationarityDiscrepancy = discrepancy, OriginalScoreDiscrepancy = score_discrepancy)
  }
  preflight <- do.call(rbind, preflight)
  write.csv(preflight, file.path(out, "validation_07_preflight.csv"), row.names = FALSE)
  cat("\nIndependent stationarity and score preflight passed:\n"); print(preflight)

  atomic_save <- function(object, path) {
    if (file.exists(path)) stop("A completed replication already exists: ", path)
    temporary <- paste0(path, ".partial"); saveRDS(object, temporary)
    if (!file.rename(temporary, path)) stop("Cannot commit replication file: ", path)
  }
  run_replication <- function(ds, rep_id) {
    current_checks <<- list(); e <- cores[[ds]]; dat <- data_for(ds, rep_id)
    X <- dat$X_train; y <- dat$y_train; n <- nrow(X); d <- ncol(X); classes <- dat$classes
    pilot_n <- ceiling(0.05 * n); total_n <- ceiling(protocol$budgets * n)
    set.seed(if (ds == "Landsat") 20260902L + 10000L * rep_id else 20260901L + 20000L * rep_id)
    pilot <- sample(seq_len(n), pilot_n, replace = FALSE); labelled <- rep(FALSE, n); labelled[pilot] <- TRUE
    candidate <- which(!labelled); counts <- table(factor(y[pilot], levels = classes))
    cat("\n", ds, ": replication ", rep_id, "/", n_rep, "; pilot counts ", paste(as.integer(counts), collapse = ","), "\n", sep = "")
    if (any(counts <= d + 1L)) stop("Insufficient pilot class counts; this replication must be reviewed, not replaced or omitted.")
    ref <- previous$records[[paste0(ds, "_", rep_id)]]
    if (!identical(pilot, ref$pilot) || !identical(dat$train_indices, ref$train_indices) ||
        !identical(dat$test_indices, ref$test_indices) ||
        max(abs(unlist(dat$preprocessing) - unlist(ref$preprocessing))) > 1e-10)
      stop("Pilot, split, or preprocessing disagrees with the audited experiment.")
    add(ds, rep_id, "input_replay", "pilot_split_preprocessing_match", TRUE)
    fit <- fit_once(e, X, y, labelled, classes, ds, rep_id, "Pilot")
    score <- e$compute_score_components(X, fit, e$build_param_index(length(classes), d))
    target_fisher <- e$weighted_E_tt(score, rep(1, n)); model <- make_model(fit)
    lj_discrepancy <- max(abs(log_joint(X, model) - score$log_joint) / (1 + abs(score$log_joint)))
    add(ds, rep_id, "model", "scaled_log_joint_discrepancy", lj_discrepancy)
    if (lj_discrepancy > 1e-9) stop("Fitted density and boundary model disagree.")
    fisher <- lapply(seq_along(total_n), function(bi) solve_one(e, score, labelled, candidate, target_fisher, total_n[bi] - pilot_n, ds))
    adaptive <- NULL; precision <- NULL; agg <- NULL; seeds <- NULL; accepted_level <- NA_integer_
    for (level in seq_along(protocol$mc_sizes)) {
      size <- protocol$mc_sizes[level]; estimates <- list(); seeds <- integer(protocol$mc_reps)
      for (r in seq_len(protocol$mc_reps)) {
        seeds[r] <- protocol$mc_seed_base + 100000L * unname(protocol$mc_case_id[ds])
        seeds[r] <- seeds[r] + 3000L + 1000L * (rep_id - 1L) + 10L * r + 200L * (level - 1L)
        cat("  Boundary n/chart=", size, "; integration replicate ", r, "/", protocol$mc_reps, "\n", sep = "")
        z <- atlas(model, size, seeds[r]); estimates[[r]] <- z$H
        add(ds, rep_id, "boundary_integrity", "maximum_scaled_root_residual", z$max_scaled_root_residual,
          detail = paste0("level=", level, "; replicate=", r, "; seed=", seeds[r]))
      }
      agg <- aggregate_H(estimates)
      add(ds, rep_id, "boundary_precision", "relative_H_MCSE", agg$relative_MCSE, detail = paste0("level=", level, "; n_chart=", size))
      if (level == 1L) {
        add(ds, rep_id, "fit_sensitivity", "relative_curvature_change_from_validation_06",
          relative(agg$H, ref$H))
      }
      if (agg$relative_MCSE > protocol$max_H_MCSE) next
      adaptive <- lapply(seq_along(total_n), function(bi) solve_one(e, score, labelled, candidate, agg$H, total_n[bi] - pilot_n, ds))
      precision <- lapply(adaptive, function(design) precision_at_design(agg, design))
      adequate <- all(vapply(precision, function(z) z$relative_white_H_MCSE <= protocol$max_white_H_MCSE &&
        z$relative_objective_MCSE <= protocol$max_objective_MCSE, logical(1)))
      for (bi in seq_along(total_n)) {
        add(ds, rep_id, "boundary_precision", "relative_information_weighted_H_MCSE", precision[[bi]]$relative_white_H_MCSE,
          "AdaptiveRisk", protocol$budgets[bi], paste0("level=", level))
        add(ds, rep_id, "boundary_precision", "relative_objective_MCSE", precision[[bi]]$relative_objective_MCSE,
          "AdaptiveRisk", protocol$budgets[bi], paste0("level=", level))
      }
      if (adequate) {accepted_level <- level; break}
    }
    if (is.na(accepted_level)) stop("Integration precision failed at the maximum effort; this replication needs review.")
    entropy <- e$score_entropy(score$tau); margin <- e$score_margin(score$tau)
    finals <- list(); raws <- list(); class_rows <- list(); design_rows <- list(); curv_rows <- list(); rr <- 1L
    for (bi in seq_along(total_n)) {
      budget <- protocol$budgets[bi]; B <- total_n[bi] - pilot_n
      set.seed(if (ds == "Landsat") 20260902L + 20000L * rep_id + bi else 20260901L + 30000L * rep_id + bi)
      selections <- list(Random = sample(candidate, B, replace = FALSE),
        Entropy = candidate[order(entropy[candidate], decreasing = TRUE)[seq_len(B)]],
        Margin = candidate[order(margin[candidate], decreasing = TRUE)[seq_len(B)]],
        Fisher = fisher[[bi]]$selected, AdaptiveRisk = adaptive[[bi]]$selected)
      for (method in METHODS) {
        selected <- selections[[method]]
        if (length(selected) != B || anyDuplicated(selected) || length(intersect(selected, pilot))) stop("Invalid acquired set.")
        lab <- labelled; lab[selected] <- TRUE
        old <- ref$final[[paste0(method, "_", budget)]]
        count <- length(setdiff(selected, old$selected))
        add(ds, rep_id, "fit_sensitivity", "acquired_points_replaced_from_validation_06", count, method, budget)
        final_fit <- fit_once(e, X, y, lab, classes, ds, rep_id, method, budget)
        pred <- e$predict_ss_qda(final_fit, dat$X_test)
        evaluation <- e$evaluate_fit(final_fit, dat$X_test, dat$y_test, classes)
        disagreements <- sum(as.character(pred$class) != as.character(old$predicted))
        add(ds, rep_id, "fit_sensitivity", "test_predictions_changed_from_validation_06", disagreements, method, budget)
        key <- paste0(method, "_", budget)
        finals[[key]] <- list(selected = selected, fit = final_fit, predicted = as.character(pred$class),
          error = evaluation$error, balanced_error = evaluation$balanced_error)
        raws[[rr]] <- data.frame(Dataset = ds, Replication = rep_id, Budget = budget, Method = method,
          NLabelled = sum(lab), Error = evaluation$error, BalancedError = evaluation$balanced_error,
          LogLik = final_fit$logLik, EMIterations = final_fit$iterations, EMConverged = final_fit$converged)
        class_rows[[rr]] <- data.frame(Dataset = ds, Replication = rep_id, Budget = budget, Method = method,
          Class = classes, ClassError = as.numeric(evaluation$class_error))
        cat("  ", method, "; budget=", budget, "; error=", sprintf("%.6f", evaluation$error), "; labels=", sum(lab), "\n", sep = "")
        rr <- rr + 1L
      }
      for (method in c("Fisher", "AdaptiveRisk")) {
        design <- if (method == "Fisher") fisher[[bi]] else adaptive[[bi]]
        design_rows[[length(design_rows) + 1L]] <- data.frame(Dataset = ds, Replication = rep_id, Budget = budget, Method = method,
          AdditionalBudget = B, Iterations = design$iterations, RelativeFWGap = design$relative_gap,
          RelaxedObjective = design$objective, RoundedObjective = design$rounded_objective,
          RoundingLossPct = 100 * (design$rounded_objective - design$objective) / design$objective,
          FractionalCount = sum(design$a > 1e-8 & design$a < 1 - 1e-8),
          MinimumInformationEigenvalue = design$minimum_information_eigenvalue)
      }
      p <- precision[[bi]]
      curv_rows[[bi]] <- data.frame(Dataset = ds, Replication = rep_id, Budget = budget,
        NPerChart = protocol$mc_sizes[accepted_level], IntegrationReplicates = protocol$mc_reps,
        RelativeHMCSE = agg$relative_MCSE, RelativeInformationWeightedHMCSE = p$relative_white_H_MCSE,
        RelativeObjectiveMCSE = p$relative_objective_MCSE, NumericalPrecisionPassed = TRUE)
    }
    compact_design <- function(z) z[c("a", "selected", "objective", "rounded_objective", "relative_gap", "iterations")]
    record <- list(Dataset = ds, Replication = rep_id, protocol_version = protocol$version,
      pilot = pilot, pilot_fit = fit, preprocessing = dat$preprocessing, train_indices = dat$train_indices, test_indices = dat$test_indices,
      H = agg$H, mc_estimates = agg$estimates, mc_norm_sq = vapply(agg$estimates, function(H) sum(H^2), numeric(1)),
      mc_seeds = seeds, mc_n = protocol$mc_sizes[accepted_level], mc_relative_se = agg$relative_MCSE,
      precision = precision, fisher = lapply(fisher, compact_design), adaptive = lapply(adaptive, compact_design),
      final = finals, raw = do.call(rbind, raws), class_errors = do.call(rbind, class_rows), designs = do.call(rbind, design_rows),
      curvature = do.call(rbind, curv_rows), checks = do.call(rbind, current_checks),
      entropy = entropy, margin = margin)
    atomic_save(record, file.path(out, "fits", sprintf("%s_replication_%03d.rds", ds, rep_id)))
    record
  }

  for (ds in datasets) for (rep_id in seq_len(n_rep)) {
    fit_path <- file.path(out, "fits", sprintf("%s_replication_%03d.rds", ds, rep_id))
    if (file.exists(fit_path)) {
      record <- readRDS(fit_path)
      if (!identical(record$protocol_version, protocol$version) || record$Dataset != ds || record$Replication != rep_id)
        stop("Resume record does not match this protocol.")
      cat("Resuming: retained completed ", ds, " replication ", rep_id, "\n", sep = "")
    } else {
      record <- tryCatch(withCallingHandlers(run_replication(ds, rep_id), warning = function(w) {
        add(ds, rep_id, "execution", "warning", conditionMessage(w))
        cat("WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning")
      }), error = function(err) {
        failures[[length(failures) + 1L]] <<- data.frame(Dataset = ds, Replication = rep_id, Reason = conditionMessage(err))
        checkpoint(); cat("\nRUN STOPPED: ", conditionMessage(err), "\nOutput: ", out, "\n", sep = "")
        stop(err)
      })
    }
    append_record(record); checkpoint(); rm(record); gc(verbose = FALSE)
  }
  merged <- merge_tables(); raw <- merged$raw
  expected_rows <- length(datasets) * n_rep * length(METHODS) * length(protocol$budgets)
  if (nrow(raw) != expected_rows || anyDuplicated(raw[, c("Dataset", "Replication", "Budget", "Method")]))
    stop("Raw results do not contain exactly the complete requested experiment.")
  summaries <- list(); contrasts <- list(); sr <- 1L; cr <- 1L
  for (ds in datasets) for (budget in protocol$budgets) for (outcome in c("Error", "BalancedError")) {
    z <- raw[raw$Dataset == ds & raw$Budget == budget, ]
    for (method in METHODS) {
      values <- z[z$Method == method, outcome]
      summaries[[sr]] <- data.frame(Dataset = ds, Budget = budget, Outcome = outcome, Method = method,
        N = length(values), Mean = mean(values), MCSE = if (length(values) > 1L) sd(values) / sqrt(length(values)) else NA_real_)
      sr <- sr + 1L
    }
    a <- z[z$Method == "AdaptiveRisk", c("Replication", outcome)]
    for (method in setdiff(METHODS, "AdaptiveRisk")) {
      b <- z[z$Method == method, c("Replication", outcome)]
      b <- b[match(a$Replication, b$Replication), ]; difference <- a[[outcome]] - b[[outcome]]
      nn <- length(difference); se <- if (nn > 1L) sd(difference) / sqrt(nn) else NA_real_
      critical <- if (nn > 1L) qt(0.975, nn - 1L) else NA_real_
      p <- if (is.finite(se) && se > 0) 2 * pt(-abs(mean(difference) / se), nn - 1L) else NA_real_
      contrasts[[cr]] <- data.frame(Dataset = ds, Budget = budget, Outcome = outcome,
        Contrast = paste0("AdaptiveRisk-minus-", method), N = nn, MeanDifference = mean(difference), MCSE = se,
        Lower95 = mean(difference) - critical * se, Upper95 = mean(difference) + critical * se, PValue = p)
      cr <- cr + 1L
    }
  }
  summary <- do.call(rbind, summaries); paired <- do.call(rbind, contrasts); paired$HolmPValue <- NA_real_
  for (ds in datasets) for (outcome in c("Error", "BalancedError")) {
    ii <- which(paired$Dataset == ds & paired$Outcome == outcome)
    paired$HolmPValue[ii] <- p.adjust(paired$PValue[ii], method = "holm")
  }
  write.csv(summary, file.path(out, "realdata_final_summary.csv"), row.names = FALSE)
  write.csv(paired, file.path(out, "realdata_final_paired.csv"), row.names = FALSE)
  status <- data.frame(Dataset = datasets, RequestedReplications = n_rep,
    CompleteReplications = vapply(datasets, function(ds) length(unique(raw$Replication[raw$Dataset == ds])), integer(1)),
    RawRows = vapply(datasets, function(ds) sum(raw$Dataset == ds), integer(1)), ExecutionComplete = TRUE)
  write.csv(status, file.path(out, "realdata_final_status.csv"), row.names = FALSE)
  review_records <- list()
  for (ds in datasets) for (rep_id in seq_len(n_rep)) {
    record <- readRDS(file.path(out, "fits", sprintf("%s_replication_%03d.rds", ds, rep_id)))
    record$mc_estimates <- NULL; record$entropy <- NULL; record$margin <- NULL
    review_records[[paste0(ds, "_", rep_id)]] <- record
  }
  saveRDS(list(protocol = protocol, status = status, preflight = preflight, tables = merged, summary = summary, paired = paired,
    bean_source_rows = bean_source_rows, records = review_records), file.path(out, "realdata_final_review.rds"))
  checkpoint(); cat("\nCompleted requested numerical run:\n"); print(status)
  cat("\nReturn these four files:\n")
  cat(file.path(out, c("realdata_final_review.rds", "realdata_final_raw.csv", "realdata_final_checks.csv", "realdata_final_log.txt")), sep = "\n")
  cat("\nKeep the fits directory and other CSV files locally for independent audits and figure reproduction.\n")
  cat("Negative paired differences favor AdaptiveRisk. Intervals are pointwise; Holm adjustment covers 12 comparisons within each dataset/outcome.\n")
  cat("Convergence checks assess local likelihood stationarity; they do not prove a global likelihood maximum.\n")
  cat("Completion does not certify manuscript claims, simulations, or figures.\n")
  invisible(out)
}
