# Regenerate simulation figures from the frozen, fully reviewed R09 run.
# Does not rerun EM, acquisition performance replications, or curvature integration.
# Performance intervals are 95%, using t with 99 degrees of freedom.

validation_10_simulation_figures <- function(project_root, validation_09_dir) {
  stopifnot(dir.exists(project_root), dir.exists(validation_09_dir))
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  validation_09_dir <- normalizePath(validation_09_dir, winslash = "/", mustWork = TRUE)
  for (package in c("ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install the package: ", package)
  }
  protected <- c("pi", "t", "chol", "eigen", "sample.int", "set.seed", "order", "max.col", "pnorm", "dnorm")
  shadows <- protected[vapply(protected, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(shadows)) stop("Remove global numerical-function shadows: ", paste(shadows, collapse = ", "))
  path <- file.path(validation_09_dir, "simulation_09_review.rds")
  if (!file.exists(path)) stop("Cannot find simulation_09_review.rds in validation_09_dir.")
  if (unname(tools::md5sum(path)) != "52b6b38539d5342dbd995d2c0c56f58b") {
    stop("The review RDS differs from the fully audited 100-replication run. Do not substitute another run.")
  }
  review <- readRDS(path)
  if (review$protocol$version != "2026-10-04-simulation-joint-v1" ||
      review$protocol$implementation_md5 != "2360ff67a4eb631bdbdce6121f376b2b" ||
      length(review$records) != 200L || any(review$status$CompleteReplications != 100L) ||
      !all(review$status$ManuscriptReplicationTargetMet)) stop("R09 protocol/completeness check failed.")
  out <- tempfile("simulation_figures_", tmpdir = project_root)
  if (!dir.create(out)) stop("Cannot create an output folder.")
  initial_rng <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  initial_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  sink(file.path(out, "simulation_10_log.txt"), split = TRUE)
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(initial_rng)))
    if (had_seed) assign(".Random.seed", initial_seed, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    sink()
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  cat("Simulation figure reproduction, revision 2; ", format(Sys.time()), "\n", sep = "")
  cat("Output: ", out, "\n", sep = "")
  cat("No original script, result, or manuscript is overwritten.\n")
  cat("Performance: all 100 replications, all budgets, all six methods; 95% t intervals.\n")
  cat("Geometry: new feature-only true-model pools, 5000 observations and 1000 acquisitions.\n")
  cat("Geometry is an illustration, not a selected performance replication.\n")
  checks <- list()
  add_check <- function(name, value, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(Check = name, Value = value, Detail = detail)
  }
  datasets <- c("Scenario1", "Scenario2")
  methods <- c("Random", "Entropy", "Margin", "Fisher", "AdaptiveRisk", "OracleRisk")
  budgets <- c(.1, .2, .3)
  raw <- do.call(rbind, lapply(review$records, function(record) record$raw))
  if (nrow(raw) != 3600L || anyDuplicated(raw[c("Dataset", "Replication", "Budget", "Method")]) ||
      any(!raw$EMConverged) || any(!is.finite(raw$Error)) || any(!is.finite(raw$BalancedError))) {
    stop("The complete raw-result integrity check failed.")
  }
  summaries <- list(); comparisons <- list(); si <- 1L; ci <- 1L
  for (ds in datasets) for (budget in budgets) for (outcome in c("Error", "BalancedError")) {
    rows <- raw[raw$Dataset == ds & raw$Budget == budget, ]
    for (method in methods) {
      group <- rows[rows$Method == method, ]
      if (nrow(group) != 100L || !all(sort(group$Replication) == 1:100)) stop("Incomplete method/budget group.")
      values <- group[[outcome]]; se <- stats::sd(values) / 10
      margin <- stats::qt(.975, 99L) * se
      summaries[[si]] <- data.frame(Dataset = ds, Budget = budget, Outcome = outcome, Method = method,
        N = 100L, Mean = mean(values), MCSE = se, Lower95 = mean(values) - margin, Upper95 = mean(values) + margin)
      si <- si + 1L
    }
    adaptive <- rows[rows$Method == "AdaptiveRisk", c("Replication", outcome)]
    for (method in setdiff(methods, "AdaptiveRisk")) {
      comparator <- rows[rows$Method == method, c("Replication", outcome)]
      comparator <- comparator[match(adaptive$Replication, comparator$Replication), ]
      delta <- adaptive[[outcome]] - comparator[[outcome]]
      se <- stats::sd(delta) / 10; margin <- stats::qt(.975, 99L) * se
      pv <- if (se > 0) 2 * stats::pt(-abs(mean(delta) / se), 99L) else
        if (all(delta == 0)) 1 else NA_real_
      comparisons[[ci]] <- data.frame(Dataset = ds, Budget = budget, Outcome = outcome,
        Contrast = paste0("AdaptiveRisk-minus-", method), N = 100L, MeanDifference = mean(delta),
        MCSE = se, Lower95 = mean(delta) - margin, Upper95 = mean(delta) + margin, PValue = pv)
      ci <- ci + 1L
    }
  }
  summary <- do.call(rbind, summaries); paired <- do.call(rbind, comparisons)
  paired$HolmPValue <- NA_real_
  for (ds in datasets) for (outcome in c("Error", "BalancedError")) {
    ii <- which(paired$Dataset == ds & paired$Outcome == outcome)
    paired$HolmPValue[ii] <- stats::p.adjust(paired$PValue[ii], "holm")
  }
  compare_tables <- function(actual, cached, keys, columns) {
    key <- function(x) do.call(paste, c(x[keys], sep = "|"))
    index <- match(key(actual), key(cached))
    if (anyNA(index) || nrow(actual) != nrow(cached)) stop("Cached summary keys do not match.")
    error <- max(abs(as.matrix(actual[columns]) - as.matrix(cached[index, columns])), na.rm = TRUE)
    if (!is.finite(error) || error > 1e-12) stop("Cached statistical reproduction failed: ", error)
    error
  }
  err <- compare_tables(summary, review$summary, c("Dataset", "Budget", "Outcome", "Method"), c("Mean", "MCSE"))
  add_check("summary_maximum_absolute_difference", err)
  err <- compare_tables(paired, review$paired, c("Dataset", "Budget", "Outcome", "Contrast"),
    c("MeanDifference", "MCSE", "Lower95", "Upper95", "PValue", "HolmPValue"))
  add_check("paired_statistics_maximum_absolute_difference", err)
  write.csv(summary, file.path(out, "simulation_10_summary.csv"), row.names = FALSE)
  write.csv(paired, file.path(out, "simulation_10_paired.csv"), row.names = FALSE)
  write.csv(raw, file.path(out, "simulation_10_raw.csv"), row.names = FALSE)
  reference <- do.call(rbind, lapply(review$scenarios, function(scenario) scenario$reference))
  if (max(abs(reference$BayesQuadratureRisk[match(datasets, reference$Dataset)] -
      c(.148657029436564, .156159498347398))) > 1e-9) stop("Population Bayes references do not match.")
  write.csv(reference, file.path(out, "simulation_10_bayes_reference.csv"), row.names = FALSE)
  cat("All 72 mean/MCSE rows and 60 paired comparisons reproduced.\n")

  sym <- function(M) (M + t(M)) / 2
  inverse <- function(M) chol2inv(chol(sym(M)))
  norm_f <- function(M) sqrt(sum(M^2))
  row_max <- function(M) do.call(pmax, lapply(seq_len(ncol(M)), function(k) M[, k]))
  model <- function(par) {
    S <- if (is.list(par$Sigma)) par$Sigma else lapply(1:3, function(k) par$Sigma[, , k])
    P <- lapply(S, inverse)
    b <- t(vapply(1:3, function(k) as.vector(P[[k]] %*% par$mu[k, ]), numeric(2)))
    cc <- vapply(1:3, function(k) log(par$pi[k]) - log(2 * pi) -
      sum(log(diag(chol(S[[k]])))) - sum(par$mu[k, ] * b[k, ]) / 2, numeric(1))
    list(pi = par$pi, mu = par$mu, S = S, P = P, b = b, cc = cc, p = 17L,
      alpha = 1:2, block = lapply(1:3, function(k) 2L + (k - 1L) * 5L + 1:5))
  }
  log_joint <- function(X, m) {
    X <- as.matrix(X)
    result <- vapply(1:3, function(k) -rowSums((X %*% m$P[[k]]) * X) / 2 +
      as.vector(X %*% m$b[k, ]) + m$cc[k], numeric(nrow(X)))
    matrix(result, nrow(X), 3L)
  }
  score_object <- function(X, m) {
    lj <- log_joint(X, m); tau <- exp(lj - row_max(lj)); tau <- tau / rowSums(tau)
    A <- matrix(-m$pi[1:2], 3, 2, byrow = TRUE); A[1, 1] <- A[1, 1] + 1; A[2, 2] <- A[2, 2] + 1
    scores <- lapply(1:3, function(k) {
      U <- sweep(X, 2, m$mu[k, ], "-") %*% m$P[[k]]
      T <- matrix(0, nrow(X), 17L); T[, 1:2] <- matrix(A[k, ], nrow(X), 2, byrow = TRUE)
      T[, m$block[[k]]] <- cbind(U, (U[, 1]^2 - m$P[[k]][1, 1]) / 2,
        U[, 1] * U[, 2] - m$P[[k]][1, 2], (U[, 2]^2 - m$P[[k]][2, 2]) / 2)
      T
    })
    bar <- Reduce(`+`, lapply(1:3, function(k) scores[[k]] * tau[, k]))
    list(tau = tau, bar = bar, residual = lapply(scores, function(T) T - bar), n = nrow(X))
  }
  information <- function(s, a) {
    M <- crossprod(s$bar) / s$n
    for (k in 1:3) M <- M + crossprod(s$residual[[k]] * sqrt(a * s$tau[, k])) / s$n
    sym(M)
  }
  value <- function(s, G) {
    psi <- numeric(s$n)
    for (k in 1:3) psi <- psi + s$tau[, k] * rowSums((s$residual[[k]] %*% G) * s$residual[[k]])
    if (min(psi) < -1e-10 * max(1, max(psi))) stop("Negative classification value beyond roundoff.")
    pmax(psi, 0)
  }
  objective <- function(H, M) sum(H * t(inverse(M)))
  solve_design <- function(s, H, B, tolerance = 1e-5, maxit = 5000L) {
    a <- rep(B / s$n, s$n); M <- information(s, a)
    for (iteration in seq_len(maxit)) {
      inv <- inverse(M); phi <- sum(H * t(inv)); psi <- value(s, sym(inv %*% H %*% inv))
      vertex <- numeric(s$n); vertex[order(psi, decreasing = TRUE)[seq_len(B)]] <- 1
      gap <- sum((vertex - a) * psi) / s$n
      if (gap < -64 * .Machine$double.eps * max(1, phi)) stop("Invalid FW gap.")
      if (max(0, gap) / phi <= tolerance) break
      if (iteration == maxit) stop("Geometry design did not converge. No figure is certified.")
      D <- information(s, vertex) - M
      derivative <- function(step) {
        C <- inverse(M + step * D)
        -sum(H * t(C %*% D %*% C))
      }
      if (derivative(0) >= 0) stop("FW direction is not descending.")
      step <- if (derivative(1) <= 0) 1 else stats::uniroot(derivative, c(0, 1), tol = 1e-12)$root
      a <- a + step * (vertex - a); M <- information(s, a)
    }
    selected <- order(a, decreasing = TRUE)[seq_len(B)]
    binary <- numeric(s$n); binary[selected] <- 1
    rounded_phi <- objective(H, information(s, binary))
    list(a = a, selected = selected, M = M, objective = phi, rounded_objective = rounded_phi,
      relative_gap = max(0, gap) / phi, rounding_loss_pct = 100 * (rounded_phi - phi) / phi,
      iterations = iteration, fractional_count = sum(a > 1e-8 & a < 1 - 1e-8))
  }
  # Check the standalone coordinate order/information arithmetic against both
  # cached true-model Oracle designs from replication 1 at the 20% budget.
  for (ds in datasets) {
    r <- review$records[[paste0(ds, "_1")]]; H <- review$scenarios[[ds]]$H
    s <- score_object(r$training$Y, model(review$scenarios[[ds]]$truth))
    cached <- r$oracle[[2L]]; a <- numeric(s$n); a[r$pilot] <- 1
    a[setdiff(seq_len(s$n), r$pilot)] <- cached$a
    M <- information(s, a); difference <- norm_f(M - cached$M) / norm_f(cached$M)
    phi_error <- abs(objective(H, M) - cached$objective) / cached$objective
    if (max(difference, phi_error) > 1e-8) stop("Standalone cached-design arithmetic check failed for ", ds)
    add_check(paste0(ds, "_cached_information_relative_difference"), difference)
    add_check(paste0(ds, "_cached_objective_relative_difference"), phi_error)
  }
  cat("Standalone Gaussian information/coordinate-order checks passed.\n")
  base_theme <- ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 11), plot.tag = ggplot2::element_text(face = "bold"))
  save_plot <- function(plot, filename, width, height) {
    ggplot2::ggsave(file.path(out, paste0(filename, ".pdf")), plot, width = width, height = height,
      units = "in", device = "pdf", bg = "white", useDingbats = FALSE)
    ggplot2::ggsave(file.path(out, paste0(filename, ".png")), plot, width = width, height = height,
      units = "in", dpi = 300, bg = "white")
  }
  plotdata <- summary[summary$Outcome == "Error", ]
  plotdata$Method <- factor(plotdata$Method, methods)
  plotdata$Scenario <- factor(plotdata$Dataset, datasets, c("Scenario 1: curved QDA", "Scenario 2: near-linear QDA"))
  reference$Scenario <- factor(reference$Dataset, datasets, levels(plotdata$Scenario))
  palette <- c(Random = "#858585", Entropy = "#E69F00", Margin = "#009E73", Fisher = "#0072B2",
    AdaptiveRisk = "#D55E00", OracleRisk = "#CC79A7")
  labels <- c("Random", "Entropy", "Margin", "Fisher", "AdaptiveRisk", "Oracle")
  performance <- ggplot2::ggplot(plotdata, ggplot2::aes(Budget, Mean, color = Method, group = Method)) +
    ggplot2::geom_hline(data = reference, ggplot2::aes(yintercept = BayesQuadratureRisk),
      linetype = "dotted", color = "#303030", linewidth = .55) +
    ggplot2::geom_line(linewidth = .65, position = ggplot2::position_dodge(width = .013)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower95, ymax = Upper95), width = .004,
      linewidth = .45, position = ggplot2::position_dodge(width = .013)) +
    ggplot2::geom_point(size = 2.1, position = ggplot2::position_dodge(width = .013)) +
    ggplot2::facet_wrap(~Scenario, nrow = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = palette, breaks = methods, labels = labels, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = budgets, labels = c("10%", "20%", "30%")) +
    ggplot2::labs(x = "Total labeled proportion", y = "Test classification error", color = NULL,
      caption = "Means and pointwise 95% t intervals over 100 replications. Dotted lines: population Bayes error.") +
    base_theme + ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(size = 9)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))
  save_plot(performance, "Figure3_performance_combined_validated", 11.5, 4.8)
  write.csv(plotdata, file.path(out, "simulation_10_performance_plotdata.csv"), row.names = FALSE)

  boundary_data <- function(m, limits, n = 1600L) {
    result <- list(); rr <- 1L
    for (j in 1:2) {
      other <- 3L - j; grid <- seq(limits[other, 1], limits[other, 2], length.out = n)
      for (k in 1:2) for (l in (k + 1L):3L) {
        A <- (m$P[[l]] - m$P[[k]]) / 2; b <- m$b[k, ] - m$b[l, ]; cc <- m$cc[k] - m$cc[l]
        aa <- A[j, j]; bb <- 2 * grid * A[other, j] + b[j]
        cv <- grid^2 * A[other, other] + grid * b[other] + cc
        if (aa == 0) {
          ii <- which(bb != 0); roots <- list(list(i = ii, z = -cv[ii] / bb[ii]))
        } else {
          disc <- bb^2 - 4 * aa * cv; ii <- which(disc > 0)
          q <- -.5 * (bb[ii] + ifelse(bb[ii] >= 0, 1, -1) * sqrt(disc[ii]))
          roots <- list(list(i = ii, z = q / aa), list(i = ii, z = cv[ii] / q))
        }
        for (branch in seq_along(roots)) {
          root <- roots[[branch]]; if (!length(root$i)) next
          X <- matrix(0, length(root$i), 2L); X[, other] <- grid[root$i]; X[, j] <- root$z
          lj <- log_joint(X, m)
          active <- is.finite(root$z) & root$z >= limits[j, 1] & root$z <= limits[j, 2] &
            (lj[, k] + lj[, l]) / 2 >= row_max(lj) - 1e-9
          if (!any(active)) next
          X <- X[active, , drop = FALSE]; index <- root$i[active]
          group <- paste(j, k, l, branch, cumsum(c(TRUE, diff(index) > 1L)), sep = "_")
          result[[rr]] <- data.frame(x = X[, 1], y = X[, 2], group = group)
          rr <- rr + 1L
        }
      }
    }
    if (!length(result)) stop("No active Bayes boundary found in the display window.")
    do.call(rbind, result)
  }
  geometry_checks <- list()
  for (di in seq_along(datasets)) {
    ds <- datasets[di]; m <- model(review$scenarios[[ds]]$truth); H <- review$scenarios[[ds]]$H
    seed <- c(2026104401L, 2026104402L)[di]; set.seed(seed)
    component <- sample.int(3L, 5000L, replace = TRUE, prob = m$pi)
    X <- matrix(0, 5000L, 2L)
    for (k in 1:3) {
      ii <- which(component == k)
      Z <- matrix(stats::rnorm(length(ii) * 2L), length(ii), 2L) %*% chol(m$S[[k]])
      X[ii, ] <- sweep(Z, 2L, m$mu[k, ], "+")
    }
    s <- score_object(X, m); M_uniform <- information(s, rep(.2, nrow(X)))
    inv <- inverse(M_uniform); G <- sym(inv %*% H %*% inv)
    psi <- value(s, G); entropy <- -rowSums(s$tau * log(pmax(s$tau, 1e-300)))
    cat("START ", ds, " geometry; seed=", seed, "; pool=5000; budget=1000\n", sep = "")
    design <- solve_design(s, H, 1000L)
    entropy_selected <- order(entropy, decreasing = TRUE)[1:1000]
    risk <- seq_len(nrow(X)) %in% design$selected
    uncertain <- seq_len(nrow(X)) %in% entropy_selected
    group <- ifelse(risk & uncertain, "Shared", ifelse(risk, "Risk only", ifelse(uncertain, "Entropy only", "Unselected")))
    q99 <- as.numeric(stats::quantile(psi, .99, names = FALSE))
    if (!is.finite(q99) || q99 <= 0) stop("Classification-value display scale is invalid.")
    sd_matrix <- t(vapply(m$S, function(S) sqrt(diag(S)), numeric(2)))
    limits <- cbind(apply(m$mu - 4.5 * sd_matrix, 2, min), apply(m$mu + 4.5 * sd_matrix, 2, max))
    grid <- expand.grid(x = seq(limits[1, 1], limits[1, 2], length.out = 220L),
      y = seq(limits[2, 1], limits[2, 2], length.out = 220L))
    grid_score <- score_object(as.matrix(grid), m)
    grid$Class <- factor(max.col(grid_score$tau, ties.method = "first"), 1:3)
    grid$Entropy <- -rowSums(grid_score$tau * log(pmax(grid_score$tau, 1e-300))) / log(3)
    grid$ClassificationValue <- pmin(value(grid_score, G) / q99, 1)
    points <- data.frame(x = X[, 1], y = X[, 2], Entropy = entropy, ClassificationValue = psi,
      RelativeClassificationValue = pmin(psi / q99, 1), RelaxedRiskWeight = design$a,
      RiskSelected = risk, EntropySelected = uncertain, Group = group)
    edges <- boundary_data(m, limits)
    frame <- function(plot) plot + ggplot2::geom_path(data = edges,
      ggplot2::aes(x, y, group = group), inherit.aes = FALSE, color = "#303030", linewidth = .35) +
      ggplot2::coord_equal(xlim = limits[1, ], ylim = limits[2, ], expand = FALSE) +
      ggplot2::labs(x = expression(Y[1]), y = expression(Y[2])) + base_theme +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.key.height = grid::unit(.35, "cm"))
    pa <- frame(ggplot2::ggplot(grid, ggplot2::aes(x, y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = Class), alpha = .6) +
      ggplot2::scale_fill_manual(values = c("#C5DDEB", "#F4DFBF", "#CFDFD3"), name = "Bayes class") +
      ggplot2::labs(title = "Bayes regions and active boundaries"))
    pb <- frame(ggplot2::ggplot(grid, ggplot2::aes(x, y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = Entropy)) +
      ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Normalized entropy") +
      ggplot2::labs(title = "Posterior uncertainty"))
    pc <- frame(ggplot2::ggplot(grid, ggplot2::aes(x, y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = ClassificationValue)) +
      ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Relative value") +
      ggplot2::labs(title = "Classification value: uniform 20% information"))
    selected <- points[points$Group != "Unselected", ]
    selected$Group <- factor(selected$Group, c("Shared", "Risk only", "Entropy only"))
    pd <- frame(ggplot2::ggplot(points, ggplot2::aes(x, y)) +
      ggplot2::geom_point(color = "#B7B7B7", alpha = .28, size = .35) +
      ggplot2::geom_point(data = selected, ggplot2::aes(color = Group), alpha = .75, size = .55) +
      ggplot2::scale_color_manual(values = c(Shared = "#303030", "Risk only" = "#D55E00", "Entropy only" = "#0072B2"),
        name = NULL, drop = FALSE) + ggplot2::labs(title = "Exact acquisition sets: 1000 labels each"))
    combined <- patchwork::wrap_plots(pa, pb, pc, pd, ncol = 2) +
      patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
    filename <- if (ds == "Scenario1") "Figure1_geometry_validated" else "Scenario2_geometry_validated"
    save_plot(combined, filename, 11.2, 9.0)
    write.csv(points, file.path(out, paste0(ds, "_geometry_points.csv")), row.names = FALSE)
    write.csv(grid, file.path(out, paste0(ds, "_geometry_surfaces.csv")), row.names = FALSE)
    write.csv(edges, file.path(out, paste0(ds, "_active_boundaries.csv")), row.names = FALSE)
    record <- list(Dataset = ds, Seed = seed, N = 5000L, Budget = 1000L, truth = review$scenarios[[ds]]$truth,
      H = H, M_uniform = M_uniform, design = design, entropy_selected = entropy_selected,
      display_q99 = q99, limits = limits, points = points, surfaces = grid, boundaries = edges)
    saveRDS(record, file.path(out, paste0(ds, "_geometry.rds")))
    geometry_checks[[di]] <- data.frame(Dataset = ds, Seed = seed, N = 5000L, Budget = 1000L,
      Iterations = design$iterations, RelativeFWGap = design$relative_gap,
      RoundingLossPct = design$rounding_loss_pct, FractionalCount = design$fractional_count,
      Shared = sum(risk & uncertain), RiskOnly = sum(risk & !uncertain), EntropyOnly = sum(uncertain & !risk),
      EntropyUniformValueSpearman = stats::cor(entropy, psi, method = "spearman"))
    cat("END ", ds, "; FW gap=", format(design$relative_gap, digits = 9),
      "; rounding loss=", format(design$rounding_loss_pct, digits = 9), "%\n", sep = "")
  }
  write.csv(do.call(rbind, geometry_checks), file.path(out, "simulation_10_geometry_checks.csv"), row.names = FALSE)
  add_check("performance_confidence_level", .95, "Student-t with 99 degrees of freedom; pointwise, conditional on fixed tests")
  write.csv(do.call(rbind, checks), file.path(out, "simulation_10_checks.csv"), row.names = FALSE)
  print(sessionInfo())
  cat("\nCompleted figure generation. ZIP this entire NEW output folder and return the ZIP:\n", out, "\n", sep = "")
  cat("Do not overwrite the completed R09 folder. Do not select a different geometry seed based on appearance.\n")
  cat("Generated figures still require visual inspection and checking in the compiled manuscript.\n")
  invisible(out)
}
