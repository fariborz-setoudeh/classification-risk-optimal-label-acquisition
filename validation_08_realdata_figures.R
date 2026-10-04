# Regenerate real-data figures from the audited validation 07 fits.
# This script does not fit models, optimize designs, or acquire new labels.
# Source the file, then call validation_08_realdata_figures(...).

validation_08_realdata_figures <- function(project_root, validation_07_dir) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  validation_07_dir <- normalizePath(validation_07_dir, winslash = "/", mustWork = TRUE)
  review_path <- file.path(validation_07_dir, "realdata_final_review.rds")
  if (!file.exists(review_path) ||
      unname(tools::md5sum(review_path)) != "a826445bf5c260963f223a53a01ace87")
    stop("Use the audited validation 07 review cache, without resaving or modifying it.")
  for (package in c("readxl", "ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Your existing ", package, " installation is required.")
  }
  shadows <- c("pi", "cov", "sd", "prcomp", "eigen", "solve", "chol", "order", "pmax", "max.col")
  shadows <- shadows[vapply(shadows, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(shadows)) stop("Global objects shadow numerical functions/constants: ", paste(shadows, collapse = ", "))
  files <- list.files(project_root, recursive = TRUE, full.names = TRUE)
  find_one <- function(names) {
    hits <- files[basename(files) %in% names]
    if (length(hits) != 1L) stop("Expected one of ", paste(names, collapse = "/"), "; found ", length(hits), ".")
    normalizePath(hits, winslash = "/", mustWork = TRUE)
  }
  paths <- list(landsat = find_one("03_landsat_analysis.R"),
    bean = find_one(c("03_drybean_analysis.R.R", "03_drybean_analysis.R")),
    trn = find_one("sat.trn"), tst = find_one("sat.tst"),
    xlsx = find_one("Dry_Bean_Dataset.xlsx"))
  expected <- c("67b711404f884efec1a3f1d51340cbf5", "d92057b86d189070ed12cc3caff12621",
    "2c5ba2900da0183cab2c41fdb279fa5b", "02c995991fecc864e809b2c4c42cd983", "a024b18787414f9881d21faae54ce9c7")
  if (!identical(unname(tools::md5sum(unlist(paths))), expected)) stop("Original source or data hashes differ from the audited inputs.")
  review <- readRDS(review_path)
  if (!identical(review$protocol$version, "2026-10-03-realdata-joint-v1") ||
      length(review$records) != 100L || any(review$status$CompleteReplications != 50L))
    stop("The complete, jointly converged validation 07 run is required.")
  raw <- review$tables$raw
  if (nrow(raw) != 1500L || anyDuplicated(raw[, c("Dataset", "Replication", "Budget", "Method")]))
    stop("Incomplete or duplicated experiment records.")
  METHODS <- c("Random", "Entropy", "Margin", "Fisher", "AdaptiveRisk")
  LABELS <- c(Random = "Random", Entropy = "Entropy", Margin = "Margin",
    Fisher = "Fisher", AdaptiveRisk = "Adaptive risk")
  COLORS <- c(Random = "#6B7280", Entropy = "#0072B2", Margin = "#E69F00",
    Fisher = "#D55E00", AdaptiveRisk = "#009E73")
  SHAPES <- c(Random = 16, Entropy = 17, Margin = 15, Fisher = 18, AdaptiveRisk = 8)
  sym <- function(M) (M + t(M)) / 2
  inverse <- function(M) chol2inv(chol(sym(M)))
  relative <- function(a, b) max(abs(a - b) / (1 + abs(b)))
  out <- tempfile("realdata_figures_", tmpdir = project_root)
  if (!dir.create(out)) stop("Cannot create figure output directory.")
  logfile <- file.path(out, "validation_08_log.txt")
  sink(logfile, split = TRUE); on.exit(sink(), add = TRUE)
  cat("Cached-fit real-data figure reproduction; ", format(Sys.time()), "\nOutput: ", out, "\n", sep = "")
  cat("All plotted acquisition sets are read from validated, saved records.\n")
  cat("Representative replications minimize distance from the median paired error-count difference; ties use the lowest replication number.\n")
  cat("Classification-value colors use a 99th-percentile cap for display only.\n")
  cat("Performance error bars are plus/minus 1.96 Monte Carlo standard errors, conditional on the supplied datasets.\n")
  print(sessionInfo())

  # Load function definitions only; no original analysis loop is evaluated.
  definitions <- function(path) {
    e <- new.env(parent = globalenv())
    list2env(list(COV_FLOOR_REL = 1e-6, INFO_EIG_FLOOR_REL = 1e-8, N_PC = 5L), e)
    for (expr in parse(file = path, keep.source = FALSE)) {
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
      rhs <- expr[[3L]]
      if (is.symbol(expr[[2L]]) && is.call(rhs) && identical(rhs[[1L]], as.name("function"))) eval(expr, e)
    }
    e$safe_inverse <- function(M, rel_floor = 1e-8) inverse(M)
    e
  }
  cores <- list(Landsat = definitions(paths$landsat), DryBean = definitions(paths$bean))
  ls_train <- as.matrix(read.table(paths$trn, header = FALSE))
  ls_test <- as.matrix(read.table(paths$tst, header = FALSE))
  bean <- as.data.frame(readxl::read_excel(paths$xlsx)); bean$Class <- factor(bean$Class)
  bean_source_rows <- which(!duplicated(bean)); bean <- bean[bean_source_rows, , drop = FALSE]
  if (!identical(bean_source_rows, review$bean_source_rows) || nrow(bean) != 13543L) stop("Dry Bean cleaning differs from the audited run.")
  bean_X <- as.matrix(bean[, setdiff(names(bean), "Class"), drop = FALSE])

  data_for <- function(ds, record) {
    e <- cores[[ds]]; prep <- record$preprocessing
    if (ds == "Landsat") {
      X <- sweep(sweep(ls_train[, 17:20, drop = FALSE], 2L, prep$center, "-"), 2L, prep$scale, "/")
      y <- factor(as.character(ls_train[, 37L]), levels = record$pilot_fit$classes)
      if (relative(colMeans(ls_train[, 17:20, drop = FALSE]), prep$center) > 1e-12 ||
          relative(apply(ls_train[, 17:20, drop = FALSE], 2L, sd), prep$scale) > 1e-12)
        stop("Landsat preprocessing does not reproduce.")
      pcs <- prcomp(X, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
      list(X = X, y = y, pcs = pcs, n_test = nrow(ls_test))
    } else {
      train <- record$train_indices; test <- record$test_indices
      fresh <- e$fit_preprocess(bean_X[train, , drop = FALSE], 5L)
      if (relative(unlist(fresh), unlist(prep)) > 1e-10) stop("Dry Bean preprocessing does not reproduce.")
      X <- e$apply_preprocess(bean_X[train, , drop = FALSE], prep)
      y <- factor(bean$Class[train], levels = record$pilot_fit$classes)
      list(X = X, y = y, pcs = X[, 1:2, drop = FALSE], n_test = length(test))
    }
  }
  representative <- function(ds, budget, n_test) {
    z <- raw[raw$Dataset == ds & raw$Budget == budget, ]
    a <- z[z$Method == "AdaptiveRisk", ]; b <- z[z$Method == "Entropy", ]
    b <- b[match(a$Replication, b$Replication), ]
    difference <- round(n_test * b$Error) - round(n_test * a$Error)
    distance <- abs(difference - median(difference))
    chosen <- min(a$Replication[distance == min(distance)])
    list(replication = chosen, median_count = median(difference), difference = difference,
      replications = a$Replication, tied = a$Replication[distance == min(distance)])
  }
  base_theme <- function() ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#EDF0F3", linewidth = 0.25),
      plot.title = ggplot2::element_text(size = 11, face = "plain"),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8), plot.margin = ggplot2::margin(6, 6, 6, 6))
  save_plot <- function(plot, filename, width, height) {
    ggplot2::ggsave(file.path(out, paste0(filename, ".png")), plot, width = width, height = height,
      units = "in", dpi = 320, bg = "white")
    ggplot2::ggsave(file.path(out, paste0(filename, ".pdf")), plot, width = width, height = height,
      units = "in", device = grDevices::pdf, useDingbats = FALSE, bg = "white")
  }

  # Recompute all summary statistics from raw replication records.
  summaries <- list(); si <- 1L
  for (ds in c("Landsat", "DryBean")) for (budget in c(.1, .2, .3))
    for (outcome in c("Error", "BalancedError")) for (method in METHODS) {
      v <- raw[raw$Dataset == ds & raw$Budget == budget & raw$Method == method, outcome]
      if (length(v) != 50L) stop("A summary does not contain all 50 replications.")
      summaries[[si]] <- data.frame(Dataset = ds, Budget = budget, Outcome = outcome,
        Method = method, N = length(v), Mean = mean(v), MCSE = sd(v) / sqrt(length(v)))
      si <- si + 1L
    }
  summary <- do.call(rbind, summaries)
  key <- function(z) paste(z$Dataset, z$Budget, z$Outcome, z$Method, sep = "|")
  ref <- review$summary[match(key(summary), key(review$summary)), ]
  if (relative(as.matrix(summary[, c("Mean", "MCSE")]), as.matrix(ref[, c("Mean", "MCSE")])) > 1e-12)
    stop("Summary statistics do not reproduce.")
  write.csv(summary, file.path(out, "validation_08_summary.csv"), row.names = FALSE)
  write.csv(review$paired, file.path(out, "validation_08_paired.csv"), row.names = FALSE)
  diagnostics <- list()
  for (ds in c("Landsat", "DryBean")) {
    budget <- if (ds == "Landsat") .2 else .1
    n_test <- if (ds == "Landsat") 2000L else 4065L
    choice <- representative(ds, budget, n_test)
    r <- review$records[[paste0(ds, "_", choice$replication)]]
    e <- cores[[ds]]; dat <- data_for(ds, r); X <- dat$X; n <- nrow(X)
    if (dat$n_test != n_test) stop("Test sample size differs.")
    pilot <- rep(FALSE, n); pilot[r$pilot] <- TRUE; candidate <- which(!pilot)
    bi <- match(budget, review$protocol$budgets); design <- r$adaptive[[bi]]
    score <- e$compute_score_components(X, r$pilot_fit, e$build_param_index(length(r$pilot_fit$pi), ncol(X)))
    M <- e$build_information_from_design(score, pilot, candidate, design$a)
    inv <- inverse(M); objective <- sum(r$H * t(inv))
    if (relative(objective, design$objective) > 1e-8) stop("Saved relaxed design objective does not reproduce.")
    psi <- e$classification_value_from_G(score, sym(inv %*% r$H %*% inv), candidate = NULL)
    if (any(!is.finite(psi)) || min(psi) < -1e-8 * max(1, max(abs(psi)))) stop("Invalid classification value.")
    psi <- pmax(psi, 0)
    cap <- unname(quantile(psi, .99, type = 7))
    if (!is.finite(cap) || cap <= 0) stop("Classification value is numerically zero.")
    entropy <- e$score_entropy(score$tau)
    selected_a <- r$final[[paste0("AdaptiveRisk_", budget)]]$selected
    selected_e <- r$final[[paste0("Entropy_", budget)]]$selected
    B <- ceiling(budget * n) - length(r$pilot)
    if (length(selected_a) != B || length(selected_e) != B ||
        length(intersect(selected_a, r$pilot)) || length(intersect(selected_e, r$pilot)))
      stop("Acquisition sets do not satisfy the budget.")
    entropy_expected <- candidate[order(entropy[candidate], decreasing = TRUE)[seq_len(B)]]
    adaptive_expected <- candidate[order(design$a, decreasing = TRUE)[seq_len(B)]]
    if (!identical(selected_a, adaptive_expected) || !identical(selected_e, entropy_expected))
      stop("Displayed acquisition sets do not reproduce.")
    plot_df <- data.frame(Row = seq_len(n), PC1 = dat$pcs[, 1L], PC2 = dat$pcs[, 2L],
      Class = dat$y, Entropy = entropy, ClassificationValue = psi,
      RelativeValue = pmin(psi / cap, 1), Pilot = pilot, Candidate = !pilot,
      AdaptiveSelected = seq_len(n) %in% selected_a, EntropySelected = seq_len(n) %in% selected_e)
    prefix <- if (ds == "Landsat") "Landsat" else "DryBean_realdata"
    write.csv(plot_df, file.path(out, paste0(prefix, "_geometry_plot_data.csv")), row.names = FALSE)
    overlap <- length(intersect(selected_a, selected_e))
    diagnostics[[ds]] <- data.frame(Dataset = ds, Replication = choice$replication, Budget = budget,
      MedianEntropyMinusAdaptiveErrorCount = choice$median_count,
      TiedReplications = paste(sort(choice$tied), collapse = ";"), AdditionalBudget = B,
      SharedAcquiredPoints = overlap, DifferentPointsPerRule = B - overlap,
      EntropyPsiSpearman = cor(entropy[candidate], psi[candidate], method = "spearman"),
      PsiColourCap99 = cap, RelaxedObjectiveDiscrepancy = abs(objective - design$objective),
      AdaptiveTestError = r$final[[paste0("AdaptiveRisk_", budget)]]$error,
      EntropyTestError = r$final[[paste0("Entropy_", budget)]]$error)
    cat("\n", ds, ": representative replication ", choice$replication,
      "; median paired error count ", choice$median_count, "; tied replications ",
      paste(sort(choice$tied), collapse = ","), "; acquired-set overlap ", overlap, "/", B, "\n", sep = "")
    class_colors <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#332288")[seq_along(levels(dat$y))]
    names(class_colors) <- levels(dat$y)
    class_labels <- if (ds == "Landsat") c("Red soil", "Cotton", "Grey soil", "Damp grey",
      "Stubble", "Very damp grey") else levels(dat$y)
    p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(PC1, PC2, colour = Class)) +
      ggplot2::geom_point(size = .55, alpha = .55) +
      ggplot2::scale_colour_manual(values = class_colors, labels = class_labels, name = "Class") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3, override.aes = list(size = 2, alpha = 1))) +
      ggplot2::labs(title = "(a) Feature geometry", x = "PC1", y = "PC2") + base_theme()
    p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(PC1, PC2, colour = Entropy)) +
      ggplot2::geom_point(size = .55, alpha = .65) +
      ggplot2::scale_colour_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
        limits = c(0, log(length(r$pilot_fit$pi))), name = "Entropy") +
      ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(3.8, "cm"))) +
      ggplot2::labs(title = "(b) Posterior uncertainty", x = "PC1", y = "PC2") + base_theme()
    p3 <- ggplot2::ggplot(plot_df, ggplot2::aes(PC1, PC2, colour = RelativeValue)) +
      ggplot2::geom_point(size = .55, alpha = .65) +
      ggplot2::scale_colour_gradientn(colours = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
        limits = c(0, 1), name = "Relative value (99% cap)") +
      ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(3.8, "cm"))) +
      ggplot2::labs(title = "(c) Classification value", x = "PC1", y = "PC2") + base_theme()
    selection <- rbind(data.frame(PC1 = plot_df$PC1, PC2 = plot_df$PC2,
      Selected = plot_df$AdaptiveSelected, Rule = "Adaptive risk"),
      data.frame(PC1 = plot_df$PC1, PC2 = plot_df$PC2,
      Selected = plot_df$EntropySelected, Rule = "Entropy"))
    selection$Rule <- factor(selection$Rule, levels = c("Adaptive risk", "Entropy"))
    p4 <- ggplot2::ggplot(selection, ggplot2::aes(PC1, PC2)) +
      ggplot2::geom_point(colour = "#CBD0D6", size = .45, alpha = .35) +
      ggplot2::geom_point(data = selection[selection$Selected, ], colour = "#111827", size = .7, alpha = .8) +
      ggplot2::facet_wrap(~ Rule, nrow = 1) +
      ggplot2::labs(title = paste0("(d) Additional labels: ", round(100 * budget), "% total budget"),
        x = "PC1", y = "PC2") + base_theme() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
        strip.text = ggplot2::element_text(size = 10))
    geometry <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) / p4 +
      patchwork::plot_layout(heights = c(1.05, 1)) +
      patchwork::plot_annotation(title = if (ds == "Landsat") "Statlog Landsat Satellite" else "Dry Bean",
        subtitle = paste0("Representative replication ", choice$replication,
          "; full fitted model used for scores and decisions"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = .5),
          plot.subtitle = ggplot2::element_text(size = 10, hjust = .5)))
    save_plot(geometry, paste0(prefix, "_geometry_acquisition"), 11.6, 7.5)
    performance <- summary[summary$Dataset == ds & summary$Outcome == "Error", ]
    performance$Method <- factor(performance$Method, levels = METHODS)
    performance$Lower <- performance$Mean - 1.96 * performance$MCSE
    performance$Upper <- performance$Mean + 1.96 * performance$MCSE
    pp <- ggplot2::ggplot(performance, ggplot2::aes(Budget, Mean, colour = Method, shape = Method, group = Method)) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = .004, linewidth = .45) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_colour_manual(values = COLORS, labels = LABELS, name = NULL) +
      ggplot2::scale_shape_manual(values = SHAPES, labels = LABELS, name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(.1, .2, .3), labels = c("10%", "20%", "30%"),
        limits = c(.085, .315)) +
      ggplot2::labs(x = "Total labeling budget", y = "Mean test classification error") + base_theme() +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
    perf_name <- if (ds == "Landsat") "Landsat_performance" else "DryBean_realdata_error"
    save_plot(pp, perf_name, 7.7, 4.9)
  }
  diagnostics <- do.call(rbind, diagnostics)
  write.csv(diagnostics, file.path(out, "validation_08_geometry_diagnostics.csv"), row.names = FALSE)
  cat("\nFigure diagnostics:\n"); print(diagnostics, row.names = FALSE)
  cat("\nReturn these seven files, or zip them together:\n")
  cat(file.path(out, c("Landsat_geometry_acquisition.png", "Landsat_performance.png",
    "DryBean_realdata_geometry_acquisition.png", "DryBean_realdata_error.png",
    "validation_08_geometry_diagnostics.csv", "validation_08_summary.csv", "validation_08_log.txt")), sep = "\n")
  cat("\nKeep both geometry_plot_data.csv files and the vector PDFs locally.\n")
  cat("No fitted model, acquisition set, or test-error record has been changed.\n")
  invisible(out)
}
