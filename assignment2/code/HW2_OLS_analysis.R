# Philadelphia block groups: reproducible Part 1 (OLS and spatial diagnostics).
# Run from code/ or the project root. The same script can be sourced by knitr.
# Reconstructed from the group's shared source data and Part 1 handoff.
run_hw2_part1 <- function() {
  root_candidates <- c(".", "..")
  root <- root_candidates[file.exists(file.path(root_candidates, "data",
    "map_RegressionData", "RegressionData.shp"))][1]
  if (is.na(root)) stop("Run this script from the project root or code/.")
  root <- normalizePath(root)
  input <- file.path(root, "data")
  out <- file.path(root, "output")
  processed <- file.path(input, "processed_ols")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed, recursive = TRUE, showWarnings = FALSE)

  shp <- sf::st_read(file.path(input, "map_RegressionData", "RegressionData.shp"),
                     quiet = TRUE, stringsAsFactors = FALSE)
  dat <- sf::st_drop_geometry(shp)
  csv <- read.csv(file.path(input, "RegressionData.csv"))
  # The CSV contains rounded AREAKEY strings; align on the unique POLY_ID.
  stopifnot(nrow(dat) == 1720L, nrow(csv) == nrow(dat),
            !anyDuplicated(dat$POLY_ID), !anyDuplicated(csv$POLY_ID),
            setequal(csv$POLY_ID, dat$POLY_ID), all(dat$MEDHVAL > 0),
            all(dat$NBelPov100 >= 0))
  csv <- csv[match(dat$POLY_ID, csv$POLY_ID), ]
  common <- c("PCTBACHMOR", "MEDHHINC", "PCTVACANT", "PCTSINGLES", "MEDHVAL")
  stopifnot(all(vapply(common, function(x) all(csv[[x]] == dat[[x]]), logical(1))),
            all(csv$NBELPOV100 == dat$NBelPov100))
  shp$LNMEDHVAL <- log1p(dat$MEDHVAL)
  shp$LNNBELPOV <- log1p(dat$NBelPov100)
  log_differences <- c(max(abs(shp$LNMEDHVAL - dat$LNMEDHVAL)),
                       max(abs(shp$LNNBELPOV - dat$LNNBELPOV)))
  dat <- sf::st_drop_geometry(shp)
  formula <- LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV
  model_vars <- all.vars(formula)
  stopifnot(all(is.finite(as.matrix(dat[model_vars]))))
  # The original course GAL carries the Queen definition. Preserve its ordered
  # neighbor indices and verify the shared GAL used in Parts 2 and 3.
  ids <- as.character(dat$POLY_ID)
  queen <- spdep::read.gal(file.path(input, "map_RegressionData", "RegressionData.gal"),
                           region.id = ids)
  shared <- file.path(input, "RegressionData_queen.gal")
  if (!file.exists(shared)) spdep::write.nb.gal(queen, shared)
  nb <- spdep::read.gal(shared, region.id = ids)
  stopifnot(identical(attr(nb, "region.id"), ids),
            identical(attr(queen, "region.id"), ids),
            all(vapply(seq_along(nb), function(i) {
              identical(queen[[i]], nb[[i]])
            }, logical(1))))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  weight_description <- "Row-standardized first-order Queen"
  cards <- spdep::card(nb)
  write_table <- function(d, name) write.csv(d, file.path(out, name), row.names = FALSE)
  write_table(data.frame(Item = c("Observations", "Missing values in model variables",
    "Maximum difference: stored LNMEDHVAL vs log1p(MEDHVAL)",
    "Maximum difference: stored LNNBELPOV vs log1p(NBelPov100)",
    "Queen neighbor links (directed)", "Mean Queen neighbors",
    "Block groups with zero neighbors", "Queen subgraphs"),
    Value = c(nrow(dat), sum(is.na(dat[model_vars])), log_differences, sum(cards),
              mean(cards), sum(cards == 0L), spdep::n.comp.nb(nb)$nc)),
    "table_01_data_and_weights_validation.csv")

  permutations <- 999L
  seed <- 671L
  global_moran <- function(values, name, filename) {
    set.seed(seed)
    mc <- spdep::moran.mc(values, lw, nsim = permutations, alternative = "greater")
    write_table(data.frame(Variable = name, Moran_I = unname(mc$statistic),
      Permutations = permutations, Pseudo_p_greater = mc$p.value,
      Weight_style = weight_description), filename)
    mc
  }
  response_mc <- global_moran(dat$LNMEDHVAL, "LNMEDHVAL", "table_02_global_moran_lnmedhval.csv")
  # Serial defaults plus a fixed seed reproduce the supplied permutation stream.
  previous_cores <- spdep::get.coresOption()
  spdep::set.coresOption(NULL)
  on.exit(spdep::set.coresOption(previous_cores), add = TRUE)
  lisa <- spdep::localmoran_perm(dat$LNMEDHVAL, lw, nsim = permutations,
                                 alternative = "two.sided", iseed = seed,
                                 no_repeat_in_row = TRUE)
  shp$LISA_I <- unname(lisa[, "Ii"])
  shp$LISA_P <- unname(lisa[, "Pr(z != E(Ii)) Sim"])
  z <- as.numeric(scale(dat$LNMEDHVAL))
  wz <- spdep::lag.listw(lw, z)
  clusters <- c("Not significant", "High-High", "Low-Low", "High-Low", "Low-High")
  cl <- ifelse(z >= 0, ifelse(wz >= 0, "High-High", "High-Low"),
                       ifelse(wz < 0, "Low-Low", "Low-High"))
  cl[shp$LISA_P > 0.05] <- "Not significant"
  shp$LISA_CL <- cl
  counts <- as.integer(table(factor(cl, levels = clusters)))
  write_table(data.frame(LISA_cluster = clusters, Block_groups = counts,
                         Share = counts / nrow(dat)), "table_03_lisa_cluster_counts.csv")

  ols <- lm(formula, data = dat)
  summary_ols <- summary(ols)
  coef_table <- function(model) {
    values <- coef(summary(model))
    data.frame(Term = rownames(values), unname(values), check.names = FALSE) |>
      setNames(c("Term", colnames(values)))
  }
  write_table(coef_table(ols), "table_04_ols_coefficients.csv")
  f <- summary_ols$fstatistic
  write_table(data.frame(Metric = c("N", "R-squared", "Adjusted R-squared",
    "Residual standard error", "F statistic", "F numerator df", "F denominator df",
    "F p-value", "Log likelihood", "AIC", "Schwarz criterion (BIC)"),
    Value = c(nobs(ols), summary_ols$r.squared, summary_ols$adj.r.squared,
      summary_ols$sigma, f, pf(f[1], f[2], f[3], lower.tail = FALSE),
      as.numeric(logLik(ols)), AIC(ols), BIC(ols))), "table_05_ols_fit_statistics.csv")
  shp$OLS_PRED <- fitted(ols)
  shp$OLS_RESID <- resid(ols)
  shp$OLS_RSTD <- rstandard(ols)
  shp$WT_RESID <- spdep::lag.listw(lw, shp$OLS_RESID)
  bp <- lmtest::bptest(ols, studentize = FALSE)
  kb <- lmtest::bptest(ols, studentize = TRUE)
  white <- lm(I(shp$OLS_RESID^2) ~ (PCTVACANT + PCTSINGLES + PCTBACHMOR +
               LNNBELPOV)^2 + I(PCTVACANT^2) + I(PCTSINGLES^2) +
               I(PCTBACHMOR^2) + I(LNNBELPOV^2), data = dat)
  white_stat <- nrow(dat) * summary(white)$r.squared
  centered <- shp$OLS_RESID - mean(shp$OLS_RESID)
  moments <- vapply(2:4, function(k) mean(centered^k), numeric(1))
  skewness <- moments[2] / moments[1]^1.5
  kurtosis <- moments[3] / moments[1]^2
  jb <- nrow(dat) / 6 * (skewness^2 + (kurtosis - 3)^2 / 4)
  stats <- unname(c(bp$statistic, kb$statistic, white_stat, jb))
  dfs <- unname(c(bp$parameter, kb$parameter, white$rank - 1L, 2L))
  write_table(data.frame(Test = c("Breusch-Pagan (classical)",
    "Koenker-Bassett (studentized BP)", "White heteroscedasticity test",
    "Jarque-Bera normality test"), Statistic = stats, Df = dfs,
    P_value = pchisq(stats, dfs, lower.tail = FALSE), Null_hypothesis =
    c(rep("Constant error variance", 3), "Normally distributed errors")),
    "table_06_ols_diagnostic_tests.csv")
  descriptive <- lm(OLS_RESID ~ WT_RESID, sf::st_drop_geometry(shp))
  write_table(coef_table(descriptive), "table_07_ols_residual_on_weighted_residual.csv")
  residual_mc <- global_moran(shp$OLS_RESID, "OLS_RESID", "table_08_ols_residual_moran.csv")

  plot_file <- function(name, draw) {
    png(file.path(out, name), width = 1200, height = 900, res = 140)
    on.exit(dev.off())
    par(mar = c(4.4, 4.7, 3.3, 1.1))
    draw()
  }
  point <- "#1F668883"
  accent <- "#AD4025"
  moran_scatter <- function(values, mc, title, label, name) plot_file(name, function() {
    zz <- as.numeric(scale(values))
    wzz <- spdep::lag.listw(lw, zz)
    plot(zz, wzz, pch = 16, col = point, cex = 0.55,
         xlab = label, ylab = "Queen-weighted neighboring value", main = title)
    abline(h = 0, v = 0, lty = 3, col = "grey40")
    abline(lm(wzz ~ zz), lwd = 2, col = accent)
    legend("topleft", bty = "n", legend = sprintf("Moran's I = %.3f", mc$statistic))
  })
  permutation_plot <- function(mc, title, name) plot_file(name, function() {
    simulated <- head(mc$res, permutations)
    hist(simulated, breaks = 35, col = "#7B9BAA", border = "white",
         xlim = range(c(simulated, mc$statistic)), xlab = "Moran's I", main = title)
    abline(v = mc$statistic, col = accent, lwd = 2)
    legend("topright", bty = "n", legend = sprintf("Observed I = %.3f\nUpper-tail p = %.3f",
           mc$statistic, mc$p.value))
  })
  moran_scatter(dat$LNMEDHVAL, response_mc, "Spatial association in log median home value",
    "Standardized log median home value", "fig_01_moran_lnmedhval_scatter.png")
  permutation_plot(response_mc, "Randomized locations of log median home value",
                   "fig_02_moran_lnmedhval_permutation.png")
  map_file <- function(groups, palette, title, filename) plot_file(filename, function() {
    plot(sf::st_geometry(shp), col = unname(palette[groups]), border = "#777777",
         lwd = 0.13, main = title)
    legend("bottomright", legend = names(palette), fill = unname(palette),
           bg = "white", box.col = "grey85", cex = 0.8)
  })
  map_file(ifelse(shp$LISA_P <= 0.05, "p <= 0.05", "p > 0.05"),
           c("p <= 0.05" = "#3C7E78", "p > 0.05" = "#E5E5E5"),
           "Local Moran permutation significance", "fig_03_lisa_significance_map.png")
  map_file(cl, setNames(c("#E5E5E5", "#B2182B", "#2166AC", "#EF8A62", "#67A9CF"),
                        clusters), "Local spatial clusters and outliers",
           "fig_04_lisa_cluster_map.png")
  plot_file("fig_05_ols_residuals_vs_fitted.png", function() {
    plot(shp$OLS_PRED, shp$OLS_RSTD, pch = 16, col = point, cex = 0.55,
         xlab = "Fitted log median home value", ylab = "Standardized OLS residual",
         main = "OLS residuals versus fitted values")
    abline(h = 0, lty = 2)
  })
  plot_file("fig_06_ols_residual_histogram.png", function() {
    hist(shp$OLS_RSTD, breaks = 40, probability = TRUE, col = "#7B9BAA",
         border = "white", main = "Standardized OLS residuals", xlab = "Standardized residual")
    curve(dnorm(x), add = TRUE, col = accent, lwd = 2)
  })
  plot_file("fig_07_ols_residual_vs_weighted_residual.png", function() {
    slope <- coef(summary(descriptive))["WT_RESID", ]
    plot(shp$WT_RESID, shp$OLS_RESID, pch = 16, col = point, cex = 0.55,
         xlab = "Queen-weighted neighboring residual", ylab = "OLS residual",
         main = "Local residual versus neighboring residuals")
    abline(descriptive, col = accent, lwd = 2)
    legend("bottomright", bty = "n", cex = 0.9,
           legend = sprintf("Slope b = %.3f; SE = %.3f\nt(1718) = %.2f; p < 0.001",
                            slope[1], slope[2], slope[3]))
  })
  moran_scatter(shp$OLS_RESID, residual_mc, "Spatial association in OLS residuals",
                "Standardized OLS residual", "fig_08_ols_residual_moran_scatter.png")
  permutation_plot(residual_mc, "Randomized locations of OLS residuals",
                   "fig_09_ols_residual_moran_permutation.png")

  # Keep column order compatible with the supplied handoff, and preserve original
  # polygon coordinates. This topological step needs no assignment of a CRS.
  fields <- c(names(dat), "OLS_PRED", "OLS_RESID", "OLS_RSTD", "WT_RESID",
              "LISA_I", "LISA_P", "LISA_CL")
  shp <- shp[, fields]
  sf::st_write(shp, file.path(processed, "RegressionFinal_OLS.shp"),
               append = FALSE, quiet = TRUE)
  write.csv(sf::st_drop_geometry(shp), file.path(processed, "RegressionFinal_OLS_attributes.csv"),
            row.names = FALSE)
  writeLines(capture.output(print(summary_ols)), file.path(out, "ols_model_summary.txt"))
  writeLines(c("Homework 2 Part 1: OLS and spatial diagnostics",
    paste("Model:", deparse(formula)), paste("Block groups:", nrow(dat)),
    sprintf("Queen links (directed): %d", sum(cards)),
    sprintf("Log home value Moran's I = %.6f (positive-tail p = %.3f)",
            response_mc$statistic, response_mc$p.value),
    sprintf("OLS R-squared = %.6f, adjusted = %.6f, AIC = %.3f, BIC = %.3f",
            summary_ols$r.squared, summary_ols$adj.r.squared, AIC(ols), BIC(ols)),
    sprintf("OLS residual Moran's I = %.6f (positive-tail p = %.3f)",
            residual_mc$statistic, residual_mc$p.value),
    sprintf("Slope in descriptive residual-on-neighbor regression = %.6f", coef(descriptive)[2]),
    sprintf("Permutation count = %d; reproducibility seed = %d", permutations, seed),
    "Local Moran p-values are two-sided and unadjusted for multiple comparisons.",
    "All four OLS slopes are significant at 0.05 with conventional OLS standard errors.",
    "Heteroscedasticity and positive residual spatial dependence warrant caution about inference.",
    "The Jarque-Bera diagnostic also rejects normal errors.",
    "Use the shared Queen GAL with the processed polygons for the next spatial analyses."),
    file.path(out, "analysis_handoff_summary.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
  paths <- c(file.path("output", list.files(out, "^(fig_0[1-9]|table_0[1-8])_")),
             file.path("output", c("ols_model_summary.txt", "analysis_handoff_summary.txt",
                                   "session_info.txt")), "data/RegressionData_queen.gal",
             file.path("data", "processed_ols", list.files(processed, "^RegressionFinal_OLS")))
  write_table(data.frame(File = paths), "output_manifest.csv")
  invisible(TRUE)
}
run_hw2_part1()
