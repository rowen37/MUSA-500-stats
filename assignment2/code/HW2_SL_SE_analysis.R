# ==========================================================================
# Homework 2 - Part 2: Spatial Lag (SL) and Spatial Error (SE) regression
#
# Authors: Qiwen Bian, Manching Ho, Fujia Zheng
#
# Model (identical to Part 1):
#   LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV
#
# Weights (identical to Part 1): row-standardised first-order Queen weights
# read from data/RegressionData_queen.gal
#
# Reproduce all Part 2 outputs:
#   cd code
#   Rscript HW2_SL_SE_analysis.R
# (the script also works from the project root)
# ==========================================================================

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(spatialreg)
  library(lmtest)
  library(ggplot2)
})

## ---- 0. Project root, model specification, plotting defaults -------------

find_project_root <- function() {
  for (p in c(".", "..")) {
    gal <- file.path(p, "data", "RegressionData_queen.gal")
    if (file.exists(gal) && dir.exists(file.path(p, "output"))) {
      return(normalizePath(p))
    }
  }
  stop("Cannot locate the project root (data/RegressionData_queen.gal).")
}

project_root <- find_project_root()
if (normalizePath(".") != project_root) setwd(project_root)

model_formula <-
  LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV

permutations <- 999
perm_seed <- 671

# Jarque-Bera's two-moment test; computed directly to avoid an optional package.
jarque.bera.test <- function(x) {
  z <- x - mean(x)
  m2 <- mean(z^2)
  jb <- length(z) * ((mean(z^3) / m2^1.5)^2 / 6 +
                        (mean(z^4) / m2^2 - 3)^2 / 24)
  list(statistic = jb, parameter = 2, p.value = pchisq(jb, 2, lower.tail = FALSE))
}

# spdep stores the observed statistic last. Count it once in each tail so
# finite Monte Carlo ranks have a positive lower bound, including at the extremes.
resid_moran <- function(x) {
  set.seed(perm_seed)
  m <- spdep::moran.mc(x, queenlist, nsim = permutations, alternative = "two.sided")
  simulated <- head(m$res, -1L)
  upper <- (1 + sum(simulated >= m$statistic)) / (permutations + 1)
  lower <- (1 + sum(simulated <= m$statistic)) / (permutations + 1)
  m$p.value <- min(1, 2 * min(upper, lower))
  m$p.greater <- upper
  m
}

write_table <- function(x, name) {
  write.csv(x, file.path("output", name), row.names = FALSE)
}

## ---- 1. Shared data and weights ------------------------------------------

# The Part 1 handoff shapefile carries the OLS fields, so Part 2 can be run
# on exactly the same block groups and attributes as Part 1.
shapefile_path <- file.path(
  "data", "processed_ols", "RegressionFinal_OLS.shp"
)
shp <- st_read(shapefile_path, quiet = TRUE)
if (nrow(shp) != 1720) {
  stop("Expected 1,720 block groups, found ", nrow(shp), ".")
}

needed <- c("POLY_ID", all.vars(model_formula))
if (length(setdiff(needed, names(shp))) > 0) {
  stop("Missing fields: ", paste(setdiff(needed, names(shp)), collapse = ", "))
}

# The supplied .prj is an undefined engineering CRS, which breaks sf's
# coordinate handling for plotting and writing. Queen contiguity and the
# spatial models do not use distance, so the undefined CRS is dropped.
shp <- st_set_crs(shp, NA)

queen <- read.gal(
  file.path("data", "RegressionData_queen.gal"),
  region.id = as.character(shp$POLY_ID)
)
queenlist <- nb2listw(queen, style = "W", zero.policy = TRUE)

model_data <- st_drop_geometry(shp)
cat("Block groups:", nrow(model_data), "\n")
cat("Neighbour cardinality:", min(card(queen)), "-", max(card(queen)),
    "| zero-neighbour units:", sum(card(queen) == 0), "\n")
cat("Weight style:", queenlist$style, "\n\n")

## ---- 2. OLS benchmark (needed for the Lagrange Multiplier and LR tests) ---

ols <- lm(model_formula, data = model_data)
ols_ll <- as.numeric(logLik(ols))
ols_aic <- AIC(ols)

## ---- 3. Lagrange Multiplier diagnostics: which spatial model? -----------

# RStests gives GeoDa's LM-lag / LM-error and their robust versions.
lm_rs <- lm.RStests(ols, queenlist, test = "all", zero.policy = TRUE)

lm_table <- data.frame(
  Test = c(
    "LM error (RSerr)",
    "LM lag (RSlag)",
    "Robust LM error (adjRSerr)",
    "Robust LM lag (adjRSlag)",
    "LM SARMA (error + lag)"
  ),
  Statistic = vapply(lm_rs, function(x) unname(x$statistic), numeric(1)),
  Df = vapply(lm_rs, function(x) unname(x$parameter), numeric(1)),
  P_value = vapply(lm_rs, function(x) x$p.value, numeric(1)),
  Null_hypothesis = c(
    "OLS residuals are spatially independent (error dependence)",
    "OLS residuals are spatially independent (lag dependence)",
    "No error dependence beyond lag dependence",
    "No lag dependence beyond error dependence",
    "No spatial dependence of either form"
  )
)
write_table(lm_table, "table_09_lm_diagnostics_ols.csv")
print(lm_table)

## ---- 4. Spatial Lag model ------------------------------------------------

sl <- lagsarlm(
  model_formula, data = model_data, listw = queenlist, zero.policy = TRUE
)
cat("\n===== Spatial Lag model =====\n")
print(summary(sl))
writeLines(
  capture.output(summary(sl)),
  file.path("output", "sl_model_summary.txt")
)

sl_coef <- as.data.frame(summary(sl)$Coef)
sl_coef <- cbind(Term = rownames(sl_coef), sl_coef)
rownames(sl_coef) <- NULL
write_table(sl_coef, "table_10_sl_coefficients.csv")

sl_rho <- sl$rho
sl_rho_se <- summary(sl)$rho.se
sl_rho_z <- sl_rho / sl_rho_se
sl_rho_p <- 2 * pnorm(-abs(sl_rho_z))

sl_lr <- LR.Sarlm(sl, ols)
sl_fit <- data.frame(
  Metric = c(
    "N",
    "Spatial autoregressive parameter (rho)",
    "Rho standard error",
    "Rho z value",
    "Rho p value (Wald test)",
    "Log likelihood",
    "AIC",
    "LR test vs OLS: statistic",
    "LR test vs OLS: df",
    "LR test vs OLS: p value",
    "OLS log likelihood (for reference)",
    "OLS AIC (for reference)"
  ),
  Value = c(
    nrow(model_data), sl_rho, sl_rho_se, sl_rho_z, sl_rho_p,
    as.numeric(logLik(sl)), AIC(sl),
    unname(sl_lr$statistic), unname(sl_lr$parameter), sl_lr$p.value,
    ols_ll, ols_aic
  )
)
write_table(sl_fit, "table_11_sl_fit_statistics.csv")
print(sl_fit)

sl_bp <- bptest.Sarlm(sl, studentize = FALSE)
sl_kb <- bptest.Sarlm(sl, studentize = TRUE)
sl_jb <- jarque.bera.test(residuals(sl))

sl_moran <- resid_moran(residuals(sl))
sl_moran_i <- unname(sl_moran$statistic)

sl_diag <- data.frame(
  Test = c(
    "Breusch-Pagan (classical)",
    "Koenker-Bassett (studentized BP)",
    "Jarque-Bera normality test",
    "Moran's I of SL residuals (999 permutations)"
  ),
  Statistic = c(unname(sl_bp$statistic), unname(sl_kb$statistic),
                unname(sl_jb$statistic), sl_moran_i),
  Df = c(unname(sl_bp$parameter), unname(sl_kb$parameter),
         unname(sl_jb$parameter), NA),
  P_value = c(sl_bp$p.value, sl_kb$p.value, sl_jb$p.value,
              sl_moran$p.value),
  Null_hypothesis = c(
    "Constant error variance",
    "Constant error variance",
    "Normally distributed errors",
    "Spatial independence of residuals"
  )
)
write_table(sl_diag, "table_12_sl_diagnostics.csv")
print(sl_diag)

## ---- 5. Spatial Error model ----------------------------------------------

se <- errorsarlm(
  model_formula, data = model_data, listw = queenlist, zero.policy = TRUE
)
cat("\n===== Spatial Error model =====\n")
print(summary(se))
writeLines(
  capture.output(summary(se)),
  file.path("output", "se_model_summary.txt")
)

se_coef <- as.data.frame(summary(se)$Coef)
se_coef <- cbind(Term = rownames(se_coef), se_coef)
rownames(se_coef) <- NULL
write_table(se_coef, "table_13_se_coefficients.csv")

se_lambda <- se$lambda
se_lambda_se <- summary(se)$lambda.se
se_lambda_z <- se_lambda / se_lambda_se
se_lambda_p <- 2 * pnorm(-abs(se_lambda_z))

se_lr <- LR.Sarlm(se, ols)
se_fit <- data.frame(
  Metric = c(
    "N",
    "Spatial autoregressive parameter (lambda)",
    "Lambda standard error",
    "Lambda z value",
    "Lambda p value (Wald test)",
    "Log likelihood",
    "AIC",
    "LR test vs OLS: statistic",
    "LR test vs OLS: df",
    "LR test vs OLS: p value",
    "OLS log likelihood (for reference)",
    "OLS AIC (for reference)"
  ),
  Value = c(
    nrow(model_data), se_lambda, se_lambda_se, se_lambda_z, se_lambda_p,
    as.numeric(logLik(se)), AIC(se),
    unname(se_lr$statistic), unname(se_lr$parameter), se_lr$p.value,
    ols_ll, ols_aic
  )
)
write_table(se_fit, "table_14_se_fit_statistics.csv")
print(se_fit)

se_bp <- bptest.Sarlm(se, studentize = FALSE)
se_kb <- bptest.Sarlm(se, studentize = TRUE)
se_jb <- jarque.bera.test(residuals(se))

se_moran <- resid_moran(residuals(se))
se_moran_i <- unname(se_moran$statistic)

se_diag <- data.frame(
  Test = c(
    "Breusch-Pagan (classical)",
    "Koenker-Bassett (studentized BP)",
    "Jarque-Bera normality test",
    "Moran's I of SE residuals (999 permutations)"
  ),
  Statistic = c(unname(se_bp$statistic), unname(se_kb$statistic),
                unname(se_jb$statistic), se_moran_i),
  Df = c(unname(se_bp$parameter), unname(se_kb$parameter),
         unname(se_jb$parameter), NA),
  P_value = c(se_bp$p.value, se_kb$p.value, se_jb$p.value,
              se_moran$p.value),
  Null_hypothesis = c(
    "Constant error variance",
    "Constant error variance",
    "Normally distributed errors",
    "Spatial independence of residuals"
  )
)
write_table(se_diag, "table_15_se_diagnostics.csv")
print(se_diag)

## ---- 6. OLS / SL / SE comparison -----------------------------------------

ols_moran <- resid_moran(residuals(ols))
ols_moran_gt <- list(p.value = ols_moran$p.greater)
sl_moran_gt <- list(p.value = sl_moran$p.greater)
se_moran_gt <- list(p.value = se_moran$p.greater)

comparison <- data.frame(
  Model = c("OLS", "Spatial Lag", "Spatial Error"),
  Log_likelihood = c(ols_ll, as.numeric(logLik(sl)), as.numeric(logLik(se))),
  AIC = c(ols_aic, AIC(sl), AIC(se)),
  BIC = c(BIC(ols), BIC(sl), BIC(se)),
  Spatial_parameter = c(NA, sl_rho, se_lambda),
  Residual_Moran_I = c(unname(ols_moran$statistic), sl_moran_i, se_moran_i),
  Moran_p_two_sided = c(ols_moran$p.value, sl_moran$p.value, se_moran$p.value),
  Moran_p_greater = c(ols_moran_gt$p.value, sl_moran_gt$p.value,
                      se_moran_gt$p.value),
  LR_vs_OLS_statistic = c(NA, unname(sl_lr$statistic), unname(se_lr$statistic)),
  LR_vs_OLS_p = c(NA, sl_lr$p.value, se_lr$p.value)
)
write_table(comparison, "table_16_model_comparison_ols_sl_se.csv")
print(comparison)

## ---- 7. Figures ----------------------------------------------------------

map_plot <- function(data, fill_var, title, subtitle) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "grey85", linewidth = 0.05) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, name = "Standardized\nresidual"
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 12) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )
}

moran_plot <- function(y, listw, title, moran_i) {
  z <- as.numeric(scale(y))
  wz <- lag(listw, z)
  ggplot(data.frame(z = z, wz = wz), aes(x = z, y = wz)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_point(alpha = 0.45, size = 0.8, colour = "#2c7fb8") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#b2182b", linewidth = 0.7) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             family = "serif", size = 3.6,
             label = sprintf("Moran's I = %.3f", moran_i)) +
    labs(
      title = title,
      subtitle = "x = standardized residual, y = spatial lag of the standardized residual",
      x = "Standardized residual", y = "Spatially lagged value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5)
    )
}

perm_plot <- function(mc, title, observed) {
  res <- data.frame(I = head(mc$res, -1L))
  ggplot(res, aes(x = I)) +
    geom_histogram(bins = 40, fill = "grey65", colour = "white") +
    geom_vline(xintercept = observed, colour = "#b2182b", linewidth = 0.9) +
    labs(
      title = title,
      subtitle = sprintf(
        "999 permutations; observed I = %.3f, two-sided pseudo p = %.3f",
        observed, mc$p.value
      ),
      x = "Moran's I from random permutations", y = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5)
    )
}

sl_sd <- sd(residuals(sl))
se_sd <- sd(residuals(se))
shp$SL_RSTD <- as.numeric(residuals(sl)) / sl_sd
shp$SE_RSTD <- as.numeric(residuals(se)) / se_sd
shp$SL_PRED <- as.numeric(fitted(sl))
shp$SE_PRED <- as.numeric(fitted(se))

save_plot <- function(p, name, width = 9, height = 7) {
  ggsave(file.path("output", name), p, width = width, height = height, dpi = 200)
}

save_plot(
  map_plot(shp, "SL_RSTD", "Standardized Spatial Lag Residuals",
           "Philadelphia block groups, n = 1,720"),
  "fig_10_sl_residual_map.png", width = 9, height = 8
)
save_plot(
  moran_plot(residuals(sl), queenlist, "Spatial Lag Residuals and Their Spatial Lag",
             sl_moran_i),
  "fig_11_sl_residual_moran_scatter.png"
)
save_plot(
  perm_plot(sl_moran, "Permutation Distribution of Moran's I: Spatial Lag Residuals",
            sl_moran_i),
  "fig_12_sl_residual_moran_permutation.png"
)
save_plot(
  map_plot(shp, "SE_RSTD", "Standardized Spatial Error Residuals",
           "Philadelphia block groups, n = 1,720"),
  "fig_13_se_residual_map.png", width = 9, height = 8
)
save_plot(
  moran_plot(residuals(se), queenlist, "Spatial Error Residuals and Their Spatial Lag",
             se_moran_i),
  "fig_14_se_residual_moran_scatter.png"
)
save_plot(
  perm_plot(se_moran, "Permutation Distribution of Moran's I: Spatial Error Residuals",
            se_moran_i),
  "fig_15_se_residual_moran_permutation.png"
)

## ---- 8. Handoff -----------------------------------------------------------

# Shapefile with the Part 2 residuals and fitted values for the GWR analyst
dir.create(file.path("data", "processed_part2"),
           showWarnings = FALSE, recursive = TRUE)
# The source shapefile carries no authoritative CRS, so it is written back
# without CRS metadata rather than with a guessed one.
st_write(
  shp,
  file.path("data", "processed_part2", "RegressionFinal_Part2.shp"),
  append = FALSE,
  quiet = TRUE
)
write.csv(
  st_drop_geometry(shp),
  file.path("data", "processed_part2", "RegressionFinal_Part2_attributes.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  File = c(
    "table_09_lm_diagnostics_ols.csv",
    "table_10_sl_coefficients.csv",
    "table_11_sl_fit_statistics.csv",
    "table_12_sl_diagnostics.csv",
    "table_13_se_coefficients.csv",
    "table_14_se_fit_statistics.csv",
    "table_15_se_diagnostics.csv",
    "table_16_model_comparison_ols_sl_se.csv",
    "sl_model_summary.txt",
    "se_model_summary.txt",
    "fig_10_sl_residual_map.png",
    "fig_11_sl_residual_moran_scatter.png",
    "fig_12_sl_residual_moran_permutation.png",
    "fig_13_se_residual_map.png",
    "fig_14_se_residual_moran_scatter.png",
    "fig_15_se_residual_moran_permutation.png"
  ),
  Content = c(
    "Lagrange multiplier diagnostics for the OLS residuals",
    "Spatial Lag coefficient table",
    "Spatial Lag fit statistics, rho and LR test against OLS",
    "Spatial Lag residual diagnostics and residual Moran's I",
    "Spatial Error coefficient table",
    "Spatial Error fit statistics, lambda and LR test against OLS",
    "Spatial Error residual diagnostics and residual Moran's I",
    "OLS, Spatial Lag and Spatial Error side-by-side",
    "Full summary() output for the Spatial Lag model",
    "Full summary() output for the Spatial Error model",
    "Map of standardized Spatial Lag residuals",
    "Moran scatter for Spatial Lag residuals",
    "Permutation distribution of Spatial Lag residual Moran's I",
    "Map of standardized Spatial Error residuals",
    "Moran scatter for Spatial Error residuals",
    "Permutation distribution of Spatial Error residual Moran's I"
  )
)
write_table(manifest, "output_manifest_part2.csv")

handoff <- c(
  "Homework 2 Part 2 handoff summary",
  "",
  "Shared model: LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV",
  "Shared weights: data/RegressionData_queen.gal, row-standardised Queen",
  "",
  sprintf("LM error statistic: %.3f (p = %.3g)",
          unname(lm_rs$RSerr$statistic), lm_rs$RSerr$p.value),
  sprintf("LM lag statistic: %.3f (p = %.3g)",
          unname(lm_rs$RSlag$statistic), lm_rs$RSlag$p.value),
  sprintf("Robust LM error statistic: %.3f (p = %.3g)",
          unname(lm_rs$adjRSerr$statistic), lm_rs$adjRSerr$p.value),
  sprintf("Robust LM lag statistic: %.3f (p = %.3g)",
          unname(lm_rs$adjRSlag$statistic), lm_rs$adjRSlag$p.value),
  "",
  sprintf("Spatial Lag rho: %.6f (z = %.3f, p = %.3g)",
          sl_rho, sl_rho_z, sl_rho_p),
  sprintf("Spatial Lag log likelihood: %.4f; AIC: %.4f", as.numeric(logLik(sl)), AIC(sl)),
  sprintf("Spatial Lag LR test vs OLS: %.3f (df = 1, p = %.3g)",
          unname(sl_lr$statistic), sl_lr$p.value),
  sprintf("Spatial Lag residual Moran's I: %.6f (999 permutations, two-sided p = %.3f, one-sided p = %.3f)",
          sl_moran_i, sl_moran$p.value, sl_moran_gt$p.value),
  "",
  sprintf("Spatial Error lambda: %.6f (z = %.3f, p = %.3g)",
          se_lambda, se_lambda_z, se_lambda_p),
  sprintf("Spatial Error log likelihood: %.4f; AIC: %.4f",
          as.numeric(logLik(se)), AIC(se)),
  sprintf("Spatial Error LR test vs OLS: %.3f (df = 1, p = %.3g)",
          unname(se_lr$statistic), se_lr$p.value),
  sprintf("Spatial Error residual Moran's I: %.6f (999 permutations, two-sided p = %.3f, one-sided p = %.3f)",
          se_moran_i, se_moran$p.value, se_moran_gt$p.value),
  "",
  sprintf("OLS log likelihood: %.4f; AIC: %.4f", ols_ll, ols_aic),
  sprintf("OLS residual Moran's I: %.6f (999 permutations, two-sided p = %.3f, one-sided p = %.3f)",
          unname(ols_moran$statistic), ols_moran$p.value, ols_moran_gt$p.value),
  "",
  "Interpretive handoff:",
  "Both spatial models dominate OLS on log likelihood and AIC.",
  "Every residual Moran's I changes from strongly positive under OLS to small",
  "and negative under both spatial models, so the positive spatial",
  "autocorrelation in the OLS residuals is absorbed by rho and lambda.",
  "The Spatial Lag model has the lower AIC and the higher log likelihood, and",
  "the robust LM diagnostics point to lag rather than error dependence, so SL",
  "is the preferred of the two Part 2 models.",
  "",
  "Reproduction note:",
  sprintf("Moran permutation tests use %d permutations with set.seed(%d).",
          permutations, perm_seed),
  "Residual Moran tests are reported two-sided because the spatial-model",
  "residuals have small negative autocorrelation; Part 1 reported a one-sided",
  "greater-than p-value for the OLS residuals. Each tail uses (count + 1)/1000.",
  "Two-sided p is twice the smaller tail probability, capped at one."
)
writeLines(handoff, file.path("output", "analysis_handoff_summary_part2.txt"))

cat("\nPart 2 outputs written.\n")
sessionInfo()
writeLines(
  capture.output(sessionInfo()),
  file.path("output", "session_info_part2.txt")
)
