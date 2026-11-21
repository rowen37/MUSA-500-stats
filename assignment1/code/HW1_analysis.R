# Reproducible analysis for CPLN 671 / MUSA 500 Homework 1
# Run from the code/ directory. Outputs are written to ../output/.

suppressPackageStartupMessages({
  library(sf)
  library(classInt)
})

data_path <- "../data/RegressionData.csv"
shape_path <- "../data/map_RegressionData/RegressionData.shp"
output_dir <- "../output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

navy <- "#173B57"
blue <- "#2C7FB8"
light_blue <- "#A6CEE3"
accent <- "#C94C24"
grid_gray <- "#D9E1E8"
map_colors <- colorRampPalette(
  c("#08306B", "#6BAED6", "#F7FBFF", "#FDBE85", "#B30000")
)(5)

open_png <- function(filename, width, height, res = 300) {
  png(file.path(output_dir, filename), width = width, height = height,
      units = "in", res = res, bg = "white")
}

set_plot_style <- function(mar = c(4.2, 4.5, 2.4, 1.2)) {
  par(
    family = "serif", fg = navy, col.axis = navy, col.lab = navy,
    col.main = navy, las = 1, bty = "l", mar = mar,
    mgp = c(2.45, 0.7, 0), tcl = -0.25
  )
}

df <- read.csv(data_path)
stopifnot(nrow(df) == 1720, !anyNA(df))

df$LNMEDHVAL <- log(df$MEDHVAL)
df$LNPCTBACHMOR <- log1p(df$PCTBACHMOR)
df$LNNBELPOV100 <- log1p(df$NBELPOV100)
df$LNPCTVACANT <- log1p(df$PCTVACANT)
df$LNPCTSINGLES <- log1p(df$PCTSINGLES)

summary_vars <- c(
  "MEDHVAL", "PCTBACHMOR", "NBELPOV100", "PCTVACANT", "PCTSINGLES"
)
summary_labels <- c(
  "Median House Value ($)", "% Bachelor's Degree or Higher",
  "Households Living in Poverty", "% Vacant Housing Units",
  "% Detached Single-Family Units"
)
summary_table <- data.frame(
  Variable = summary_labels,
  Mean = sapply(df[summary_vars], mean),
  SD = sapply(df[summary_vars], sd),
  row.names = NULL
)
write.csv(summary_table, file.path(output_dir, "table_01_summary_statistics.csv"),
          row.names = FALSE)

predictors <- c("PCTVACANT", "PCTSINGLES", "PCTBACHMOR", "LNNBELPOV100")
cor_matrix <- cor(df[predictors])
write.csv(cor_matrix, file.path(output_dir, "table_02_correlations.csv"))

fit <- lm(
  LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV100,
  data = df
)
fit_summary <- summary(fit)
coef_table <- as.data.frame(coef(fit_summary))
coef_table$Term <- rownames(coef_table)
rownames(coef_table) <- NULL
coef_table <- coef_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
write.csv(coef_table, file.path(output_dir, "table_03_ols_coefficients.csv"),
          row.names = FALSE)

anova_raw <- anova(fit)
reg_ss <- sum(anova_raw$`Sum Sq`[seq_len(length(predictors))])
err_ss <- deviance(fit)
total_ss <- reg_ss + err_ss
anova_table <- data.frame(
  Source = c("Regression", "Error", "Total"),
  Df = c(length(predictors), df.residual(fit), nrow(df) - 1),
  Sum_Sq = c(reg_ss, err_ss, total_ss),
  Mean_Sq = c(reg_ss / length(predictors), err_ss / df.residual(fit), NA),
  F = c(unname(fit_summary$fstatistic[1]), NA, NA),
  p_value = c(pf(fit_summary$fstatistic[1], length(predictors),
                 df.residual(fit), lower.tail = FALSE), NA, NA)
)
write.csv(anova_table, file.path(output_dir, "table_04_anova.csv"),
          row.names = FALSE)

drop_table <- drop1(fit, test = "none")
drop_table$Model <- rownames(drop_table)
rownames(drop_table) <- NULL
write.csv(drop_table, file.path(output_dir, "table_05_stepwise_aic.csv"),
          row.names = FALSE)

# DAAG::CVlm assigns folds as sample(n) %% m + 1. This local implementation
# gives the same predictions without installing packages during analysis.
cv_predictions <- function(formula, data, m = 5, seed = 671) {
  set.seed(seed)
  n <- nrow(data)
  fold <- sample(n) %% m + 1
  predictions <- numeric(n)
  for (i in sort(unique(fold))) {
    train <- fold != i
    model_i <- lm(formula, data = data[train, , drop = FALSE])
    predictions[!train] <- predict(model_i, newdata = data[!train, , drop = FALSE])
  }
  list(predictions = predictions, fold = fold)
}

cv_full <- cv_predictions(
  LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV100,
  df
)
cv_simple <- cv_predictions(LNMEDHVAL ~ PCTVACANT + MEDHHINC, df)
rmse_full <- sqrt(mean((df$LNMEDHVAL - cv_full$predictions)^2))
rmse_simple <- sqrt(mean((df$LNMEDHVAL - cv_simple$predictions)^2))
cv_table <- data.frame(
  Model = c("Original four-predictor model", "PCTVACANT + MEDHHINC model"),
  RMSE = c(rmse_full, rmse_simple)
)
write.csv(cv_table, file.path(output_dir, "table_06_cross_validation.csv"),
          row.names = FALSE)

df$fitted <- fitted(fit)
df$residual <- residuals(fit)
df$rstandard <- rstandard(fit)
write.csv(
  df[, c("POLY_ID", "LNMEDHVAL", "fitted", "residual", "rstandard")],
  file.path(output_dir, "regression_diagnostics.csv"), row.names = FALSE
)

# Figure 1: paired original and transformed distributions
original_vars <- c("MEDHVAL", "PCTBACHMOR", "NBELPOV100", "PCTVACANT", "PCTSINGLES")
log_vars <- c("LNMEDHVAL", "LNPCTBACHMOR", "LNNBELPOV100", "LNPCTVACANT", "LNPCTSINGLES")
pretty_names <- c("Median House Value", "Bachelor's Degree or Higher",
                  "Households in Poverty", "Vacant Housing Units",
                  "Detached Single-Family Units")
open_png("fig_01_distributions.png", 8.0, 11.0)
par(mfrow = c(5, 2), oma = c(0.2, 0.2, 1.0, 0.2))
for (i in seq_along(original_vars)) {
  set_plot_style(c(3.3, 3.8, 2.0, 0.8))
  hist(df[[original_vars[i]]], breaks = 28, col = light_blue, border = "white",
       main = paste(pretty_names[i], "- original"), xlab = original_vars[i])
  set_plot_style(c(3.3, 3.8, 2.0, 0.8))
  hist(df[[log_vars[i]]], breaks = 28, col = blue, border = "white",
       main = paste(pretty_names[i], "- transformed"), xlab = log_vars[i])
}
mtext("Original and Log-Transformed Distributions", outer = TRUE,
      side = 3, line = -0.1, font = 2, col = navy, cex = 1.15)
dev.off()

# Prepare spatial data without duplicated derived columns.
phila <- st_read(shape_path, quiet = TRUE)
idx <- match(phila$POLY_ID, df$POLY_ID)
stopifnot(!anyNA(idx))
phila$LNMEDHVAL <- df$LNMEDHVAL[idx]
phila$LNNBELPOV100 <- df$LNNBELPOV100[idx]
phila$rstandard <- df$rstandard[idx]

map_var <- function(sf_obj, var, title, legend_title = var) {
  intervals <- classIntervals(sf_obj[[var]], n = 5, style = "quantile")
  breaks <- intervals$brks
  classes <- cut(sf_obj[[var]], breaks = breaks, include.lowest = TRUE,
                 dig.lab = 4)
  plot(st_geometry(sf_obj), col = map_colors[classes], border = NA,
       main = title, col.main = navy)
  legend(
    "bottomright", legend = levels(classes), fill = map_colors,
    title = legend_title, cex = 0.62, bty = "o", bg = "white",
    box.col = grid_gray, x.intersp = 0.6, y.intersp = 0.85
  )
}

# Figure 2: response map
open_png("fig_02_map_response.png", 7.4, 8.8)
set_plot_style(c(0.8, 0.8, 2.0, 0.8))
map_var(phila, "LNMEDHVAL", "Log Median House Value", "LNMEDHVAL")
dev.off()

# Figure 3: predictor maps
open_png("fig_03_map_predictors.png", 10.0, 8.5)
par(mfrow = c(2, 2), oma = c(0.2, 0.2, 0.8, 0.2))
set_plot_style(c(0.5, 0.5, 1.8, 0.5))
map_var(phila, "PCTVACANT", "Vacant Housing Units", "PCTVACANT")
set_plot_style(c(0.5, 0.5, 1.8, 0.5))
map_var(phila, "PCTSINGLES", "Detached Single-Family Units", "PCTSINGLES")
set_plot_style(c(0.5, 0.5, 1.8, 0.5))
map_var(phila, "PCTBACHMOR", "Bachelor's Degree or Higher", "PCTBACHMOR")
set_plot_style(c(0.5, 0.5, 1.8, 0.5))
map_var(phila, "LNNBELPOV100", "Log Households in Poverty", "LNNBELPOV100")
mtext("Neighborhood Predictors", outer = TRUE, side = 3, line = -0.2,
      font = 2, col = navy, cex = 1.15)
dev.off()

# Figure 4: response-predictor scatterplots with a common style.
scatter_x <- c("PCTVACANT", "PCTSINGLES", "PCTBACHMOR", "LNNBELPOV100")
scatter_labels <- c("Vacant Housing Units (%)", "Detached Single-Family Units (%)",
                    "Bachelor's Degree or Higher (%)", "Log Households in Poverty")
open_png("fig_04_scatterplots.png", 9.2, 7.6)
par(mfrow = c(2, 2), oma = c(0.2, 0.2, 0.8, 0.2))
for (i in seq_along(scatter_x)) {
  set_plot_style(c(4.0, 4.4, 1.6, 0.8))
  x <- df[[scatter_x[i]]]
  plot(x, df$LNMEDHVAL, pch = 19, cex = 0.42,
       col = adjustcolor(blue, alpha.f = 0.30),
       xlab = scatter_labels[i], ylab = "LNMEDHVAL", main = "")
  grid(col = adjustcolor(grid_gray, alpha.f = 0.75), lty = 1)
  points(x, df$LNMEDHVAL, pch = 19, cex = 0.38,
         col = adjustcolor(blue, alpha.f = 0.28))
  lines(lowess(x, df$LNMEDHVAL, f = 2 / 3), col = accent, lwd = 2.2)
  legend("bottomright", legend = "LOWESS", col = accent, lwd = 2.2,
         bty = "n", cex = 0.75)
}
mtext("Log Median House Value and Model Predictors", outer = TRUE,
      side = 3, line = -0.2, font = 2, col = navy, cex = 1.15)
dev.off()

# Figure 5: standardized residual distribution with normal reference.
open_png("fig_05_residual_histogram.png", 7.6, 5.2)
set_plot_style()
h <- hist(df$rstandard, breaks = 45, probability = TRUE,
          col = light_blue, border = "white",
          main = "Distribution of Standardized Residuals",
          xlab = "Standardized residual")
curve(dnorm(x), add = TRUE, col = accent, lwd = 2.4)
abline(v = 0, col = navy, lwd = 1)
legend("bottomright", legend = c("Normal reference", "Zero"),
       col = c(accent, navy), lwd = c(2.4, 1), bty = "n", cex = 0.85)
dev.off()

# Figure 6: standardized residuals versus fitted values.
open_png("fig_06_residuals_vs_fitted.png", 7.6, 5.2)
set_plot_style()
plot(df$fitted, df$rstandard, pch = 19, cex = 0.48,
     col = adjustcolor(blue, alpha.f = 0.38),
     main = "Standardized Residuals versus Fitted Values",
     xlab = "Fitted LNMEDHVAL", ylab = "Standardized residual")
grid(col = adjustcolor(grid_gray, alpha.f = 0.75), lty = 1)
points(df$fitted, df$rstandard, pch = 19, cex = 0.43,
       col = adjustcolor(blue, alpha.f = 0.34))
abline(h = 0, col = navy, lwd = 1.2)
abline(h = c(-2, 2), col = accent, lty = 2, lwd = 1.5)
legend("bottomright", legend = c("Zero", "+/- 2 reference"),
       col = c(navy, accent), lty = c(1, 2), lwd = c(1.2, 1.5),
       bty = "n", cex = 0.82)
dev.off()

# Figure 7: standardized residual map.
open_png("fig_07_residual_map.png", 7.4, 8.8)
set_plot_style(c(0.8, 0.8, 2.0, 0.8))
map_var(phila, "rstandard", "Standardized Regression Residuals", "rstandard")
dev.off()

cat("Analysis complete.\n")
cat(sprintf("R2 = %.6f; adjusted R2 = %.6f; F = %.4f\n",
            fit_summary$r.squared, fit_summary$adj.r.squared,
            unname(fit_summary$fstatistic[1])))
cat(sprintf("RMSE: full = %.6f; simple = %.6f\n", rmse_full, rmse_simple))
cat(sprintf("|standardized residual| > 2: %d; > 3: %d\n",
            sum(abs(df$rstandard) > 2), sum(abs(df$rstandard) > 3)))
