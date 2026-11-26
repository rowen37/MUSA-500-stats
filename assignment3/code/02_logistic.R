# HW3 Task 2: multivariable logistic regression, classification performance,
# ROC/AUC, the reduced model, and the AIC comparison.
# Run from the assignment3 folder: Rscript code/02_logistic.R
# Or supply the CSV path: Rscript code/02_logistic.R "/path/to/Logistic Regression Data.csv"
# Required course package (install once if needed):
# install.packages("pROC")

# 1. Locate the assignment folder and import the original data.
# Task 1 does not fit a logistic model; both tasks start from this same CSV.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
  project_dir <- dirname(dirname(script_path))
} else {
  project_dir <- getwd()
}
data_path <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  file.path(project_dir, "Logistic Regression Data.csv")
setwd(project_dir)
if (!file.exists(data_path)) {
  stop("Original CSV not found. Put 'Logistic Regression Data.csv' in the assignment folder, or pass its path to Rscript.")
}
mydata <- read.csv(data_path, stringsAsFactors = FALSE)
binary_vars <- c("FATAL_OR_M", "OVERTURNED", "CELL_PHONE", "SPEEDING",
                 "AGGRESSIVE", "DRIVER1617", "DRIVER65PLUS")
continuous_vars <- c("PCTBACHMOR", "MEDHHINC")
predictors <- c(binary_vars, continuous_vars)
variable_labels <- c(
  FATAL_OR_M = "Fatality or major injury",
  OVERTURNED = "Overturned vehicle",
  CELL_PHONE = "Driver using a cell phone",
  SPEEDING = "Speeding vehicle",
  AGGRESSIVE = "Aggressive driving",
  DRIVER1617 = "At least one driver aged 16 or 17",
  DRIVER65PLUS = "At least one driver aged 65 or older",
  PCTBACHMOR = "Block-group adults aged 25+ with at least a bachelor's degree",
  MEDHHINC = "Block-group median household income"
)
required_vars <- c("DRINKING_D", predictors)
missing_vars <- setdiff(required_vars, names(mydata))
if (length(missing_vars)) stop("Missing columns: ", paste(missing_vars, collapse = ", "))
for (v in c("DRINKING_D", binary_vars)) {
  if (any(!is.na(mydata[[v]]) & !mydata[[v]] %in% c(0, 1))) {
    stop(v, " must contain only 0, 1, or missing values.")
  }
}
if (!all(c(0, 1) %in% mydata$DRINKING_D)) stop("Both drinking-status groups are needed.")

# Complete-case analysis: every model below uses the same rows.
model_data <- mydata[complete.cases(mydata[required_vars]), ]
cat("HW3 Task 2\n", R.version.string, "\n")
cat("Rows in source CSV:", nrow(mydata), "\n")
cat("Rows used in the models (complete cases on DRINKING_D and nine predictors):",
    nrow(model_data), "\n")
cat("DUI cases in the model data:", sum(model_data$DRINKING_D == 1), "of",
    nrow(model_data), sprintf("(%.2f%%)\n", 100 * mean(model_data$DRINKING_D == 1)))

dir.create("output", showWarnings = FALSE)
write_table <- function(x, name) {
  write.csv(x, file.path("output", name), row.names = FALSE, na = "")
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Install the required package first: install.packages('pROC')")
}
library(pROC)
cat("pROC version:", as.character(packageVersion("pROC")), "\n")

display_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
number <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")
p_text <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.4f", p))

# 2. Full model: DRINKING_D regressed on all nine predictors (step 3.a).
full_formula <- DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING +
  AGGRESSIVE + DRIVER1617 + DRIVER65PLUS + PCTBACHMOR + MEDHHINC
mylogit1 <- glm(full_formula, data = model_data, family = binomial(link = "logit"))
cat("\n================ Full model: summary() ================\n")
print(summary(mylogit1))

# Odds ratios and 95% confidence intervals, merged with the coefficient matrix
# (step 3.a.ii). confint() on a glm uses the profile-likelihood method shown in
# the course slides; the Wald intervals are exported as a check.
coef_mat <- summary(mylogit1)$coefficients
or <- exp(coef(mylogit1))
ci_profile <- exp(confint(mylogit1))
ci_wald <- exp(confint.default(mylogit1))
merged1 <- data.frame(
  term = rownames(coef_mat),
  beta = coef_mat[, "Estimate"],
  std_error = coef_mat[, "Std. Error"],
  z_value = coef_mat[, "z value"],
  p_value = coef_mat[, "Pr(>|z|)"],
  p_display = display_p(coef_mat[, "Pr(>|z|)"]),
  odds_ratio = as.numeric(or),
  or_ci_low = ci_profile[, 1],
  or_ci_high = ci_profile[, 2],
  or_ci_low_wald = ci_wald[, 1],
  or_ci_high_wald = ci_wald[, 2],
  significant_05 = coef_mat[, "Pr(>|z|)"] < 0.05,
  row.names = NULL
)
write_table(merged1, "table_05_logit_full_model.csv")
cat("\n================ Full model: beta, statistic, p, OR and 95% CI ================\n")
print(data.frame(
  term = merged1$term,
  beta = round(merged1$beta, 4),
  se = round(merged1$std_error, 4),
  z = round(merged1$z_value, 3),
  p = merged1$p_display,
  OR = round(merged1$odds_ratio, 4),
  OR_low = round(merged1$or_ci_low, 4),
  OR_high = round(merged1$or_ci_high, 4),
  row.names = NULL
))

# Fitted values and global fit statistics.
model_data$p_hat <- predict(mylogit1, type = "response")
loglik <- as.numeric(logLik(mylogit1))
null_dev <- mylogit1$null.deviance
resid_dev <- mylogit1$deviance
fit_summary <- data.frame(
  n = nrow(model_data),
  events = sum(model_data$DRINKING_D == 1),
  null_deviance = null_dev,
  residual_deviance = resid_dev,
  mcfadden_pseudo_r2 = 1 - loglik / as.numeric(logLik(update(mylogit1, . ~ 1))),
  aic = AIC(mylogit1),
  bic = BIC(mylogit1)
)
write_table(fit_summary, "table_06_full_model_fit.csv")
cat("\nFull model fit summary:\n")
print(fit_summary, row.names = FALSE, digits = 5)

# 3. Classification performance across the assigned cut-offs (step 3.a.iii).
cutoffs <- c(0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.50)
classify <- function(y, p, cutoff) {
  yhat <- as.integer(p >= cutoff)
  tp <- sum(y == 1 & yhat == 1)
  fn <- sum(y == 1 & yhat == 0)
  tn <- sum(y == 0 & yhat == 0)
  fp <- sum(y == 0 & yhat == 1)
  data.frame(cutoff = cutoff, true_positive = tp, false_negative = fn,
             true_negative = tn, false_positive = fp,
             sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
             misclassification_rate = (fp + fn) / length(y))
}
cutoff_table <- do.call(rbind, lapply(cutoffs, function(c0)
  classify(model_data$DRINKING_D, model_data$p_hat, c0)))
cutoff_table$lowest_misclassification <- cutoff_table$misclassification_rate ==
  min(cutoff_table$misclassification_rate)
write_table(cutoff_table, "table_07_cutoff_performance.csv")
cat("\n================ Classification by cut-off (full model) ================\n")
print(data.frame(cutoff = cutoff_table$cutoff,
                 sensitivity = round(cutoff_table$sensitivity, 4),
                 specificity = round(cutoff_table$specificity, 4),
                 misclassification = round(cutoff_table$misclassification_rate, 4),
                 lowest = cutoff_table$lowest_misclassification,
                 row.names = NULL))
best_cutoff_min_error <- cutoff_table$cutoff[which.min(cutoff_table$misclassification_rate)]
cat("Cut-off with the lowest misclassification rate:",
    best_cutoff_min_error,
    sprintf("(%.2f%% misclassified)\n", 100 * min(cutoff_table$misclassification_rate)))

# 4. ROC curve, optimal cut-off and AUC (steps 3.a.iv and 3.a.v).
# The optimal cut-off is the one that minimises the distance to the upper-left
# corner of the ROC curve, as required by the report outline.
roc_full <- roc(model_data$DRINKING_D, model_data$p_hat, levels = c(0, 1),
                direction = "<", quiet = TRUE)
roc_points <- coords(roc_full, x = "all", ret = c("threshold", "sensitivity", "specificity"),
                     transpose = FALSE)
distance_ul <- sqrt((1 - roc_points$sensitivity)^2 + (1 - roc_points$specificity)^2)
optimal_row <- roc_points[which.min(distance_ul), ]
auc_full <- as.numeric(auc(roc_full))
auc_ci <- as.numeric(ci.auc(roc_full))
cat("\n================ ROC and AUC (full model) ================\n")
cat("AUC:", round(auc_full, 4), " 95% CI (DeLong):", round(auc_ci[1], 4), "-",
    round(auc_ci[3], 4), "\n")
cat("Optimal cut-off (minimum distance to the upper-left corner):",
    round(optimal_row$threshold, 4), "\n")
cat("  sensitivity =", round(optimal_row$sensitivity, 4),
    " specificity =", round(optimal_row$specificity, 4), "\n")
optimal_performance <- classify(model_data$DRINKING_D, model_data$p_hat, optimal_row$threshold)
optimal_table <- data.frame(
  method = "Minimum distance to the upper-left corner of the ROC curve",
  cutoff = optimal_row$threshold,
  sensitivity = optimal_row$sensitivity,
  specificity = optimal_row$specificity,
  misclassification_rate = optimal_performance$misclassification_rate,
  auc = auc_full, auc_ci_low = auc_ci[1], auc_ci_high = auc_ci[3]
)
# Youden's J is exported only as a comparison, not as the reported choice.
youden_row <- coords(roc_full, x = "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
optimal_table <- rbind(
  optimal_table,
  data.frame(method = "Youden index (comparison only)", cutoff = youden_row$threshold,
             sensitivity = youden_row$sensitivity, specificity = youden_row$specificity,
             misclassification_rate = classify(model_data$DRINKING_D, model_data$p_hat,
                                               youden_row$threshold)$misclassification_rate,
             auc = auc_full, auc_ci_low = auc_ci[1], auc_ci_high = auc_ci[3])
)
write_table(optimal_table, "table_08_roc_optimal_cutoff.csv")

# Draw the curve directly from the coordinates so that the horizontal axis is
# 1 - specificity measured left to right from 0 to 1 (pROC's own axes differ).
png(file.path("output", "fig_01_roc_curve.png"), width = 1500, height = 1250, res = 200)
par(mar = c(4.6, 4.6, 3.2, 1.8))
fpr <- 1 - roc_points$specificity
tpr <- roc_points$sensitivity
ord <- order(fpr, tpr)
plot(fpr[ord], tpr[ord], type = "l", col = "#1f4e79", lwd = 2,
     xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", las = 1,
     main = sprintf("ROC curve: full logistic model (AUC = %.3f)", auc_full),
     xlab = "1 - Specificity (false positive rate)", ylab = "Sensitivity (true positive rate)")
abline(0, 1, col = "grey55", lty = 2)
points(1 - optimal_row$specificity, optimal_row$sensitivity, pch = 19, col = "#c00000", cex = 1.3)
text(1 - optimal_row$specificity, optimal_row$sensitivity,
     labels = sprintf(" optimal cut-off = %.3f", optimal_row$threshold),
     pos = 4, col = "#c00000", cex = 0.85, offset = 0.6)
box()
dev.off()
cat("ROC image written to output/fig_01_roc_curve.png\n")

# 5. Reduced model without the two block-group predictors (step 3.b).
reduced_formula <- DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING +
  AGGRESSIVE + DRIVER1617 + DRIVER65PLUS
mylogit2 <- glm(reduced_formula, data = model_data, family = binomial(link = "logit"))
cat("\n================ Reduced model (no PCTBACHMOR, no MEDHHINC): summary() ================\n")
print(summary(mylogit2))

coef_mat2 <- summary(mylogit2)$coefficients
merged2 <- data.frame(
  term = rownames(coef_mat2),
  beta = coef_mat2[, "Estimate"],
  std_error = coef_mat2[, "Std. Error"],
  z_value = coef_mat2[, "z value"],
  p_value = coef_mat2[, "Pr(>|z|)"],
  p_display = display_p(coef_mat2[, "Pr(>|z|)"]),
  odds_ratio = exp(coef(mylogit2)),
  or_ci_low = exp(confint(mylogit2))[, 1],
  or_ci_high = exp(confint(mylogit2))[, 2],
  or_ci_low_wald = exp(confint.default(mylogit2))[, 1],
  or_ci_high_wald = exp(confint.default(mylogit2))[, 2],
  significant_05 = coef_mat2[, "Pr(>|z|)"] < 0.05,
  row.names = NULL
)
write_table(merged2, "table_09_logit_reduced_model.csv")
cat("\n================ Reduced model: beta, statistic, p, OR and 95% CI ================\n")
print(data.frame(
  term = merged2$term,
  beta = round(merged2$beta, 4),
  se = round(merged2$std_error, 4),
  z = round(merged2$z_value, 3),
  p = merged2$p_display,
  OR = round(merged2$odds_ratio, 4),
  OR_low = round(merged2$or_ci_low, 4),
  OR_high = round(merged2$or_ci_high, 4),
  row.names = NULL
))

# Same classification diagnostics for the reduced model, for a fair comparison.
model_data$p_hat2 <- predict(mylogit2, type = "response")
cutoff_table2 <- do.call(rbind, lapply(cutoffs, function(c0)
  classify(model_data$DRINKING_D, model_data$p_hat2, c0)))
cutoff_table2$lowest_misclassification <- cutoff_table2$misclassification_rate ==
  min(cutoff_table2$misclassification_rate)
write_table(cutoff_table2, "table_10_cutoff_performance_reduced.csv")
roc_reduced <- roc(model_data$DRINKING_D, model_data$p_hat2, levels = c(0, 1),
                   direction = "<", quiet = TRUE)
cat("\nReduced model AUC:", round(as.numeric(auc(roc_reduced)), 4),
    " lowest-misclassification cut-off:", cutoff_table2$cutoff[which.min(cutoff_table2$misclassification_rate)], "\n")

# 6. AIC comparison (step 3.c). Lower AIC = better model.
aic_table <- data.frame(
  model = c("Full model (nine predictors)", "Reduced model (seven binary predictors only)"),
  formula = c("DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING + AGGRESSIVE + DRIVER1617 + DRIVER65PLUS + PCTBACHMOR + MEDHHINC",
              "DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING + AGGRESSIVE + DRIVER1617 + DRIVER65PLUS"),
  n = nrow(model_data),
  parameters = c(length(coef(mylogit1)), length(coef(mylogit2))),
  log_likelihood = c(loglik, as.numeric(logLik(mylogit2))),
  aic = c(AIC(mylogit1), AIC(mylogit2)),
  bic = c(BIC(mylogit1), BIC(mylogit2)),
  delta_aic = c(0, AIC(mylogit2) - AIC(mylogit1)),
  stringsAsFactors = FALSE
)
aic_table$better_model <- aic_table$aic == min(aic_table$aic)
write_table(aic_table, "table_11_aic_comparison.csv")
cat("\n================ AIC comparison ================\n")
cat("AIC(full) =", round(AIC(mylogit1), 4), "  AIC(reduced) =", round(AIC(mylogit2), 4),
    "  difference =", round(AIC(mylogit2) - AIC(mylogit1), 4), "\n")
cat("Lower AIC indicates the better model:", aic_table$model[which.min(aic_table$aic)], "\n")

# 7. Which coefficients changed sign or significance between the two models?
comparison <- merge(
  merged1[merged1$term != "(Intercept)", c("term", "beta", "p_value", "odds_ratio")],
  merged2[merged2$term != "(Intercept)", c("term", "beta", "p_value", "odds_ratio")],
  by = "term", suffixes = c("_full", "_reduced")
)
comparison$beta_change <- comparison$beta_reduced - comparison$beta_full
sig_full <- comparison$p_value_full < 0.05
sig_reduced <- comparison$p_value_reduced < 0.05
comparison$significance_change <- ifelse(
  sig_full,
  ifelse(sig_reduced, "significant in both", "significant only in full model"),
  ifelse(sig_reduced, "significant only in reduced model", "not significant in either")
)
write_table(comparison, "table_12_model_comparison.csv")
cat("\n================ Coefficient comparison across models ================\n")
print(data.frame(term = comparison$term,
                 beta_full = round(comparison$beta_full, 4),
                 beta_reduced = round(comparison$beta_reduced, 4),
                 p_full = display_p(comparison$p_value_full),
                 p_reduced = display_p(comparison$p_value_reduced),
                 status = comparison$significance_change,
                 row.names = NULL))

# 8. Handoff to Task 3 / the report writer.
p_label <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.4f", p))
notes <- c(
  "HW3 Task 2: logistic regression handoff",
  paste("Source CSV:", basename(data_path)),
  paste("Rows in source data:", nrow(mydata), "; rows used:", nrow(model_data)),
  "Both models start from the same original CSV used by Task 1 and use the same complete-case rows.",
  "Full model: DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING + AGGRESSIVE + DRIVER1617 + DRIVER65PLUS + PCTBACHMOR + MEDHHINC.",
  "Reduced model: the same minus PCTBACHMOR and MEDHHINC.",
  "Link: logit (binomial). Odds ratios = exp(beta). CIs: profile-likelihood confint(); Wald CIs also exported for a check.",
  "Model fit: McFadden pseudo-R2, null and residual deviance, AIC, BIC exported in table_06.",
  sprintf("Cut-off with the lowest misclassification rate: %.2f (%.2f%% misclassified).",
          best_cutoff_min_error, 100 * min(cutoff_table$misclassification_rate)),
  sprintf("ROC optimal cut-off (minimum distance to upper-left corner): %.4f (sensitivity %.4f, specificity %.4f).",
          optimal_row$threshold, optimal_row$sensitivity, optimal_row$specificity),
  sprintf("AUC = %.4f (95%% CI %.4f-%.4f, DeLong).", auc_full, auc_ci[1], auc_ci[3]),
  sprintf("AIC: full = %.2f, reduced = %.2f; the %s has the lower AIC.",
          AIC(mylogit1), AIC(mylogit2),
          ifelse(AIC(mylogit2) < AIC(mylogit1), "reduced model", "full model")),
  "",
  "Significant predictors in the full model (alpha = 0.05):",
  paste(merged1$term[merged1$term != "(Intercept)" & merged1$significant_05],
        p_label(merged1$p_value[merged1$term != "(Intercept)" & merged1$significant_05]),
        sep = " | "),
  "Not significant in the full model:",
  paste(merged1$term[merged1$term != "(Intercept)" & !merged1$significant_05],
        p_label(merged1$p_value[merged1$term != "(Intercept)" & !merged1$significant_05]),
        sep = " | "),
  "",
  sprintf("Reduced model AUC = %.4f, so dropping the two block-group predictors also costs a little discrimination.",
          as.numeric(auc(roc_reduced))),
  "",
  "Note: these are conditional associations in one model, not causal effects.",
  "The block-group variables describe the crash location, not the driver; their coefficients are not driver-level effects.",
  "DRINKING_D is rare (about 5.7%), so classification at any cut-off above the observed prevalence is weak; discuss this in the report.",
  paste("pROC", as.character(packageVersion("pROC")), "| R", R.version.string)
)
writeLines(notes, file.path("output", "task2_handoff.txt"))

# 9. Readable tables and a draft Results text for the report.
markdown_table <- function(headers, rows) {
  c(paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |"),
    apply(as.matrix(rows), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
model_rows <- data.frame(
  merged1$term,
  number(merged1$beta, 4),
  number(merged1$std_error, 4),
  number(merged1$z_value, 3),
  merged1$p_display,
  number(merged1$odds_ratio, 3),
  sprintf("%s - %s", number(merged1$or_ci_low, 3), number(merged1$or_ci_high, 3))
)
cutoff_rows <- data.frame(
  sprintf("%.2f", cutoff_table$cutoff),
  number(cutoff_table$sensitivity, 4),
  number(cutoff_table$specificity, 4),
  sprintf("%s%s", number(100 * cutoff_table$misclassification_rate, 2),
          ifelse(cutoff_table$lowest_misclassification, " % (lowest)", " %"))
)
reduced_rows <- data.frame(
  merged2$term,
  number(merged2$beta, 4),
  number(merged2$std_error, 4),
  number(merged2$z_value, 3),
  merged2$p_display,
  number(merged2$odds_ratio, 3),
  sprintf("%s - %s", number(merged2$or_ci_low, 3), number(merged2$or_ci_high, 3))
)
auc_grade <- function(a) {
  if (a >= 0.9) "excellent" else if (a >= 0.8) "good" else
    if (a >= 0.7) "fair" else if (a >= 0.6) "poor" else "failing"
}
aic_lower <- if (AIC(mylogit2) < AIC(mylogit1)) "reduced" else "full"
aic_sentence <- if (aic_lower == "reduced") {
  sprintf("The reduced model has the lower AIC (%.2f vs %.2f; Delta AIC = %.2f), so dropping PCTBACHMOR and MEDHHINC is supported by this criterion.", AIC(mylogit2), AIC(mylogit1), AIC(mylogit1) - AIC(mylogit2))
} else {
  sprintf("The full model has the lower AIC (%.2f vs %.2f; Delta AIC = %.2f), only %.2f units better, so AIC does not justify dropping PCTBACHMOR and MEDHHINC. BIC penalises the two extra parameters more heavily and is lower for the reduced model (%.2f vs %.2f), so the two criteria disagree and the AIC advantage of the full model is too small to call a meaningful improvement.",
          AIC(mylogit1), AIC(mylogit2), AIC(mylogit2) - AIC(mylogit1),
          AIC(mylogit2) - AIC(mylogit1), BIC(mylogit2), BIC(mylogit1))
}
or_note <- paste0(
  "Note on the two continuous predictors: MEDHHINC is in dollars, so its odds ratio is per $1 and prints as 1.000. ",
  sprintf("Per $1,000 of median household income the odds ratio is %.3f (coefficient re-scaled to thousands = %.4f, %s). ",
          exp(1000 * merged1$beta[merged1$term == "MEDHHINC"]), 1000 * merged1$beta[merged1$term == "MEDHHINC"],
          p_text(merged1$p_value[merged1$term == "MEDHHINC"])),
  sprintf("For PCTBACHMOR, a 10-percentage-point increase multiplies the odds by exp(10 * %.6f) = %.4f.", merged1$beta[merged1$term == "PCTBACHMOR"], exp(10 * merged1$beta[merged1$term == "PCTBACHMOR"]))
)
p_text <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.4f", p))
significant_full <- merged1$term[merged1$term != "(Intercept)" & merged1$significant_05]
nonsignificant_full <- merged1$term[merged1$term != "(Intercept)" & !merged1$significant_05]
results <- c(
  "# HW3 Task 2: Logistic Regression", "",
  "## Data and methods", "",
  sprintf("The models use the same %s Philadelphia crash records as Task 1; no rows are added or dropped beyond complete cases on DRINKING_D and the nine predictors.", number(nrow(model_data), 0)),
  "Fitted by maximum likelihood with `glm(..., family = binomial)`. Odds ratios are exp(beta); 95% intervals use the profile-likelihood method, and Wald intervals are exported as a check. Classification thresholds are applied to fitted probabilities p-hat.", "",
  "## Table 5. Full model: coefficients, odds ratios and 95% CI", "",
  markdown_table(c("Term", "beta", "SE", "z", "p", "OR", "95% CI for OR"), model_rows), "",
  or_note, "",
  sprintf("McFadden pseudo-R2 = %.4f; AIC = %.2f; null deviance = %.1f; residual deviance = %.1f.",
          fit_summary$mcfadden_pseudo_r2, fit_summary$aic, fit_summary$null_deviance,
          fit_summary$residual_deviance), "",
  sprintf("At alpha = 0.05 the full model flags %s. %s are not significant.",
          ifelse(length(significant_full), paste(significant_full, collapse = ", "), "no predictor"),
          ifelse(length(nonsignificant_full), paste(nonsignificant_full, collapse = ", "), "No predictor")), "",
  "## Table 6. Sensitivity, specificity and misclassification rate by cut-off", "",
  markdown_table(c("Cut-off", "Sensitivity", "Specificity", "Misclassification rate"), cutoff_rows), "",
  sprintf("The lowest misclassification rate occurs at cut-off %.2f (%.2f%% misclassified). Because only %.2f%% of the crashes involve a drinking driver, a cut-off that high assigns almost every crash to the non-DUI class: at 0.50 only %s of %s crashes are classified as DUI (the highest fitted probability in the sample is %.3f). The minimum error rate is therefore close to the error rate of simply predicting the majority class (%.2f%%), and it does not indicate good discrimination. The error rate also jumps sharply between cut-offs 0.03 and 0.07, where the predicted probabilities are concentrated.",
          best_cutoff_min_error, 100 * min(cutoff_table$misclassification_rate),
          100 * mean(model_data$DRINKING_D),
          number(sum(model_data$p_hat >= best_cutoff_min_error), 0), number(nrow(model_data), 0),
          max(model_data$p_hat), 100 * mean(model_data$DRINKING_D == 0)), "",
  "## ROC curve, optimal cut-off and AUC", "",
  "![ROC curve for the full model](fig_01_roc_curve.png)", "",
  sprintf("AUC = %.4f (95%% CI %.4f-%.4f, DeLong), which is %s under the usual 0.5-0.6 failing / 0.6-0.7 poor / 0.7-0.8 fair / 0.8-0.9 good / above 0.9 excellent scale. The optimal cut-off that minimises the distance to the upper-left corner is %.4f, giving sensitivity %.4f and specificity %.4f; the corresponding misclassification rate is %.2f%%. Youden's J would choose %.4f instead (comparison only). Note that the cut-off that minimises the distance to the upper-left corner and the cut-off that minimises misclassification are not the same: the former balances sensitivity and specificity, the latter favours the majority class.",
          auc_full, auc_ci[1], auc_ci[3], auc_grade(auc_full), optimal_row$threshold,
          optimal_row$sensitivity, optimal_row$specificity,
          100 * optimal_performance$misclassification_rate, youden_row$threshold), "",
  "## Table 7. Reduced model (no PCTBACHMOR, no MEDHHINC)", "",
  markdown_table(c("Term", "beta", "SE", "z", "p", "OR", "95% CI for OR"), reduced_rows), "",
  "## Table 8. AIC comparison", "",
  markdown_table(c("Model", "Parameters", "Log-likelihood", "AIC", "BIC", "Delta AIC"),
                 data.frame(aic_table$model, sprintf("%d", aic_table$parameters),
                            number(aic_table$log_likelihood, 2), number(aic_table$aic, 2),
                            number(aic_table$bic, 2), number(aic_table$delta_aic, 2))), "",
  aic_sentence, "",
  "## Predictors whose significance changed between models", "",
  paste0("- **", comparison$term, ":** beta ", number(comparison$beta_full, 4), " -> ",
         number(comparison$beta_reduced, 4), "; ", comparison$significance_change, "."), "",
  "## Draft text for the Results section (Section 3c)", "",
  sprintf("The full logistic model estimated on %s crashes is significant overall relative to the intercept-only model (residual deviance %.1f vs null deviance %.1f on %d degrees of freedom). McFadden's pseudo-R2 is %.3f, so the nine predictors explain only a modest share of the variation in the drinking-driver indicator.",
          number(nrow(model_data), 0), fit_summary$residual_deviance, fit_summary$null_deviance,
          mylogit1$df.residual, fit_summary$mcfadden_pseudo_r2), "",
  sprintf("In the full model, %s. The remaining coefficients are not statistically distinguishable from zero at the 5%% level, so their odds ratios should be read as imprecise rather than as evidence of no relationship.",
          paste0(sprintf("%s is significant (OR = %.2f, 95%% CI %.2f-%.2f; %s)",
                         significant_full,
                         merged1$odds_ratio[match(significant_full, merged1$term)],
                         merged1$or_ci_low[match(significant_full, merged1$term)],
                         merged1$or_ci_high[match(significant_full, merged1$term)],
                         p_text(merged1$p_value[match(significant_full, merged1$term)])),
                 collapse = "; ")), "",
  sprintf("Because alcohol-involved crashes are only %.2f%% of the records, a high predicted probability is rare. Across the assigned cut-offs the misclassification rate is minimised at %.2f (%.2f%% misclassified), while sensitivity and specificity trade off against each other. The ROC curve reaches an AUC of %.3f, which is %s on the standard scale, and the cut-off closest to the upper-left corner is %.4f.",
          100 * mean(model_data$DRINKING_D), best_cutoff_min_error,
          100 * min(cutoff_table$misclassification_rate), auc_full, auc_grade(auc_full), optimal_row$threshold), "",
  sprintf("Dropping PCTBACHMOR and MEDHHINC leaves the seven binary coefficients essentially unchanged (largest change in beta = %.4f) and no predictor switches significance. %s",
          max(abs(comparison$beta_change)), aic_sentence), "",
  "## Sources", "",
  "- Logistic Regression Data.csv, provided with the assignment; one row per crash.",
  "- HW 3 - Instructions.pdf, pp. 7-8 (logistic regression, cut-offs, ROC, AUC, reduced model, AIC) and pp. 11-12 (Results requirements).",
  "- output/task1_handoff.txt (same CSV, same nine predictors, no bivariate screening).",
  "- pROC package documentation for `roc`, `coords`, `auc` and `ci.auc`."
)
writeLines(results, file.path("output", "task2_results.md"), useBytes = TRUE)
cat("\nTables 5-12, fig_01_roc_curve.png, task2_results.md and task2_handoff.txt written to:",
    file.path(getwd(), "output"), "\n")
