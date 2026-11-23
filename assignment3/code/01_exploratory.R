# HW3 Task 1: data preparation and exploratory analysis
# Run from the assignment3 folder: Rscript code/01_exploratory.R
# Or supply the CSV path: Rscript code/01_exploratory.R "/path/to/Logistic Regression Data.csv"
# Required course package (install once if needed):
# install.packages("gmodels")

# 1. Locate the assignment folder and import the original data.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
  project_dir <- dirname(dirname(script_path))
} else {
  # In RStudio, set the working directory to the assignment3 folder first.
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
if (!all(vapply(mydata[required_vars], is.numeric, logical(1)))) {
  stop("Required columns must be numeric. Check the original CSV before analyzing it.")
}
for (v in c("DRINKING_D", binary_vars)) {
  if (any(!is.na(mydata[[v]]) & !mydata[[v]] %in% c(0, 1))) {
    stop(v, " must contain only 0, 1, or missing values.")
  }
}
if (!all(c(0, 1) %in% mydata$DRINKING_D)) stop("Both drinking-status groups are needed.")
dir.create("output", showWarnings = FALSE)
write_table <- function(x, name) {
  write.csv(x, file.path("output", name), row.names = FALSE, na = "")
}
if (!requireNamespace("gmodels", quietly = TRUE)) {
  stop("Install the required package first: install.packages('gmodels')")
}
library(gmodels)
cat("HW3 Task 1\n", R.version.string, "\ngmodels version:",
    as.character(packageVersion("gmodels")), "\n")

# 2. Outcome distribution. Percentages exclude missing DRINKING_D.
dv_counts <- table(factor(mydata$DRINKING_D, levels = c(0, 1)))
dv_table <- data.frame(
  DRINKING_D = c(0, 1),
  group = c("Non-DUI", "DUI"),
  n = as.integer(dv_counts),
  proportion = as.numeric(prop.table(dv_counts)),
  percent = 100 * as.numeric(prop.table(dv_counts))
)
write_table(dv_table, "table_01_dv_distribution.csv")
cat("\nDependent-variable distribution:\n")
print(dv_table, row.names = FALSE)

# 3. Binary predictors: count predictor=1 in each outcome column.
# All percentages below use the outcome-column total as the denominator.
# The Pearson chi-square test uses the ENTIRE 2 x 2 table, with no Yates correction.
binary_results <- lapply(binary_vars, function(v) {
  d <- mydata[complete.cases(mydata[c(v, "DRINKING_D")]), c(v, "DRINKING_D")]
  cat("\nPredictor:", v, "\n")
  cross_table <- gmodels::CrossTable(d[[v]], d$DRINKING_D,
                                    prop.r = FALSE, prop.t = FALSE,
                                    prop.chisq = FALSE, chisq = TRUE)
  tab <- cross_table$t
  pct <- 100 * cross_table$prop.col
  # CrossTable prints both versions. The assignment requires the uncorrected one.
  test <- cross_table$chisq   # Do not use cross_table$chisq.corr (Yates).
  data.frame(
    variable = v, description = unname(variable_labels[v]), reported_value = 1,
    non_dui_n = unname(tab["1", "0"]), non_dui_pct = unname(pct["1", "0"]),
    dui_n = unname(tab["1", "1"]), dui_pct = unname(pct["1", "1"]),
    total_n = sum(tab["1", ]), total_pct = 100 * sum(tab["1", ]) / sum(tab),
    non_dui_denominator = sum(tab[, "0"]), dui_denominator = sum(tab[, "1"]),
    analysis_n = sum(tab), chi_square = unname(test$statistic),
    df = unname(test$parameter), p_value = test$p.value,
    p_display = format.pval(test$p.value, digits = 4, eps = 2.2e-16),
    significant_05 = test$p.value < 0.05,
    conclusion = if (test$p.value < 0.05) "Reject H0: associated with DUI status"
                 else "Do not reject H0: insufficient evidence of association"
  )
})
binary_table <- do.call(rbind, binary_results)
write_table(binary_table, "table_02_crosstab_chisq.csv")

# 4. Continuous predictors: Welch two-sample t-tests (R's default).
# Positive t means the mean in DRINKING_D=0 is higher than in DRINKING_D=1.
continuous_results <- lapply(continuous_vars, function(v) {
  d <- mydata[complete.cases(mydata[c(v, "DRINKING_D")]), c(v, "DRINKING_D")]
  means <- tapply(d[[v]], d$DRINKING_D, mean)
  sds <- tapply(d[[v]], d$DRINKING_D, sd)
  ns <- table(d$DRINKING_D)
  test <- t.test(reformulate("DRINKING_D", response = v), data = d)
  cat("\nContinuous predictor:", v, "\nGroup means:\n")
  print(means)
  cat("Group standard deviations:\n")
  print(sds)
  print(test)
  data.frame(
    variable = v, description = unname(variable_labels[v]),
    unit = if (v == "PCTBACHMOR") "percent" else "USD (2000 Census variable)",
    non_dui_n = unname(ns["0"]),
    non_dui_mean = unname(means["0"]), non_dui_sd = unname(sds["0"]),
    dui_n = unname(ns["1"]), dui_mean = unname(means["1"]), dui_sd = unname(sds["1"]),
    mean_difference_0_minus_1 = unname(means["0"] - means["1"]),
    t = unname(test$statistic), df = unname(test$parameter), p_value = test$p.value,
    p_display = format.pval(test$p.value, digits = 4, eps = 2.2e-16),
    significant_05 = test$p.value < 0.05,
    conclusion = if (test$p.value < 0.05) "Reject H0: group means differ"
                 else "Do not reject H0: insufficient evidence of a difference"
  )
})
continuous_table <- do.call(rbind, continuous_results)
write_table(continuous_table, "table_03_continuous_ttest.csv")

# 5. Pearson correlation matrix of the nine predictors (not the outcome).
# Use the same complete rows for every matrix entry.
cor_data <- mydata[complete.cases(mydata[predictors]), predictors]
cor_matrix <- cor(cor_data, method = "pearson")
cat("\nPearson correlation matrix:\n")
print(round(cor_matrix, 6))
cor_table <- data.frame(variable = rownames(cor_matrix), cor_matrix, check.names = FALSE)
write_table(cor_table, "table_04_correlation_matrix.csv")
pair_indices <- which(upper.tri(cor_matrix), arr.ind = TRUE)
strongest_pair <- pair_indices[which.max(abs(cor_matrix[pair_indices])), ]
strongest_note <- sprintf("The largest absolute correlation is between %s and %s (r = %.4f).",
                          rownames(cor_matrix)[strongest_pair[1]],
                          colnames(cor_matrix)[strongest_pair[2]],
                          cor_matrix[strongest_pair[1], strongest_pair[2]])
high_pairs <- which(upper.tri(cor_matrix) & abs(cor_matrix) > 0.8, arr.ind = TRUE)
if (!nrow(high_pairs)) {
  collinearity_note <- "No predictor pair has |r| > 0.8; this screen does not flag severe pairwise collinearity."
} else {
  pair_notes <- apply(high_pairs, 1, function(ij) {
    sprintf("%s / %s: r = %.4f", rownames(cor_matrix)[ij[1]],
            colnames(cor_matrix)[ij[2]], cor_matrix[ij[1], ij[2]])
  })
  collinearity_note <- paste("Pairs with |r| > 0.8:", paste(pair_notes, collapse = "; "))
}

# A short handoff to Tasks 2 and 3; these are exploratory associations, not causal effects.
p_label <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.4f", p))
notes <- c(
  "HW3 Task 1: exploratory analysis handoff",
  paste("Source CSV:", basename(data_path)),
  paste("Rows in source data:", nrow(mydata)),
  paste("Rows with missing DRINKING_D:", sum(is.na(mydata$DRINKING_D))),
  paste("Missing values across the ten analysis variables:", sum(is.na(mydata[required_vars]))),
  sprintf("DUI: %d of %d non-missing outcomes (%.2f%%).", dv_table$n[2],
          sum(dv_table$n), dv_table$percent[2]),
  "",
  "Binary tables: predictor=1 counts; percentages are within DRINKING_D columns.",
  "Total percentages use all complete cases for that predictor and outcome.",
  "Tests: Pearson chi-square WITHOUT Yates correction; alpha = 0.05.",
  paste(binary_table$variable, binary_table$conclusion,
        p_label(binary_table$p_value), sep = " | "),
  "",
  "Continuous variables: two-sided Welch t-tests; difference = non-DUI mean minus DUI mean.",
  paste(continuous_table$variable, continuous_table$conclusion,
        p_label(continuous_table$p_value), sep = " | "),
  "",
  paste("Complete cases in the nine-predictor correlation matrix:", nrow(cor_data)),
  strongest_note,
  collinearity_note,
  "For two binary variables, Pearson r equals the phi coefficient and is affected by category prevalence.",
  "This correlation screen does not prove that all regression assumptions hold.",
  "",
  "Missing values: each univariate test excludes missing predictor/outcome pairs.",
  "The original data are not edited; continuous predictors remain in their original units.",
  "PCTBACHMOR and MEDHHINC describe the crash location's 2000 Census block group, not the individual driver.",
  "Their group means are computed across crash records, so block groups are weighted by their observed crash counts.",
  "Task 2 should start from the same original CSV. Task 1 does not fit a logistic model.",
  "Keep all nine specified predictors for the assigned full model; do not select predictors using these bivariate p-values.",
  paste("CrossTable output generated with gmodels", as.character(packageVersion("gmodels"))),
  "Binary CSV counts, column percentages and chi-square tests come directly from CrossTable return values.",
  "The exported test is CrossTable$chisq (uncorrected), not CrossTable$chisq.corr (Yates).",
  "p_display is formatted for reading; p_value retains the numerical result.",
  paste("R version:", R.version.string)
)
writeLines(notes, file.path("output", "task1_handoff.txt"))

# 6. Readable tables and brief interpretations for the group handoff.
# These use the same unrounded results exported to the four CSV files above.
number <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")
display_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
markdown_table <- function(headers, rows) {
  c(paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |"),
    apply(as.matrix(rows), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
binary_interpretations <- vapply(seq_len(nrow(binary_table)), function(i) {
  r <- binary_table[i, ]
  direction <- if (r$dui_pct > r$non_dui_pct) "more" else "less"
  conclusion <- if (r$significant_05)
    paste0("Reject H0; the characteristic was ", direction, " common in alcohol-involved crashes.")
    else "Do not reject H0; the test does not show a significant association."
  sprintf("- **%s:** %s%% of alcohol-involved versus %s%% of non-alcohol crashes; %s. %s",
          r$variable, number(r$dui_pct), number(r$non_dui_pct), p_label(r$p_value), conclusion)
}, character(1))
continuous_interpretations <- vapply(seq_len(nrow(continuous_table)), function(i) {
  r <- continuous_table[i, ]
  conclusion <- if (r$significant_05) "Reject H0: the group means differ."
    else "Do not reject H0: there is insufficient evidence that the group means differ."
  sprintf("- **%s:** t(%s) = %s, %s. %s", r$variable, number(r$df),
          number(r$t, 3), p_label(r$p_value), conclusion)
}, character(1))
results <- c(
  "# HW3 Task 1: Exploratory Analysis", "",
  "## Data and methods", "",
  sprintf("The supplied CSV contains %s Philadelphia crash records. The assignment describes crashes in residential block groups during 2008-2012, with neighborhood characteristics from the 2000 Census.", number(nrow(mydata), 0)),
  sprintf("There are %s missing values across DRINKING_D and the nine predictors. Tests use alpha = 0.05.", number(sum(is.na(mydata[required_vars])), 0)),
  "Binary predictors are compared with Pearson chi-square tests without Yates correction. Continuous predictors are compared with two-sided Welch t-tests. Pearson correlations describe the nine predictors, with |r| > 0.8 used as the task's screening threshold.", "",
  "## Table 1. Drinking-driver indicator", "",
  markdown_table(c("DRINKING_D", "Group", "N", "%"),
                 data.frame(dv_table$DRINKING_D, dv_table$group, number(dv_table$n, 0), number(dv_table$percent))), "",
  sprintf("Alcohol-involved crashes account for %s of %s records (%.2f%%).",
          number(dv_table$n[2], 0), number(sum(dv_table$n), 0), dv_table$percent[2]), "",
  "## Table 2. Binary predictors and chi-square tests", "",
  "Each row reports the predictor=1 category. No alcohol = DRINKING_D=0; alcohol = DRINKING_D=1. Percentages are within each drinking-status column, so they do not add to 100% across a row.", "",
  markdown_table(c("Variable", "No alcohol N", "No alcohol %", "Alcohol N", "Alcohol %", "Total N", "Chi-square", "df", "p"),
                 data.frame(binary_table$variable, number(binary_table$non_dui_n, 0), number(binary_table$non_dui_pct),
                            number(binary_table$dui_n, 0), number(binary_table$dui_pct), number(binary_table$total_n, 0),
                            number(binary_table$chi_square, 3), binary_table$df, display_p(binary_table$p_value))), "",
  "Variable definitions:", "",
  paste0("- **", binary_vars, ":** ", unname(variable_labels[binary_vars]), "."), "",
  "Test conclusions:", "", binary_interpretations, "",
  "## Table 3. Continuous predictors and Welch t-tests", "",
  markdown_table(c("Variable", "No alcohol mean", "No alcohol SD", "Alcohol mean", "Alcohol SD", "t", "df", "p"),
                 data.frame(continuous_table$variable, number(continuous_table$non_dui_mean), number(continuous_table$non_dui_sd),
                            number(continuous_table$dui_mean), number(continuous_table$dui_sd), number(continuous_table$t, 3),
                            number(continuous_table$df), display_p(continuous_table$p_value))), "",
  "PCTBACHMOR is the percentage of adults aged 25+ with at least a bachelor's degree in the crash location's block group. Its mean and SD are on the 0-100 percentage scale. MEDHHINC is that block group's median household income, in the original Census dollar units.",
  "These are neighborhood characteristics, not the driver's own education or income. Means are calculated across crashes, so they are weighted by the number of crash records in each block group. The t-statistic uses mean(DRINKING_D=0) minus mean(DRINKING_D=1).", "",
  continuous_interpretations, "",
  "## Table 4. Pearson correlations among predictors", "",
  markdown_table(c("Variable", predictors), data.frame(variable = predictors,
    apply(cor_matrix, 2, function(x) sprintf("%.3f", x)), check.names = FALSE)), "",
  strongest_note, collinearity_note,
  "For two binary predictors, Pearson r is the phi coefficient. Its magnitude depends on the category proportions, so a low r does not by itself prove there is no relationship. This matrix is a pairwise collinearity screen, not a test of every regression assumption.", "",
  "## Handoff", "",
  "The CSV tables retain numerical precision; these display tables round counts and statistics for reading. The script does not change the source data. All associations above are unadjusted; they should not be interpreted as causal effects or as substitutes for the assigned multivariable logistic models. Task 2 should use all nine specified predictors in the full model.",
  paste0("Cross-tabulations were generated with gmodels::CrossTable (gmodels ",
         as.character(packageVersion("gmodels")), "). Counts and column percentages were exported directly from its returned tables. The reported chi-square statistic, degrees of freedom and p-value come from the uncorrected chisq result, not the Yates-corrected chisq.corr result."), "",
  "## Sources", "",
  "- Logistic Regression Data.csv, provided with the assignment; one row per crash.",
  "- HW 3 - Instructions.pdf, pp. 1-2 (data definitions), pp. 3-6 (exploratory analyses), and p. 11 (interpretation requirements).",
  "- HW3_分工方案(1).md, Task 1 (deliverables and the |r| > 0.8 threshold)."
)
writeLines(results, file.path("output", "task1_results.md"), useBytes = TRUE)
cat("\n", paste(notes, collapse = "\n"), "\n", sep = "")
cat("\nFour CSV tables, task1_results.md and task1_handoff.txt written to:", file.path(getwd(), "output"), "\n")
