# HW3 Task 1: Exploratory Analysis

## Data and methods

The supplied CSV contains 43,364 Philadelphia crash records. The assignment describes crashes in residential block groups during 2008-2012, with neighborhood characteristics from the 2000 Census.
There are 0 missing values across DRINKING_D and the nine predictors. Tests use alpha = 0.05.
Binary predictors are compared with Pearson chi-square tests without Yates correction. Continuous predictors are compared with two-sided Welch t-tests. Pearson correlations describe the nine predictors, with |r| > 0.8 used as the task's screening threshold.

## Table 1. Drinking-driver indicator

| DRINKING_D | Group | N | % |
| --- | --- | --- | --- |
| 0 | Non-DUI | 40,879 | 94.27 |
| 1 | DUI | 2,485 | 5.73 |

Alcohol-involved crashes account for 2,485 of 43,364 records (5.73%).

## Table 2. Binary predictors and chi-square tests

Each row reports the predictor=1 category. No alcohol = DRINKING_D=0; alcohol = DRINKING_D=1. Percentages are within each drinking-status column, so they do not add to 100% across a row.

| Variable | No alcohol N | No alcohol % | Alcohol N | Alcohol % | Total N | Chi-square | df | p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FATAL_OR_M | 1,181 | 2.89 | 188 | 7.57 | 1,369 | 167.562 | 1 | <0.001 |
| OVERTURNED | 612 | 1.50 | 110 | 4.43 | 722 | 122.788 | 1 | <0.001 |
| CELL_PHONE | 426 | 1.04 | 28 | 1.13 | 454 | 0.162 | 1 | 0.6873 |
| SPEEDING | 1,261 | 3.08 | 260 | 10.46 | 1,521 | 376.781 | 1 | <0.001 |
| AGGRESSIVE | 18,522 | 45.31 | 916 | 36.86 | 19,438 | 67.602 | 1 | <0.001 |
| DRIVER1617 | 674 | 1.65 | 12 | 0.48 | 686 | 20.452 | 1 | <0.001 |
| DRIVER65PLUS | 4,237 | 10.36 | 119 | 4.79 | 4,356 | 80.605 | 1 | <0.001 |

Variable definitions:

- **FATAL_OR_M:** Fatality or major injury.
- **OVERTURNED:** Overturned vehicle.
- **CELL_PHONE:** Driver using a cell phone.
- **SPEEDING:** Speeding vehicle.
- **AGGRESSIVE:** Aggressive driving.
- **DRIVER1617:** At least one driver aged 16 or 17.
- **DRIVER65PLUS:** At least one driver aged 65 or older.

Test conclusions:

- **FATAL_OR_M:** 7.57% of alcohol-involved versus 2.89% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was more common in alcohol-involved crashes.
- **OVERTURNED:** 4.43% of alcohol-involved versus 1.50% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was more common in alcohol-involved crashes.
- **CELL_PHONE:** 1.13% of alcohol-involved versus 1.04% of non-alcohol crashes; p = 0.6873. Do not reject H0; the test does not show a significant association.
- **SPEEDING:** 10.46% of alcohol-involved versus 3.08% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was more common in alcohol-involved crashes.
- **AGGRESSIVE:** 36.86% of alcohol-involved versus 45.31% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was less common in alcohol-involved crashes.
- **DRIVER1617:** 0.48% of alcohol-involved versus 1.65% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was less common in alcohol-involved crashes.
- **DRIVER65PLUS:** 4.79% of alcohol-involved versus 10.36% of non-alcohol crashes; p < 0.001. Reject H0; the characteristic was less common in alcohol-involved crashes.

## Table 3. Continuous predictors and Welch t-tests

| Variable | No alcohol mean | No alcohol SD | Alcohol mean | Alcohol SD | t | df | p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PCTBACHMOR | 16.57 | 18.21 | 16.61 | 18.72 | -0.108 | 2,777.54 | 0.9137 |
| MEDHHINC | 31,483.05 | 16,930.10 | 31,998.75 | 17,810.50 | -1.405 | 2,763.87 | 0.1600 |

PCTBACHMOR is the percentage of adults aged 25+ with at least a bachelor's degree in the crash location's block group. Its mean and SD are on the 0-100 percentage scale. MEDHHINC is that block group's median household income, in the original Census dollar units.
These are neighborhood characteristics, not the driver's own education or income. Means are calculated across crashes, so they are weighted by the number of crash records in each block group. The t-statistic uses mean(DRINKING_D=0) minus mean(DRINKING_D=1).

- **PCTBACHMOR:** t(2,777.54) = -0.108, p = 0.9137. Do not reject H0: there is insufficient evidence that the group means differ.
- **MEDHHINC:** t(2,763.87) = -1.405, p = 0.1600. Do not reject H0: there is insufficient evidence that the group means differ.

## Table 4. Pearson correlations among predictors

| Variable | FATAL_OR_M | OVERTURNED | CELL_PHONE | SPEEDING | AGGRESSIVE | DRIVER1617 | DRIVER65PLUS | PCTBACHMOR | MEDHHINC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FATAL_OR_M | 1.000 | 0.033 | 0.002 | 0.082 | -0.011 | -0.003 | -0.013 | -0.015 | -0.018 |
| OVERTURNED | 0.033 | 1.000 | -0.001 | 0.059 | 0.016 | 0.004 | -0.020 | 0.009 | 0.028 |
| CELL_PHONE | 0.002 | -0.001 | 1.000 | -0.004 | -0.026 | 0.001 | -0.003 | -0.001 | 0.002 |
| SPEEDING | 0.082 | 0.059 | -0.004 | 1.000 | 0.212 | 0.016 | -0.033 | -0.001 | 0.012 |
| AGGRESSIVE | -0.011 | 0.016 | -0.026 | 0.212 | 1.000 | 0.028 | 0.015 | 0.027 | 0.043 |
| DRIVER1617 | -0.003 | 0.004 | 0.001 | 0.016 | 0.028 | 1.000 | -0.021 | -0.003 | 0.023 |
| DRIVER65PLUS | -0.013 | -0.020 | -0.003 | -0.033 | 0.015 | -0.021 | 1.000 | 0.026 | 0.050 |
| PCTBACHMOR | -0.015 | 0.009 | -0.001 | -0.001 | 0.027 | -0.003 | 0.026 | 1.000 | 0.478 |
| MEDHHINC | -0.018 | 0.028 | 0.002 | 0.012 | 0.043 | 0.023 | 0.050 | 0.478 | 1.000 |

The largest absolute correlation is between PCTBACHMOR and MEDHHINC (r = 0.4779).
No predictor pair has |r| > 0.8; this screen does not flag severe pairwise collinearity.
For two binary predictors, Pearson r is the phi coefficient. Its magnitude depends on the category proportions, so a low r does not by itself prove there is no relationship. This matrix is a pairwise collinearity screen, not a test of every regression assumption.

## Handoff

The CSV tables retain numerical precision; these display tables round counts and statistics for reading. The script does not change the source data. All associations above are unadjusted; they should not be interpreted as causal effects or as substitutes for the assigned multivariable logistic models. Task 2 should use all nine specified predictors in the full model.
Cross-tabulations were generated with gmodels::CrossTable (gmodels 2.19.1). Counts and column percentages were exported directly from its returned tables. The reported chi-square statistic, degrees of freedom and p-value come from the uncorrected chisq result, not the Yates-corrected chisq.corr result.

## Sources

- Logistic Regression Data.csv, provided with the assignment; one row per crash.
- HW 3 - Instructions.pdf, pp. 1-2 (data definitions), pp. 3-6 (exploratory analyses), and p. 11 (interpretation requirements).
- HW3_分工方案(1).md, Task 1 (deliverables and the |r| > 0.8 threshold).
