Homework 3 — Logistic Regression: Predictors of Alcohol-Involved Crashes
Philadelphia crash records, 2008-2012 (43,364 crashes in residential block groups)

Authors: Qiwen Bian, Manching Ho, and Fujia Zheng.

Reading the results
-------------------
HW3_report.pdf is the final report. HW3_LaTeX.zip contains its LaTeX
source and the ROC figure. The output/ folder contains the R tables,
the ROC figure, and the handoff notes from the exploratory (Task 1)
and modeling (Task 2) stages.

All models estimate the drinking-driver indicator from nine predictors:

    DRINKING_D ~ FATAL_OR_M + OVERTURNED + CELL_PHONE + SPEEDING
               + AGGRESSIVE + DRIVER1617 + DRIVER65PLUS
               + PCTBACHMOR + MEDHHINC

Main results
------------
Alcohol was involved in 2,485 of 43,364 crashes (5.73%). Six of the
seven crash characteristics are associated with alcohol involvement in
chi-square tests. In the full logistic model, speeding (OR = 4.66),
overturning (OR = 2.53), and fatality or major injury (OR = 2.26) raise
the odds of alcohol involvement, while aggressive driving (OR = 0.55),
a driver aged 16-17 (OR = 0.28), and a driver aged 65+ (OR = 0.46)
lower them. Cell-phone use and block-group education are not
significant; block-group median household income is barely significant
(OR = 1.003 per $1,000). The ROC AUC is 0.640 (poor discrimination).
The cut-off minimizing the distance to the upper-left corner of the
ROC curve is 0.064; the cut-off minimizing the misclassification rate
(0.50, 5.73% errors) reflects the base rate rather than discrimination.
Dropping PCTBACHMOR and MEDHHINC changes no coefficient's significance
and yields nearly the same AIC (18,360.5 vs. 18,359.6).

Reproduce
---------
  cd code
  Rscript 01_exploratory.R
  Rscript 02_logistic.R

Compile the report: unzip HW3_LaTeX.zip, then run
  tectonic HW3_report.tex

Data: "Logistic Regression Data.csv", one row per crash, provided with
the assignment (PennDOT crashes spatially joined to Census 2000 block
groups).
