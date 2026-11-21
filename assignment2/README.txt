Homework 2 — Spatial Autocorrelation and Spatial Regression
Philadelphia median home values in 1,720 Census block groups

Authors: Qiwen Bian, Manching Ho, and Fujia Zheng.
Fujia Zheng completed the geographically weighted regression (GWR) and the
comparison of all four models using the group's shared OLS and spatial
regression work.

Reading the results
-------------------
HW2_report.pdf is the final report. HW2_LaTeX.zip contains its LaTeX
source and all required figures. The output/ folder contains the R
tables, figures, and diagnostics used to check the results. Table 21
compares the four models under both AIC conventions.

All models estimate log median home value from the same four predictors:

    LNMEDHVAL ~ PCTVACANT + PCTSINGLES + PCTBACHMOR + LNNBELPOV

LNMEDHVAL is log(MEDHVAL + 1). LNNBELPOV is log(NBelPov100 + 1). The first
three predictors are percent vacant housing, percent detached single-family housing, and percent
with a bachelor's degree or higher, as coded in the course data.

Main results
------------
Log median home value clusters geographically (global Moran's I = 0.794,
one-sided permutation p = 0.001). OLS explains 66.2% of the sample
variation. Its residuals retain positive spatial autocorrelation (I = 0.313).

The spatial lag estimate is rho = 0.651 and the spatial error estimate is
lambda = 0.815. Their residual Moran's I values are -0.082 and -0.095,
respectively (two-sided p = 0.002 for both). All residual comparisons use
the same row-standardized Queen weights and 999 permutations, with seed 671.

GWR uses an adaptive Gaussian kernel. At the lowest AICc among the tested
bandwidths (662.25), the bandwidth is the distance to the 14th other
centroid, a median radius of 723 metres. Observations farther away still
receive smaller positive Gaussian weights. The in-sample quasi-global R²
is 0.848. GWR has the residual Moran's I closest to zero (0.034), with
two-sided p = 0.020. Some positive spatial dependence remains.

Under the common likelihood-based AIC comparison, the spatial lag model
has the lowest AIC (525.12), followed by GWR (570.16), spatial error
(758.99), and OLS (1,434.99). The GWR penalty here uses twice its effective
parameter count, 2 × [trace(S) + 1], where S is the GWR smoothing matrix.
The course's spgwr equation 4.22 convention gives 310.36 for the same GWR
fit, using -2 × log likelihood + trace(S); under that convention GWR ranks
first. Table 18 records both GWR AIC values and the bandwidth-selection
AICc. The report and table 21 explain the cross-model ranking.

Files and sources
-----------------
code/HW2_analysis.Rmd contains the narrative, code and figures for Parts
1–3. The original lecture reference is code/Lecture 14 模版.Rmd.

data/RegressionData.csv and data/map_RegressionData/ are the course inputs.
The CSV is aligned to the source shapefile by POLY_ID. The shared
data/RegressionData_queen.gal preserves the neighbor lists from the supplied
first-order Queen GAL file. This graph defines the spatial regressions and
the global, local and residual Moran diagnostics.

The course geometry has no documented native coordinate system. GWR uses
official U.S. Census Bureau Philadelphia County Census 2000 block groups
from the 2010 TIGER/Line release. Each course AREAKEY matches one official
BKGPIDFP00. The matched polygons are projected to NAD83 / UTM zone 18N
(EPSG:26918) before calculating centroid distances in metres. The original
archive, retrieval URL, checksum and detailed provenance are under
data/crs_reference/. data/gwr_reference_centroids.csv holds the 1,720
matched centroids. The corresponding polygons are stored in
data/crs_reference/Philadelphia_analysis_bg.gpkg.

Part 1: code/HW2_OLS_analysis.R generates figures 01–09, tables 01–08 and
data/processed_ols/ (fitted values, residuals, neighbor-averaged residuals
and local Moran classifications). Global Moran tests in this part use the
positive alternative. Local classifications use 999 conditional
permutations without replacement and a two-sided 0.05 threshold.

Part 2: code/HW2_SL_SE_analysis.R generates figures 10–15, tables 09–16 and
data/processed_part2/ (the Part 1 fields, spatial predictions and residuals).
Table 09 contains the OLS score diagnostics for lag and error alternatives.
The full OLS, lag and error summaries are also available as text in output/.

Part 3: code/HW2_prepare_geography.R rebuilds the projected Census subset
and centroids. code/HW2_GWR_analysis.R generates figures 16–19 and tables
17–20 for bandwidth selection, fit, local coefficients, local R² and
residual spatial autocorrelation. data/processed_gwr/GWR_local_results.csv
contains one row per block group, including local estimates, standard errors,
coefficient-to-standard-error ratios, local R², fitted values and residuals.
Use POLY_ID and AREAKEY to match these results to the Census polygons.

Reproducing the analysis
------------------------
Run the following command from the package directory that contains this
README, keeping the code/, data/ and output/ folders together:

    Rscript code/render_report.R

The analysis was rerun with R 4.5.3. The analysis requires sf, spdep, spatialreg,
lmtest, ggplot2 and patchwork. Rendering additionally uses knitr, rmarkdown
and Pandoc to regenerate the companion HTML and Markdown. After
installing any missing dependencies, the included data
allow the analysis to run offline. The three output/session_info*.txt files
record the R and package versions.
