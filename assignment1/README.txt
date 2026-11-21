Homework 1 submission package

Authors: Qiwen Bian, Manching Ho, Fujia Zheng

Contents
- code/HW1_analysis.R: complete reproducible analysis and figure generation
- data/: regression CSV and shapefile components
- output/HW1_report.tex: LaTeX report source
- output/HW1_report.pdf: compiled report
- output/fig_*.png: report figures
- output/*.csv: computed tables and diagnostics

Reproduce the figures and tables:
  cd code
  Rscript HW1_analysis.R

Compile the report:
  cd output
  tectonic HW1_report.tex

The analysis script uses only local course data and writes results to output/.
