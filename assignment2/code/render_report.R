# Run from the assignment root with Rscript code/render_report.R.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
setwd(root)
packages <- c("rmarkdown", "knitr", "sf", "spdep", "lmtest", "spatialreg", "ggplot2", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the following R packages first: ", paste(missing, collapse = ", "))
if (!rmarkdown::pandoc_available()) {
  for (folder in c("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
                   "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64")) {
    if (file.exists(file.path(folder, "pandoc"))) {
      Sys.setenv(RSTUDIO_PANDOC = folder)
      break
    }
  }
}
if (!rmarkdown::pandoc_available()) stop("Install Pandoc or render from an RStudio session.")
source("code/HW2_prepare_geography.R")
prepare_gwr_geography(root)
dir.create("output", showWarnings = FALSE)
rmarkdown::render("code/HW2_analysis.Rmd", output_dir = file.path(root, "output"),
                  envir = new.env(parent = globalenv()))
# The HTML embeds figures. The companion Markdown uses paths relative to its own folder.
md <- file.path(root, "output", "HW2_analysis.md")
lines <- readLines(md, warn = FALSE)
lines <- gsub(paste0(root, "/output/"), "", lines, fixed = TRUE)
writeLines(lines, md)
