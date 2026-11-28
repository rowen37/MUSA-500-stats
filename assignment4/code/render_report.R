# From repository root: Rscript assignment4/code/render_report.R
# From assignment4/:    Rscript code/render_report.R
# Add --all to regenerate A and B as well; default preserves their checked-in outputs.
args <- commandArgs(trailingOnly=TRUE)
file_arg <- grep('^--file=',commandArgs(),value=TRUE)
if(!length(file_arg)) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(sub('^--file=','',file_arg[1])),'..'))
setwd(root)
required <- c('sf','spatstat.geom','spatstat.explore','spatstat.random','ggplot2','rmarkdown','knitr')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install packages first: install.packages(c(',paste(sprintf('"%s"',missing),collapse=','),'))')
suppressPackageStartupMessages(library(sf))
# RStudio bundles Pandoc; this does not require a separate installation on macOS.
if(!rmarkdown::pandoc_available()) {
  candidates <- c('/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64',
    '/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64')
  pd <- candidates[file.exists(file.path(candidates,'pandoc'))]
  if(length(pd)) Sys.setenv(RSTUDIO_PANDOC=pd[1])
}
if('--all' %in% args) {
  source('code/HW4_prepare_data.R'); prepare_hw4_data('.')
  source('code/HW4_NNA_analysis.R')
  source('code/HW4_Kfunction_analysis.R')
}
if(!('--report-only' %in% args)) {
  source('code/HW4_Kinhom_analysis.R'); run_hw4_part_c('.')
}
stopifnot(file.exists('output/table_01_nna_summary.csv'),file.exists('output/table_02_kfunction_bands.csv'))
rmarkdown::render('code/HW4_analysis.Rmd',output_format='html_document',
  output_file='HW4_report.html',output_dir=file.path(root,'output'),quiet=TRUE,
  envir=new.env(parent=globalenv()))
if(!('--html-only' %in% args)) {
  engine <- if(nzchar(Sys.which('tectonic'))) 'tectonic' else 'pdflatex'
  if(!nzchar(Sys.which(engine))) stop('HTML completed. For PDF, install TinyTeX with tinytex::install_tinytex(), or install Tectonic; then rerun.')
  if(engine=='tectonic') {
    rmarkdown::render('code/HW4_analysis.Rmd',output_format=rmarkdown::latex_document(),
      output_file='HW4_report.tex',output_dir=file.path(root,'output'),quiet=TRUE,
      envir=new.env(parent=globalenv()))
    status <- system2(Sys.which('tectonic'),c('-k',shQuote(file.path(root,'output','HW4_report.tex'))))
    if(status!=0L) stop('Tectonic failed; inspect the compiler output.')
  } else {
    rmarkdown::render('code/HW4_analysis.Rmd',output_format=rmarkdown::pdf_document(
      latex_engine=engine,keep_tex=TRUE,number_sections=FALSE),
      output_file='HW4_report.pdf',output_dir=file.path(root,'output'),quiet=TRUE,
      envir=new.env(parent=globalenv()))
  }
  # Package the generated LaTeX and referenced figures with portable paths.
  bundle <- tempfile('HW4_LaTeX_'); dir.create(bundle)
  dir.create(file.path(bundle,'figures'))
  tex_path <- 'output/HW4_report.tex'
  tex <- readLines(tex_path,warn=FALSE)
  figs <- c('fig_04_kfunction_envelope.png','fig_07_lambda_intensity_map.png',
            'fig_08_kinhom_curves.png','fig_09_income_map_markets.png')
  for(f in figs) {
    file.copy(file.path('output',f),file.path(bundle,'figures',f))
    tex <- gsub(paste0('../output/',tools::file_path_sans_ext(f)),paste0('figures/',tools::file_path_sans_ext(f)),tex,fixed=TRUE)
    tex <- gsub(paste0(root,'/output/',tools::file_path_sans_ext(f)),paste0('figures/',tools::file_path_sans_ext(f)),tex,fixed=TRUE)
  }
  writeLines(tex,file.path(bundle,'HW4_report.tex'))
  writeLines('Compile HW4_report.tex with Tectonic or pdfLaTeX (two passes). Figures are included.',file.path(bundle,'README.txt'))
  old <- getwd(); setwd(bundle)
  archive <- file.path(root,'output','HW4_LaTeX.zip')
  if(file.exists(archive)) unlink(archive)
  utils::zip(archive,list.files('.',recursive=TRUE),flags='-q')
  setwd(old)
}
writeLines(capture.output(sessionInfo()),'output/session_info_hw4.txt')
cat('Report complete in assignment4/output/.\n')
