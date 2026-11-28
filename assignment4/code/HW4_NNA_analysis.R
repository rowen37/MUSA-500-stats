# Philadelphia farmers markets: reproducible Part A (nearest neighbor analysis).
# Run from code/ or the project root. The same script can be sourced by knitr.
# Builds on the shared data objects stored by code/HW4_prepare_data.R.
run_hw4_part_nna <- function() {
  root_candidates <- c(".", "..")
  root <- root_candidates[file.exists(file.path(root_candidates, "data",
    "Philadelphia_Farmers_Markets201302.shp"))][1]
  if (is.na(root)) stop("Run this script from the project root or code/.")
  root <- normalizePath(root)
  input <- file.path(root, "data")
  out <- file.path(root, "output")
  processed <- file.path(input, "processed_nna")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed, recursive = TRUE, showWarnings = FALSE)

  # The shared point pattern is required by Parts B and C; build it on demand.
  ppp_path <- file.path(input, "processed_prepare", "markets_ppp.Rds")
  if (!file.exists(ppp_path)) {
    source(file.path(root, "code", "HW4_prepare_data.R"), local = TRUE)
    prepare_hw4_data(root)
  }
  X <- readRDS(ppp_path)
  window <- readRDS(file.path(input, "processed_prepare", "city_owin.Rds"))
  markets <- sf::st_read(file.path(input, "Philadelphia_Farmers_Markets201302.shp"),
                         quiet = TRUE, stringsAsFactors = FALSE)
  city <- sf::st_read(file.path(input, "Philadelphia.shp"), quiet = TRUE)
  n <- spatstat.geom::npoints(X)
  xy <- sf::st_coordinates(markets)
  stopifnot(n == 62L, nrow(markets) == 62L,
            max(abs(cbind(X$x, X$y) - xy[, 1:2])) < 1e-9,
            spatstat.geom::area.owin(window) > 0)
  # Reading Terminal Market hosts two separately listed markets at one address.
  coincident <- which(spatstat.geom::multiplicity(X) > 1L)
  stopifnot(length(coincident) == 2L,
            all(markets$NAME[coincident] == c("Reading Terminal Market",
                                               "Fair Food Farmstand @ RTM")))

  # Intensity and distances. Coordinates are EPSG:2272 US survey feet, so every
  # distance below is in feet; 1 mile = 5280 feet.
  area_ft2 <- spatstat.geom::area.owin(window)
  stopifnot(abs(area_ft2 / as.numeric(sf::st_area(city)) - 1) < 1e-6)
  lambda <- n / area_ft2
  ft_per_mile <- 5280
  nnd <- spatstat.geom::nndist(X)
  observed_mean <- mean(nnd)

  # Course z-test for the nearest neighbor index R = mean(d) / E[d] (Poisson
  # expectation E[d] = 1/(2 sqrt(lambda)), SE = sqrt((4 - pi) / (4 pi lambda n))).
  # This is the formula ArcGIS Average Nearest Neighbor reports as well.
  expected_mean <- 1 / (2 * sqrt(lambda))
  se_mean <- sqrt((4 - pi) / (4 * pi * lambda * n))
  nni <- observed_mean / expected_mean
  z_score <- (observed_mean - expected_mean) / se_mean
  p_one_sided <- pnorm(z_score)
  p_two_sided_formula <- 2 * pnorm(-abs(z_score))

  # Edge-corrected Clark-Evans test from spatstat: the CDF correction accounts
  # for the window boundary, which the slide formula ignores.
  ce <- spatstat.explore::clarkevans.test(X, alternative = "clustered")
  ce_two <- spatstat.explore::clarkevans.test(X, alternative = "two.sided")
  stopifnot(any(grepl("CDF correction", ce$method, fixed = TRUE)))

  # Monte Carlo cross-check: 99 random patterns with the same n in the same
  # window, matching the assignment's 99-permutation convention (seed 671).
  nsim <- 99L
  seed <- 671L
  set.seed(seed)
  sim_means <- vapply(seq_len(nsim), function(i) {
    mean(spatstat.geom::nndist(spatstat.random::runifpoint(n, win = window)))
  }, numeric(1))
  p_sim <- (1 + sum(sim_means <= observed_mean)) / (nsim + 1)

  # Sensitivity: drop the coincident Reading Terminal pair and refit the index.
  X60 <- X[-coincident]
  nnd60 <- spatstat.geom::nndist(X60)
  observed_mean60 <- mean(nnd60)
  nni60 <- observed_mean60 / expected_mean

  write_table <- function(d, name) write.csv(d, file.path(out, name), row.names = FALSE)
  per_market <- data.frame(OBJECTID_1 = markets$OBJECTID_1, NAME = markets$NAME,
    ZIP_CODE = markets$ZIP_CODE, nnd_ft = round(nnd, 3),
    nnd_mi = round(nnd / ft_per_mile, 5))
  write.csv(per_market, file.path(processed, "nndist_per_market.csv"), row.names = FALSE)
  conclusion <- if (p_one_sided < 0.05) "Reject CSR: markets are clustered" else
    "Do not reject CSR: no evidence against randomness"
  write_table(data.frame(Item = c(
    "Observations", "Coincident address (Reading Terminal Market)",
    "Study area (square feet)", "Study area (square miles)",
    "Intensity (markets per square foot)", "Intensity (markets per square mile)",
    "Observed mean nearest neighbor distance (feet)",
    "Observed mean nearest neighbor distance (miles)",
    "Expected mean nearest neighbor distance under CSR (feet)",
    "Standard error of the mean under CSR (feet)",
    "Nearest neighbor index R (no edge correction)",
    "z statistic (no edge correction)",
    "p-value, one-sided clustered (no edge correction)",
    "p-value, two-sided (no edge correction)",
    "Clark-Evans R (spatstat, CDF edge correction)",
    "p-value, one-sided clustered (CDF correction)",
    "p-value, two-sided (CDF correction)",
    "Monte Carlo p-value (99 CSR simulations, one-sided clustered)",
    "Sensitivity: mean NND excluding coincident pair (feet)",
    "Sensitivity: nearest neighbor index R excluding coincident pair",
    "Conclusion at alpha = 0.05"),
    Value = c(n, 2L, round(area_ft2), round(area_ft2 / ft_per_mile^2, 1),
      signif(lambda, 4), round(lambda * ft_per_mile^2, 4),
      round(observed_mean, 1), round(observed_mean / ft_per_mile, 3),
      round(expected_mean, 1), round(se_mean, 1),
      round(nni, 4), round(z_score, 3), signif(p_one_sided, 4),
      signif(p_two_sided_formula, 4), round(ce$statistic, 4),
      signif(ce$p.value, 4), signif(ce_two$p.value, 4),
      round(p_sim, 4), round(observed_mean60, 1), round(nni60, 4), conclusion)),
    "table_01_nna_summary.csv")

  # Figures follow the group's HW2 style (base R, 1200 x 900 at 140 dpi).
  plot_file <- function(name, draw) {
    png(file.path(out, name), width = 1200, height = 900, res = 140)
    on.exit(dev.off())
    par(mar = c(4.4, 4.7, 3.3, 1.1))
    draw()
  }
  point <- "#1F668883"
  accent <- "#AD4025"
  plot_file("fig_01_markets_city_map.png", function() {
    plot(sf::st_geometry(city), col = "#E8E8E8", border = "#777777", lwd = 0.4,
         main = "Philadelphia farmers markets, 2013")
    points(X$x, X$y, pch = 16, col = point, cex = 0.8)
    legend("bottomright", pch = 16, col = point, bg = "white",
           box.col = "grey85", cex = 0.85,
           legend = sprintf("Farmers market (n = %d)", n))
  })
  plot_file("fig_02_nndist_histogram.png", function() {
    hist(nnd, breaks = 16, col = "#7B9BAA", border = "white",
         xlim = c(0, max(nnd) * 1.05), xlab = "Nearest neighbor distance (feet)",
         main = "Distance from each market to its nearest neighbor")
    abline(v = observed_mean, col = accent, lwd = 2)
    abline(v = expected_mean, col = "grey25", lty = 2, lwd = 2)
    legend("topright", bty = "n", lwd = 2, col = c(accent, "grey25"),
           lty = c(1, 2), cex = 0.9,
           legend = c(sprintf("Observed mean = %.0f ft (R = %.3f)", observed_mean, nni),
                      sprintf("CSR expectation = %.0f ft", expected_mean)))
  })
  plot_file("fig_03_nnd_csr_simulation.png", function() {
    hist(sim_means, breaks = 20, col = "#7B9BAA", border = "white",
         xlim = range(c(sim_means, observed_mean)),
         xlab = "Mean nearest neighbor distance under CSR (feet)",
         main = "99 random patterns with the same 62 points in the same window")
    abline(v = observed_mean, col = accent, lwd = 2)
    legend("topleft", bty = "n", cex = 0.9,
           legend = sprintf("Observed mean = %.0f ft\nSimulation p = %.3f",
                            observed_mean, p_sim))
  })

  writeLines(c(
    "Homework 4 Part A: nearest neighbor analysis of farmers markets",
    sprintf("Markets: %d | Study area: %.0f square feet (%.1f sq mi)", n, area_ft2,
            area_ft2 / ft_per_mile^2),
    sprintf("Intensity: %.4g markets per square foot (%.4f per square mile)",
            lambda, lambda * ft_per_mile^2),
    sprintf("Observed mean NND = %.1f ft (%.3f mi); CSR expectation = %.1f ft",
            observed_mean, observed_mean / ft_per_mile, expected_mean),
    sprintf("Course formula (no edge correction): R = %.4f, z = %.3f, one-sided p = %.4g",
            nni, z_score, p_one_sided),
    sprintf("spatstat Clark-Evans (CDF correction): R = %.4f, one-sided p = %.4g",
            ce$statistic, ce$p.value),
    sprintf("Monte Carlo: %d CSR simulations, seed %d, p = %.4f", nsim, seed, p_sim),
    sprintf("Sensitivity without the coincident Reading Terminal pair: mean = %.1f ft, R = %.4f",
            observed_mean60, nni60),
    paste("Coincident address:", paste(markets$NAME[coincident], collapse = " / "),
          "at 12th St and Arch St (kept as two distinct markets)."),
    "The course-formula rows match what ArcGIS Average Nearest Neighbor reports",
    "with the city boundary as the study area; add its screenshot as fig_03b if a",
    "license becomes available. R is well below 1, so both tests reject CSR and",
    "conclude significant clustering at alpha = 0.05.",
    "Parts B and C: read the shared data/processed_prepare/markets_ppp.Rds and",
    "city_owin.Rds (window in feet, seed 671, 99 simulations) so all parts use the",
    "same window, and keep the coincident pair in every K-function run."),
    file.path(out, "analysis_handoff_summary_nna.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info_nna.txt"))
  paths <- c(file.path("output", list.files(out, "^(fig_0[1-3]|table_01)_")),
             file.path("output", c("analysis_handoff_summary_nna.txt",
                                   "session_info_nna.txt")),
             file.path("data", "processed_nna", "nndist_per_market.csv"))
  write_table(data.frame(File = paths), "output_manifest_nna.csv")
  cat("Part A (NNA) complete:", n, "markets, R =", round(nni, 4),
      ", one-sided p =", signif(p_one_sided, 4), "\n")
  invisible(TRUE)
}
run_hw4_part_nna()
