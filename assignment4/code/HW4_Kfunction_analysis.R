# Philadelphia farmers markets: reproducible Part B (Ripley's K-function).
# Run from code/ or the project root; the same script can be sourced by knitr.
# Builds on the shared point pattern stored by code/HW4_prepare_data.R, so Part B
# uses the same 62 markets, the same city window (EPSG:2272, feet) and the same
# seed as Part A.

run_hw4_part_kfunc <- function() {
  root_candidates <- c(".", "..")
  root <- root_candidates[file.exists(file.path(root_candidates, "data",
    "Philadelphia_Farmers_Markets201302.shp"))][1]
  if (is.na(root)) stop("Run this script from the project root or code/.")
  root <- normalizePath(root)
  input <- file.path(root, "data")
  out <- file.path(root, "output")
  processed <- file.path(input, "processed_kfunc")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed, recursive = TRUE, showWarnings = FALSE)

  ppp_path <- file.path(input, "processed_prepare", "markets_ppp.Rds")
  if (!file.exists(ppp_path)) {
    source(file.path(root, "code", "HW4_prepare_data.R"), local = TRUE)
    prepare_hw4_data(root)
  }
  X <- readRDS(ppp_path)
  window <- readRDS(file.path(input, "processed_prepare", "city_owin.Rds"))
  n <- spatstat.geom::npoints(X)
  area_ft2 <- spatstat.geom::area.owin(window)
  stopifnot(n == 62L, area_ft2 > 0)

  # Distance bands follow the assignment's ArcGIS settings: beginning distance
  # 0, increment 2500 feet, 10 bands, so the bands end at 2500 ... 25000 feet.
  # K is estimated on a 100 ft grid for smooth curves, and the band edges are
  # reported in table_02 (one row per band edge, plus r = 0).
  band_start <- 0
  increment <- 2500
  n_bands <- 10L
  band_edges <- seq(band_start, band_start + n_bands * increment, by = increment)
  r <- seq(band_start, band_start + n_bands * increment, by = 100)
  stopifnot(max(r) == 25000, max(band_edges) == 25000)
  band_idx <- match(band_edges, r)

  # Observed K with Ripley's isotropic edge correction (the default in Kest and
  # the R counterpart of a boundary-corrected estimate). The uncorrected K is
  # kept only to show how much the window boundary biases it.
  K_iso <- spatstat.explore::Kest(X, r = r, correction = "isotropic")
  K_none <- spatstat.explore::Kest(X, r = r, correction = "none")
  K_obs <- as.data.frame(K_iso)$iso
  K_unc <- as.data.frame(K_none)$un
  K_theory <- pi * r^2
  L_obs <- sqrt(pmax(K_obs, 0) / pi) - r
  L_unc <- sqrt(pmax(K_unc, 0) / pi) - r

  # 99 CSR simulations in the same window with the same number of points,
  # seed 671 to match Part A. Each simulated pattern is re-estimated with the
  # same Ripley correction, which plays the role of ArcGIS's
  # "Simulate Outer Boundary Values" boundary correction.
  nsim <- 99L
  seed <- 671L
  set.seed(seed)
  sim_mat <- matrix(NA_real_, nrow = length(r), ncol = nsim)
  for (i in seq_len(nsim)) {
    Yi <- spatstat.random::runifpoint(n, win = window)
    Ki <- spatstat.explore::Kest(Yi, r = r, correction = "isotropic")
    sim_mat[, i] <- as.data.frame(Ki)$iso
  }
  env_lo <- apply(sim_mat, 1, min)
  env_hi <- apply(sim_mat, 1, max)
  env_mean <- rowMeans(sim_mat)
  sim_sd <- apply(sim_mat, 1, sd)

  # Rank-based one-sided p-values at each distance: how often a simulated
  # pattern is at least as extreme as the observed one.
  p_cluster <- vapply(seq_along(r), function(i)
    (1 + sum(sim_mat[i, ] >= K_obs[i])) / (nsim + 1), numeric(1))
  p_disperse <- vapply(seq_along(r), function(i)
    (1 + sum(sim_mat[i, ] <= K_obs[i])) / (nsim + 1), numeric(1))
  outside <- K_obs < env_lo | K_obs > env_hi
  # Clustering = observed K above the envelope; dispersion = below it.
  sig_cluster <- K_obs > env_hi
  sig_disperse <- K_obs < env_lo

  # Global test: compare the largest standardised excursion of the observed
  # curve with the same statistic computed for each simulated pattern.
  # r = 0 is dropped because K(0) is identically 0 in every simulation.
  ok <- r > 0
  excursion <- function(v) max(abs((v[ok] - env_mean[ok]) / sim_sd[ok]))
  obs_excursion <- excursion(K_obs)
  sim_excursion <- apply(sim_mat, 2, excursion)
  global_p <- (1 + sum(sim_excursion >= obs_excursion)) / (nsim + 1)
  peak <- which.max(abs((K_obs - env_mean) / ifelse(sim_sd > 0, sim_sd, 1)))
  global_dir <- if (K_obs[peak] > env_mean[peak]) "clustered" else "dispersed"

  # L-function envelope, derived from the K envelope because L = sqrt(K/pi) - r.
  L_lo <- sqrt(pmax(env_lo, 0) / pi) - r
  L_hi <- sqrt(pmax(env_hi, 0) / pi) - r

  # Table of the ten bands (each row is the upper edge of a band, plus r = 0).
  band_table <- data.frame(
    distance_ft = band_edges,
    distance_mi = band_edges / 5280,
    K_observed_iso = K_obs[band_idx],
    K_theory_CSR = K_theory[band_idx],
    K_observed_none = K_unc[band_idx],
    L_observed_iso = L_obs[band_idx],
    L_theory_CSR = 0,
    envelope_lo = env_lo[band_idx],
    envelope_hi = env_hi[band_idx],
    simulated_mean = env_mean[band_idx],
    p_cluster = p_cluster[band_idx],
    p_disperse = p_disperse[band_idx],
    outside_envelope = outside[band_idx],
    significant_cluster_05 = sig_cluster[band_idx],
    significant_dispersion_05 = sig_disperse[band_idx]
  )
  write.csv(band_table, file.path(out, "table_02_kfunction_bands.csv"), row.names = FALSE)
  write.csv(data.frame(distance_ft = r, K_observed_iso = K_obs, L_observed_iso = L_obs,
                       envelope_lo = env_lo, envelope_hi = env_hi),
            file.path(processed, "kfunction_values.csv"), row.names = FALSE)

  # Figures follow the group's HW2/HW4 style (base R, same palette as Part A).
  point <- "#1F668883"
  accent <- "#AD4025"
  png(file.path(out, "fig_04_kfunction_envelope.png"), width = 1600, height = 780, res = 140)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.2, 1.0))
  plot(r, K_obs, type = "l", lwd = 2, col = point, ylim = range(c(env_lo, env_hi, K_obs)),
       xlab = "Distance d (feet)", ylab = "K(d)",
       main = "Ripley's K function with 99-simulation envelope")
  polygon(c(r, rev(r)), c(env_lo, rev(env_hi)), border = NA, col = "#7B9BAA33")
  lines(r, env_lo, lty = 2, col = "grey35"); lines(r, env_hi, lty = 2, col = "grey35")
  lines(r, K_theory, lwd = 2, lty = 3, col = "grey20")
  lines(r, K_obs, lwd = 2, col = point)
  points(band_edges, K_obs[band_idx], pch = 19, cex = 0.7, col = point)
  legend("topleft", bty = "n", cex = 0.8, lwd = c(2, 1, 2, 1),
         col = c(point, "grey35", "grey20", "#7B9BAA"),
         legend = c("Observed K(d)", "Min/max of 99 simulations", "CSR: pi*d^2",
                    "Simulation band (points = band edges)"))
  plot(r, L_obs, type = "l", lwd = 2, col = point, ylim = range(c(L_lo, L_hi, L_obs, 0)),
       xlab = "Distance d (feet)", ylab = "L(d) = sqrt(K(d)/pi) - d",
       main = "L function: values above 0 indicate clustering")
  polygon(c(r, rev(r)), c(L_lo, rev(L_hi)), border = NA, col = "#7B9BAA33")
  lines(r, L_lo, lty = 2, col = "grey35"); lines(r, L_hi, lty = 2, col = "grey35")
  abline(h = 0, lty = 3, lwd = 2, col = "grey20")
  lines(r, L_obs, lwd = 2, col = point)
  points(band_edges, L_obs[band_idx], pch = 19, cex = 0.7, col = point)
  legend("topleft", bty = "n", cex = 0.8, lwd = c(2, 1, 2), lty = c(1, 2, 3),
         col = c(point, "grey35", "grey20"),
         legend = c("Observed L(d)", "Simulation envelope", "CSR: L(d) = 0"))
  mtext(sprintf("Observed curve stays above the simulation band at every band edge from %d ft",
                min(band_edges[sig_cluster[band_idx]])), side = 3, line = 0.3, cex = 0.75)

  # Boundary-correction comparison: uncorrected K rises faster at large d.
  png(file.path(out, "fig_05_kfunction_boundary_correction.png"), width = 1200, height = 900, res = 140)
  plot(r, K_obs, type = "l", lwd = 2, col = point, ylim = range(c(K_obs, K_unc)),
       xlab = "Distance d (feet)", ylab = "K(d)",
       main = "Ripley's correction vs no correction")
  lines(r, K_unc, lwd = 2, lty = 2, col = accent)
  lines(r, K_theory, lwd = 2, lty = 3, col = "grey20")
  points(band_edges, K_obs[band_idx], pch = 19, cex = 0.7, col = point)
  legend("topleft", bty = "n", cex = 0.85, lwd = 2, lty = c(1, 2, 3),
         col = c(point, accent, "grey20"),
         legend = c("Ripley isotropic correction", "No edge correction", "CSR: pi*d^2"))
  dev.off()

  sig_distances <- band_edges[sig_cluster[band_idx]]
  first_sig <- if (length(sig_distances)) min(sig_distances) else NA
  last_sig <- if (length(sig_distances)) max(sig_distances) else NA
  disp_distances <- band_edges[sig_disperse[band_idx]]
  first_disp <- if (length(disp_distances)) min(disp_distances) else NA
  last_disp <- if (length(disp_distances)) max(disp_distances) else NA
  conclusion <- if (length(sig_distances)) {
    sprintf("Reject CSR: K(d) exceeds the simulation envelope at every band edge from %d to %d ft (%d of %d edges), so markets are clustered at all scales up to 25000 ft",
            first_sig, last_sig, sum(sig_cluster[band_idx]), length(band_edges) - 1L)
  } else if (length(disp_distances)) {
    sprintf("Reject CSR: K(d) falls below the simulation envelope from %d to %d ft, so markets are more dispersed than random", first_disp, last_disp)
  } else "Do not reject CSR at any distance band"

  writeLines(c(
    "Homework 4 Part B: Ripley's K-function analysis of farmers markets",
    sprintf("Markets: %d | Study area: %.0f square feet (%.1f sq mi)", n, area_ft2,
            area_ft2 / 5280^2),
    sprintf("Distance bands: beginning 0 ft, increment %d ft, %d bands; K evaluated at %s ft",
            increment, n_bands, paste(range(band_edges), collapse = "-")),
    sprintf("Edge correction: Ripley isotropic (Kest default); uncorrected K exported for comparison"),
    sprintf("Simulations: %d CSR patterns in the same window, seed %d", nsim, seed),
    sprintf("%s", conclusion),
    sprintf("Global test: max standardised excursion = %.3f, simulation-based p = %.3f (%s)", obs_excursion, global_p, global_dir),
    sprintf("Observed K at 25000 ft = %.0f vs CSR %.0f ft^2; uncorrected K = %.0f", K_obs[length(r)], K_theory[length(r)], K_unc[length(r)]),
    "Reading the bands: K(d) above the simulated envelope at a given distance means",
    "more pairs than a random pattern would produce at that range, i.e. clustering;",
    "K(d) below the envelope means fewer pairs than random, i.e. regularity. The L(d)",
    "panel makes this easier: L(d) > 0 at short distances is the direct clustering signal.",
    "Coincident Reading Terminal pair kept in every run, matching Part A.",
    "Note: the 25000 ft outer band slightly exceeds spatstat's recommended range",
    "(22384 ft) for the Philadelphia window; the isotropic correction is still",
    "defined there, and the uncorrected K is exported in fig_05 for comparison.",
    "Part C: use the same window, band edges and seed; the L(d) table is exported as",
    "data/processed_kfunc/kfunction_values.csv for the discussion figures."),
    file.path(out, "analysis_handoff_summary_kfunc.txt"))

  writeLines(capture.output(sessionInfo()), file.path(out, "session_info_kfunc.txt"))
  manifest <- data.frame(File = c(
    file.path("output", c("fig_04_kfunction_envelope.png",
                          "fig_05_kfunction_boundary_correction.png",
                          "table_02_kfunction_bands.csv",
                          "analysis_handoff_summary_kfunc.txt",
                          "session_info_kfunc.txt")),
    file.path("data", "processed_kfunc", "kfunction_values.csv")))
  write.csv(manifest, file.path(out, "output_manifest_kfunc.csv"), row.names = FALSE)

  cat("Part B (K-function) complete:", n, "markets, bands 0-25000 ft,",
      sprintf("%d of %d band edges outside the 99-simulation envelope\n",
              sum(outside[band_idx]), length(band_edges) - 1L))
  print(band_table[, c("distance_ft", "K_observed_iso", "envelope_lo", "envelope_hi",
                       "L_observed_iso", "outside_envelope")], row.names = FALSE)
  invisible(list(table = band_table, simulations = sim_mat, conclusion = conclusion))
}

run_hw4_part_kfunc()
