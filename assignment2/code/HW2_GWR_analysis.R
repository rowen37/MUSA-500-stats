# Homework 2, Part 3 | Fujia Zheng
# Run from the package root or code/: Rscript code/HW2_GWR_analysis.R
# Adaptive Gaussian geographically weighted regression with one local WLS
# regression per block group. Base R handles the fit; sf, spdep and ggplot2
# read/map the handoff and test residuals on the team's original Queen graph.

run_gwr <- function(root) {
  suppressPackageStartupMessages({
    library(sf)
    library(spdep)
    library(ggplot2)
  })
  out <- file.path(root, "output")
  data <- file.path(root, "data")
  local_dir <- file.path(data, "processed_gwr")
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  save_table <- function(value, filename) write.csv(value, file.path(out, filename), row.names = FALSE)
  shp <- st_read(file.path(data, "processed_part2", "RegressionFinal_Part2.shp"), quiet = TRUE)
  # The reference package matches every source feature to a projected
  # Census block-group polygon. Source shapefile coordinates are unused.
  fields <- c("PCTVACANT", "PCTSINGLES", "PCTBACHMOR", "LNNBELPOV")
  formula <- reformulate(fields, "LNMEDHVAL")
  observations <- st_drop_geometry(shp)
  X <- model.matrix(formula, observations)
  y <- observations$LNMEDHVAL
  n <- nrow(X)
  ids <- as.character(shp$POLY_ID)
  stopifnot(n == 1720L, !anyDuplicated(ids), all(is.finite(X)), all(is.finite(y)))

  # Metric locations are produced by the companion CRS verification step.
  # A coordinate file must match BOTH stable identifiers before use.
  ref <- read.csv(file.path(data, "gwr_reference_centroids.csv"), colClasses = c(AREAKEY = "character"))
  stopifnot(nrow(ref) == n, !anyDuplicated(ref$POLY_ID), !anyDuplicated(ref$AREAKEY))
  ref <- ref[match(ids, as.character(ref$POLY_ID)), , drop = FALSE]
  stopifnot(identical(ref$AREAKEY, as.character(shp$AREAKEY)))
  coords <- as.matrix(ref[c("x_m", "y_m")])
  stopifnot(all(is.finite(coords)), !anyDuplicated(as.data.frame(coords)))
  distances <- as.matrix(dist(coords))
  sorted <- apply(distances, 2L, sort)
  distances_sq <- distances^2
  rm(distances)
  queen <- read.gal(file.path(data, "RegressionData_queen.gal"), region.id = ids)
  wt <- nb2listw(queen, style = "W", zero.policy = TRUE)
  stopifnot(identical(attr(queen, "region.id"), ids),
            max(abs(lag.listw(wt, observations$OLS_RESID, zero.policy = TRUE) - observations$WT_RESID)) < 1e-6)

  # Gaussian w_ij = exp[-d_ij^2/(2 h_i^2)]. The radius h_i is the
  # distance to the k-th OTHER centroid. Points beyond h_i still contribute.
  # At each site, beta_i = (X'W_iX)^(-1)X'W_i y. Only hat diagonals
  # are needed for bandwidth selection (no dense n by n hat matrix).
  fit_at_k <- function(k, diagnostics = FALSE) {
    radius <- sorted[k + 1L, ]
    p <- ncol(X)
    b <- matrix(NA_real_, n, p, dimnames = list(NULL, c("Intercept", fields)))
    h <- numeric(n)
    if (diagnostics) {
      raw_variance <- matrix(NA_real_, n, p)
      condition <- numeric(n)
    }
    for (i in seq_len(n)) {
      w <- exp(-distances_sq[, i] / (2 * radius[i]^2))
      wx <- w * X
      inv <- chol2inv(chol(crossprod(X, wx)))
      b[i, ] <- inv %*% crossprod(wx, y)
      h[i] <- sum((X[i, ] %*% inv) * X[i, ])
      if (diagnostics) {
        # Diagonal of C_i C_i', where C_i = (X'W_iX)^(-1)X'W_i.
        raw_variance[i, ] <- colSums((wx %*% inv)^2)
        centered <- sweep(X[, -1, drop = FALSE], 2, colSums(wx[, -1, drop = FALSE]) / sum(w))
        weighted <- sqrt(w) * centered
        scaled <- sweep(weighted, 2, sqrt(colSums(weighted^2)), "/")
        condition[i] <- kappa(scaled, exact = TRUE)
      }
    }
    pred <- rowSums(b * X)
    residual <- y - pred
    rss <- sum(residual^2)
    trace <- sum(h)
    ll <- -0.5 * n * (1 + log(2 * pi) + log(rss / n))
    aicc <- -2 * ll + 2 * n * (trace + 1) / (n - trace - 2)
    if (!diagnostics) return(aicc)
    # SE convention: RSS/(n - trace(S)); same convention as spgwr's _se.
    se <- sqrt(raw_variance * rss / (n - trace))
    local_r2 <- vapply(seq_len(n), function(i) {
      w <- exp(-distances_sq[, i] / (2 * radius[i]^2))
      local_mean <- weighted.mean(y, w)
      1 - sum(w * residual^2) / sum(w * (y - local_mean)^2)
    }, numeric(1))
    fit <- data.frame(N = n, Bandwidth_k = k, Bandwidth_k_over_n = k/n,
                      Radius_min_m = min(radius), Radius_median_m = median(radius), Radius_max_m = max(radius),
                      Trace_S = trace, Effective_DF = n - trace, RSS = rss,
                      R2 = 1 - rss / sum((y - mean(y))^2), Log_likelihood = ll,
                      AIC = -2 * ll + 2 * (trace + 1), AICc = aicc,
                      AIC_spgwr_eq4_22 = -2 * ll + trace)
    list(b = b, se = se, ratios = b / se, local_r2 = local_r2, pred = pred,
         residual = residual, leverage = h, radius = radius, condition = condition, fit = fit)
  }
  # Evaluate a log-spaced coarse grid from 7 to all other centroids, then
  # refine the winning interval to integer resolution. Save all candidates.
  grid <- unique(round(exp(seq(log(ncol(X) + 2), log(n - 1), length.out = 22))))
  scores <- setNames(vapply(grid, fit_at_k, numeric(1)), grid)
  best <- which.min(scores)
  lo <- grid[max(1L, best - 1L)]
  hi <- grid[min(length(grid), best + 1L)]
  if (best == 1L) lo <- ncol(X)
  # A wide winning bracket is subdivided before the final sweep so the
  # candidate search remains practical even for a nearly global process.
  while (hi - lo > 60) {
    refine <- unique(round(seq(lo, hi, length.out = 11)))
    todo <- setdiff(refine, as.integer(names(scores)))
    scores <- c(scores, setNames(vapply(todo, fit_at_k, numeric(1)), todo))
    best <- which.min(scores[as.character(refine)])
    lo <- refine[max(1L, best - 1L)]
    hi <- refine[min(length(refine), best + 1L)]
  }
  todo <- setdiff(seq.int(lo, hi), as.integer(names(scores)))
  scores <- c(scores, setNames(vapply(todo, fit_at_k, numeric(1)), todo))
  k <- as.integer(names(scores)[which.min(scores)])
  search <- data.frame(Bandwidth_k = as.integer(names(scores)), AICc = unname(scores))
  search$Selected <- search$Bandwidth_k == k
  save_table(search[order(search$Bandwidth_k), ], "table_17_gwr_bandwidth.csv")
  g <- fit_at_k(k, TRUE)
  stopifnot(isTRUE(all.equal(g$fit$AICc, min(scores))), all(is.finite(g$se)),
            all(g$se > 0), all(is.finite(g$local_r2)), all(g$leverage >= 0 & g$leverage <= 1.00001))
  # Independently verify the normal-equation solution with R's QR-based
  # weighted least squares at three pre-specified locations.
  for (i in c(1L, n %/% 2L, n)) {
    w <- exp(-distances_sq[, i] / (2 * g$radius[i]^2))
    stopifnot(max(abs(lm.wfit(X, y, w)$coefficients - g$b[i, ])) < 1e-7)
  }
  save_table(g$fit, "table_18_gwr_fit.csv")
  local_summary <- function(z, metric, term) {
    q <- unname(quantile(z, c(0, .25, .5, .75, 1)))
    data.frame(Metric = metric, Term = term, Min = q[1], Q1 = q[2], Median = q[3],
               Q3 = q[4], Max = q[5], Positive_n = sum(z > 0), Negative_n = sum(z < 0))
  }
  summary_table <- do.call(rbind, lapply(seq_len(ncol(X)), function(j) {
    rbind(local_summary(g$b[, j], "Coefficient", colnames(g$b)[j]),
          local_summary(g$ratios[, j], "Coefficient / SE", colnames(g$b)[j]))
  }))
  summary_table <- rbind(summary_table, local_summary(g$local_r2, "Local R2", "All predictors"),
                         local_summary(g$condition, "Condition number", "All predictors"))
  save_table(summary_table, "table_19_gwr_local_summary.csv")
  local <- data.frame(POLY_ID = ids, AREAKEY = ref$AREAKEY, x_m = ref$x_m, y_m = ref$y_m,
                      Bandwidth_m = g$radius, Observed = y, Fitted = g$pred,
                      Residual = g$residual, Local_R2 = g$local_r2,
                      Condition_number = g$condition, Leverage = g$leverage)
  for (j in seq_len(ncol(X))) {
    term <- colnames(g$b)[j]
    local[[paste0(term, "_beta")]] <- g$b[, j]
    local[[paste0(term, "_se")]] <- g$se[, j]
    local[[paste0(term, "_ratio")]] <- g$ratios[, j]
  }
  write.csv(local, file.path(local_dir, "GWR_local_results.csv"), row.names = FALSE)

  # Hold geometry and the GAL weights fixed and permute the residuals.
  # spdep returns the observed statistic as the LAST entry in mc$res.
  set.seed(671)
  mc <- moran.mc(g$residual, wt, 999, alternative = "greater", zero.policy = TRUE)
  observed <- unname(mc$statistic)
  nulls <- head(mc$res, -1L)
  stopifnot(length(nulls) == 999L, abs(tail(mc$res, 1L) - observed) < 1e-10)
  upper <- (1 + sum(nulls >= observed)) / 1000
  lower <- (1 + sum(nulls <= observed)) / 1000
  moran <- data.frame(I = observed, P_two_sided = min(1, 2 * min(upper, lower)),
                      P_greater = upper, P_less = lower, Permutations = length(nulls), Seed = 671)
  save_table(moran, "table_20_gwr_moran.csv")

  map_theme <- theme_void(base_size = 11) + theme(
    plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 10),
    strip.text = element_text(face = "bold"), legend.text = element_text(size = 9),
    legend.title = element_text(size = 10), plot.margin = margin(12, 12, 12, 12))
  polygons <- st_read(file.path(data, "crs_reference", "Philadelphia_analysis_bg.gpkg"),
                      layer = "block_groups", quiet = TRUE)
  stopifnot(nrow(polygons) == n, !anyDuplicated(polygons$POLY_ID),
            isTRUE(st_crs(polygons)$epsg == 26918))
  polygons <- polygons[match(ids, as.character(polygons$POLY_ID)), c("POLY_ID", "AREAKEY")]
  stopifnot(identical(as.character(polygons$AREAKEY), ref$AREAKEY))
  outlines <- function() geom_sf(colour = NA)
  save_figure <- function(p, filename, width, height) {
    ggsave(file.path(out, filename), p, width = width, height = height, units = "in", dpi = 160,
           bg = "white", device = "png")
  }
  # Shared breaks keep magnitude and direction comparable across four panels.
  panel_labels <- c("Vacant housing (%)", "Detached single-family housing (%)",
                    "Bachelor's degree or higher (%)", "Log poverty count")
  panel <- do.call(rbind, lapply(seq_along(fields), function(j) {
    z <- polygons
    z$Variable <- panel_labels[j]
    z$Ratio <- g$ratios[, j + 1L]
    z
  }))
  panel$Variable <- factor(panel$Variable, levels = panel_labels)
  bins <- c("Below -2", "-2 to below 0", "0 to 2", "Above 2")
  panel$Bin <- factor(findInterval(panel$Ratio, c(-2, 0, 2)) + 1L,
                      levels = seq_along(bins), labels = bins)
  panel$Bin[panel$Ratio == 2] <- "0 to 2"
  map_ratio <- ggplot(panel, aes(fill = Bin)) + outlines() + facet_wrap(~Variable, ncol = 2) +
    scale_fill_manual(values = c("#9d1832", "#f2b6be", "#afcfe6", "#20578e"),
                      drop = FALSE, name = "Coefficient / SE") +
    labs(title = "Local coefficients relative to their standard errors") + map_theme
  save_figure(map_ratio, "fig_16_gwr_coefficient_se.png", 9, 8.5)
  polygons$Local_R2 <- g$local_r2
  map_fit <- ggplot(polygons, aes(fill = Local_R2)) + outlines() +
    scale_fill_viridis_c(name = expression(R^2), option = "C") +
    labs(title = "Local model fit", subtitle = "Higher values indicate a larger share of locally weighted variation explained") + map_theme
  save_figure(map_fit, "fig_17_gwr_local_r2.png", 7.1, 6.2)
  polygons$Residual <- g$residual
  map_res <- ggplot(polygons, aes(fill = Residual)) + outlines() +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", name = "Log residual") +
    labs(title = "GWR residuals", subtitle = "Positive values indicate observed log home value above the local prediction") + map_theme
  save_figure(map_res, "fig_18_gwr_residual_map.png", 7.1, 6.2)
  z <- as.numeric(scale(g$residual))
  lag_z <- lag.listw(wt, z, zero.policy = TRUE)
  scatter <- ggplot(data.frame(z, lag_z), aes(z, lag_z)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") + geom_point(size = .55, alpha = .5, colour = "#226388") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#ac3c39", linewidth = .8) +
    labs(title = "Residual Moran scatterplot", subtitle = sprintf("Queen weights; I = %.3f", observed),
         x = "Standardized GWR residual", y = "Spatially lagged standardized residual") + theme_minimal(base_size = 10)
  permutation_plot <- ggplot(data.frame(I = nulls), aes(I)) +
    geom_histogram(bins = 30, fill = "#879eaa", colour = "white") + geom_vline(xintercept = observed, colour = "#ac3c39") +
    labs(title = "Permutation reference", subtitle = sprintf("999 randomizations; two-sided p = %.3f", moran$P_two_sided),
         x = "Moran's I (red = observed)", y = "Count") + theme_minimal(base_size = 10)
  save_figure(patchwork::wrap_plots(scatter, permutation_plot, ncol = 2),
              "fig_19_gwr_residual_moran.png", 10, 4.4)
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info_part3.txt"))
  print(g$fit)
  print(moran)
  invisible(g$fit)
}

root <- if (file.exists("data/RegressionData_queen.gal")) "." else ".."
run_gwr(normalizePath(root))
