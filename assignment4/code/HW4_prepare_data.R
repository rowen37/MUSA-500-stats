# Build the shared spatstat objects for every HW4 analysis (Parts A, B and C).
# Reads the three course shapefiles, validates them, and stores the farmers
# market point pattern, the study-area window and the zip polygons under
# data/processed_prepare/. Run from code/ or the project root.
prepare_hw4_data <- function(root = c("..", ".")) {
  root <- root[file.exists(file.path(root, "data",
    "Philadelphia_Farmers_Markets201302.shp"))][1]
  if (is.na(root)) stop("Run this script from the project root or code/.")
  root <- normalizePath(root)
  input <- file.path(root, "data")
  processed <- file.path(input, "processed_prepare")
  dir.create(processed, recursive = TRUE, showWarnings = FALSE)

  markets <- sf::st_read(file.path(input, "Philadelphia_Farmers_Markets201302.shp"),
                         quiet = TRUE, stringsAsFactors = FALSE)
  city <- sf::st_read(file.path(input, "Philadelphia.shp"), quiet = TRUE)
  zip <- sf::st_read(file.path(input, "Philadelphia_ZipCodes.shp"), quiet = TRUE)
  # The course data are delivered in NAD83 StatePlane Pennsylvania South (ft),
  # EPSG:2272. The foot is the native distance unit for the whole assignment.
  stopifnot(nrow(markets) == 62L, nrow(city) == 1L, nrow(zip) == 47L,
            all(sf::st_geometry_type(markets) == "POINT"),
            all(sf::st_geometry_type(city) %in% c("POLYGON", "MULTIPOLYGON")),
            all(sf::st_geometry_type(zip) %in% c("POLYGON", "MULTIPOLYGON")),
            all(vapply(list(markets, city, zip),
                       function(x) sf::st_crs(x)$epsg == 2272L, logical(1))))
  # Every market must sit inside the study-area polygon that defines the window.
  stopifnot(all(lengths(sf::st_within(markets, city)) == 1L))

  # as.owin() understands sf objects in current spatstat.geom; fall back to the
  # raw ring coordinates if the method is unavailable in the installed version.
  sf_to_owin <- function(poly_sf) {
    win <- tryCatch(spatstat.geom::as.owin(poly_sf), error = function(e) NULL)
    if (!is.null(win)) return(win)
    rings <- list()
    for (g in sf::st_geometry(poly_sf)) {
      pieces <- if (inherits(g, c("multipolygon", "MULTIPOLYGON"))) g else list(g)
      for (part in pieces) for (ring in part) {
        rings[[length(rings) + 1L]] <- ring[, 1:2, drop = FALSE]
      }
    }
    anticlockwise <- function(m) {
      if (identical(m[1, ], m[nrow(m), ])) m <- m[-nrow(m), , drop = FALSE]
      nxt <- c(seq_len(nrow(m) - 1L) + 1L, 1L)
      if (sum(m[, 1L] * m[nxt, 2L] - m[nxt, 1L] * m[, 2L]) < 0) {
        m <- m[nrow(m):1, , drop = FALSE]
      }
      m
    }
    spatstat.geom::owin(poly = lapply(rings, anticlockwise))
  }
  window <- sf_to_owin(city)
  area_sf <- as.numeric(sf::st_area(city))
  # The window must reproduce the shapefile area to within rounding tolerance.
  stopifnot(abs(spatstat.geom::area.owin(window) - area_sf) / area_sf < 1e-6)
  spatstat.geom::unitname(window) <- "feet"

  xy <- sf::st_coordinates(markets)
  stopifnot(all(is.finite(xy)))
  # Two separately listed markets share one address: Reading Terminal Market
  # and Fair Food Farmstand @ RTM (12th St and Arch St). Keep both points and
  # verify no other coordinates coincide.
  coincident <- duplicated(xy[, 1:2, drop = FALSE]) |
    duplicated(xy[, 1:2, drop = FALSE], fromLast = TRUE)
  stopifnot(sum(coincident) == 2L,
            all(markets$NAME[coincident] == c("Reading Terminal Market",
                                               "Fair Food Farmstand @ RTM")))
  X <- spatstat.geom::ppp(xy[, 1L], xy[, 2L], window = window)
  spatstat.geom::unitname(X) <- "feet"
  stopifnot(spatstat.geom::npoints(X) == 62L,
            all(spatstat.geom::inside.owin(X$x, X$y, window)))

  saveRDS(X, file.path(processed, "markets_ppp.Rds"))
  saveRDS(window, file.path(processed, "city_owin.Rds"))
  sf::st_write(zip, file.path(processed, "zip_sf.gpkg"),
               layer = "zip_codes", delete_layer = TRUE, quiet = TRUE)
  utils::write.csv(data.frame(OBJECTID_1 = markets$OBJECTID_1,
                               NAME = markets$NAME, ZIP_CODE = markets$ZIP_CODE,
                               x_ft = xy[, 1L], y_ft = xy[, 2L]),
                   file.path(processed, "markets_points.csv"), row.names = FALSE)
  invisible(X)
}

if (sys.nframe() == 0L) prepare_hw4_data()
