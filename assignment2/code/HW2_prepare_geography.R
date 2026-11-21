# Rebuild the independently georeferenced GWR locations from the bundled
# US Census Bureau 2000 block-group boundaries. Run from the code directory.
prepare_gwr_geography <- function(root = "..") {
  root <- normalizePath(root, mustWork = TRUE)
  ref <- file.path(root, "data", "crs_reference")
  src <- file.path(ref, "tl_2010_42101_bg00.zip")
  checksum <- "e7cd1443c3ae2a3d1e307474332392af"
  stopifnot(identical(unname(tools::md5sum(src)), checksum))

  input <- sf::st_read(file.path(root, "data/map_RegressionData/RegressionData.shp"),
                       quiet = TRUE)
  target <- sf::st_drop_geometry(input)[, c("POLY_ID", "AREAKEY")]
  official <- sf::st_read(paste0("/vsizip/", normalizePath(src)), quiet = TRUE)
  stopifnot(nrow(target) == 1720L, nrow(official) == 1816L,
            !anyNA(target$AREAKEY), !anyDuplicated(target$AREAKEY),
            !anyDuplicated(target$POLY_ID), !anyDuplicated(official$BKGPIDFP00),
            sf::st_crs(official)$epsg == 4269L)
  idx <- match(target$AREAKEY, official$BKGPIDFP00)
  stopifnot(!anyNA(idx), identical(target$AREAKEY, official$BKGPIDFP00[idx]))
  polygons <- sf::st_sf(target, geometry = sf::st_geometry(official)[idx])
  polygons <- sf::st_transform(polygons, 26918) # NAD83 / UTM zone 18N, metres
  centroids <- sf::st_centroid(sf::st_geometry(polygons))
  xy <- sf::st_coordinates(centroids)
  stopifnot(all(is.finite(xy)), !anyDuplicated(as.data.frame(xy)),
            sf::st_crs(polygons)$units_gdal == "metre")
  utils::write.csv(data.frame(target, x_m = xy[, 1L], y_m = xy[, 2L]),
                   file.path(root, "data/gwr_reference_centroids.csv"),
                   row.names = FALSE)
  sf::st_write(polygons, file.path(ref, "Philadelphia_analysis_bg.gpkg"),
               layer = "block_groups", delete_layer = TRUE, quiet = TRUE)
  invisible(polygons)
}

if (sys.nframe() == 0L) prepare_gwr_geography()
