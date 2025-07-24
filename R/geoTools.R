#' @title Transform a spatial feature object to csv.
#' @description Get palaeomaps in .csv format to create binary map files in BayesTraits.
#' @author Jorge Avaria-Llautureo (j.l.avaria@reading.ac.uk).
#' @param age Target age in million years. Default 150. Based on the \code{\link[rgplates]{reconstruct}} function.
#' @param model The map reconstruction model. Default is "PALEOMAP". Based on the \code{\link[rgplates]{reconstruct}} function - see full list of options there.
#' @param grid The spatial grid to create a mask of land/ocean.
#' @param extent A numeric vector defining the mask bounding box. Default is the global extent c(-180, 180, -90, 90). More information about this can be found in the help pages for the \code{\link[terra]{ext}} function from the \code{terra} package.
#' @param resolution A numeric value defining the grid cell size in degrees for the WGS84 crs. Default is 2 degrees. There are no explicit bounds for this but values of <0.01 can result in slow compute times.
#' @param outputfile A character string defining the filename (and full filepath) where the csv file should be saved. By default, the file will be saved simply as "masked.csv" in the present working directory.If left as blank (""), no file will be output.
#' @return A \code{data frame} with longitude, latitude, and mask value for land/ocean. This table is also written to the filepath specified by \code{outputfile}.
#' @importFrom terra rast crs rasterize values
#' @importFrom rgplates reconstruct
#' @export
#' @examples
#' df <- sf_to_csv(age = 200)
#' df <- sf_to_csv(age = 200, extent = c(-80, 30, -90, 90), resolution = 1)
#'

# function
sf_to_csv <- function(age = 150, grid = NULL, extent = c(-180, 180, -90, 90),
                      resolution = 2, model = "PALEOMAP", outputfile = "./masked.csv") {
  # Validate inputs
  if (!is.numeric(age) || age < 0) stop("Age must be a non-negative numeric value.")
  if (is.null(grid) && !is.numeric(extent)) stop("Extent must be a numeric vector of length 4 (xmin, xmax, ymin, ymax).")
  if (is.null(grid) && !is.numeric(resolution)) stop("Resolution must be numeric.")

  # Reconstruct coastlines
  map <- reconstruct("coastlines", age = age, model = model)

  # Remove "coastlines" class to ensure compatibility with terra
  class(map) <- c("sf", "data.frame")

  # Create or validate raster template
  if (is.null(grid)) {
    grid <- rast(xmin = extent[1], xmax = extent[2], ymin = extent[3], ymax = extent[4],
                        res = resolution, crs = (map))
  } else if (!inherits(grid, "SpatRaster")) {
    stop("Grid must be a SpatRaster object.")
  }

  # Rasterize to binary mask (1 for land, 0 for ocean)
  map_rast <- rasterize(map, grid)
  terra::values(map_rast) <- ifelse(terra::values(map_rast) > 0, 1,
                                    terra::values(map_rast))
  terra::values(map_rast) <- ifelse(is.na(values(map_rast)), 0,
                                    terra::values(map_rast))

  # Convert to data frame with lon, lat, mask
  df <- as.data.frame(map_rast, xy = TRUE, na.rm = FALSE)
  colnames(df) <- c("lon", "lat", "mask")

  if(outputfile != "") write.csv(df, file = outputfile,  row.names = F, quote = F)

  # Return data frame
  return(df)
}



