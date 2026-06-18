# ----------------------------
#
#   make depth spatial points layer
#
#     * use baltsem data to create a spatial points layer of depths
#
# ----------------------------

# load packages etc.
header("data")

# Define paths
inputPath <<- file.path("OxygenDebt/Input", assessmentPeriod)
outputPath <<- file.path("OxygenDebt/Output", assessmentPeriod)

# read raw depth points
bathy <- read.csv(file.path(inputPath, "BALTIC_BATHY_BALTSEM.csv"))
names(bathy) <- cleanColumnNames(names(bathy))
bathy <- dplyr::rename(bathy, depth = dybde)
bathy <- bathy[c("x", "y", "depth")]

# make into spatial points object (note implicit utm34 in BALTIC_BATHY_BALTSEM.csv)
bathy <- sf::st_as_sf(bathy, coords = c("x", "y"), crs = 32634 )

# trim to extent of assessment units
helcom <- sf::st_read(dsn = outputPath, layer = "oxy_areas", quiet = TRUE)
bathy <- sf::st_as_sf(bathy, coords = c("x", "y"), crs = 32634)
keep <- lengths(sf::st_intersects(bathy, sf::st_union(helcom))) > 0
bathy <- bathy[keep, ]

# join points with new helcom polygons
bathy$Basin <- sf::st_join(bathy, helcom["Basin"])$Basin

# check
if (FALSE) {
  plot(bathy["depth"], pch = ".", col = gplots::rich.colors(50, alpha=0.5)[cut(bathy$depth, 50)])
  plot(sf::st_geometry(helcom), border = "red", add = TRUE)
}

# write
sf::st_write(bathy[c("depth", "Basin")], dsn = outputPath, layer = "oxy_bathymetry", driver = "ESRI Shapefile", delete_layer = TRUE)

# add to zip
zip(file.path(outputPath, "oxy_bathymetry.zip"), file.path(outputPath, dir(outputPath, pattern = "^oxy_bathymetry*")))
