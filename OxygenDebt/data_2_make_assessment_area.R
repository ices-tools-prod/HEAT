# ----------------------------
#
#   make assessment area
#
#     * merge helcom assessment areas and balsem areas
#
# ----------------------------

# load packages etc.
header("data")

# Define paths
inputPath <<- file.path("OxygenDebt/Input", assessmentPeriod)
outputPath <<- file.path("OxygenDebt/Output", assessmentPeriod)

# ----------------------------
#
#  Create new shapefile for helcom areas
#
# ----------------------------

# read helcom and drop non SEA areas
helcom <- sf::st_read(dsn = inputPath, layer = "AssessmentUnits", quiet = TRUE)
helcom <- helcom[grep("^SEA-", helcom$Code), ]
# transform to utm34
helcom <- sf::st_transform(helcom, 32634)

# merge areas (need to buffer a bit for a clean merge)
helcom_balsem <- sf::st_union(sf::st_buffer(helcom, dist = 10))

# read baltsem, and cut over helcom
baltsem <- sf::st_read(dsn = inputPath, layer = "Baltsem_utm34", quiet = TRUE)
helcom_balsem <- sf::st_intersection(baltsem, helcom_balsem)
if (!all(c("Bo_Basin", "Basin") %in% names(helcom_balsem))) {
  sample_points <- sf::st_point_on_surface(helcom_balsem)
  sample_points <- sf::st_join(sample_points, baltsem)
  helcom_balsem <- cbind(helcom_balsem, sample_points[, setdiff(names(sample_points), names(helcom_balsem)), drop = FALSE])
}

# fix names
helcom_balsem$Bo_Basin <- gsub("Ö", "oe", helcom_balsem$Bo_Basin)
helcom_balsem$Bo_Basin <- iconv(helcom_balsem$Bo_Basin, "UTF-8", "ASCII", sub = "")
helcom_balsem$Bo_Basin <- gsub("oeresund", "Oeresund", helcom_balsem$Bo_Basin)
helcom_balsem$Basin <- helcom_balsem$Bo_Basin

# keep only certain areas
helcom_balsem <-
  helcom_balsem[helcom_balsem$Basin %in% c("Arkona Basin",
                                           "Baltic Proper",
                                           "Bornholm Basin",
                                           "Bothnian Bay",
                                           "Bothnian Sea",
                                           "Gulf of Finland"), ]

# merge Gulf of Finland with Baltic Proper
tmp <- sf::st_union(helcom_balsem[helcom_balsem$Basin %in% c("Baltic Proper", "Gulf of Finland"), ])
merged <- sf::st_as_sf(data.frame(Basin = "Baltic Proper", stringsAsFactors = FALSE), geometry = sf::st_geometry(tmp))
helcom_balsem <- rbind(
  helcom_balsem[helcom_balsem$Basin != "Baltic Proper" & helcom_balsem$Basin != "Gulf of Finland", c("Basin", "geometry")],
  merged
)

# check
if (FALSE) {
  plot(helcom_balsem["Basin"], col = gplots::rich.colors(nrow(helcom_balsem), alpha=0.5))
  text(sf::st_coordinates(sf::st_centroid(sf::st_geometry(helcom_balsem))), as.character(helcom_balsem$Basin), cex = 0.7)
}

# write
sf::st_write(helcom_balsem["Basin"], dsn = outputPath, layer = "oxy_areas", driver = "ESRI Shapefile", delete_layer = TRUE)

# add to zip
zip(file.path(outputPath, "oxy_areas.zip"), file.path(outputPath, dir(outputPath, pattern = "^oxy_areas*")))
