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
#helcom <- rgdal::readOGR(inputPath, "AssessmentUnits", verbose = FALSE)
helcom <- st_read(file.path(inputPath, "AssessmentUnits.shp"))
helcom <- helcom[grep("^SEA-", helcom$Code),]

# transform to utm34
#helcom <- sp::spTransform(helcom, sp::CRS("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
helcom <- st_transform(helcom, sp::CRS("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# merge areas (need to buffer a bit for a clean merge)
#helcom_balsem <- rgeos::gUnaryUnion(rgeos::gBuffer(helcom, byid = TRUE, width = 10))
helcom_balsem <- st_union(st_buffer(helcom, 10))

# read baltsem, and cut over helcom
#baltsem <- rgdal::readOGR(inputPath, "Baltsem_utm34", verbose = FALSE)
baltsem <- st_read(file.path(inputPath, "Baltsem_utm34.shp"))
#helcom_balsem <- rgeos::gIntersection(baltsem, helcom_balsem, byid = TRUE)
helcom_balsem <- st_intersects(baltsem, helcom_balsem, sparse = TRUE)
for (i in 1:length(helcom_balsem)) helcom_balsem@polygons[[i]]@ID <- paste(i)
data <-
  do.call(rbind,
          lapply(1:length(helcom_balsem),
                 function(i) sp::over(sp::spsample(helcom_balsem[i,], 1, type = "random"), baltsem)))
helcom_balsem <- sp::SpatialPolygonsDataFrame(helcom_balsem, data)

# fix names
helcom_balsem$Bo_Basin <- gsub("Ã", "oe", helcom_balsem$Bo_Basin)
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
                                           "Gulf of Finland"),]

# merge Gulf of Finland with Baltic Proper
tmp <- rgeos::gUnaryUnion(helcom_balsem[helcom_balsem$Basin %in% c("Baltic Proper", "Gulf of Finland"),])
helcom_balsem@polygons[[which(helcom_balsem$Basin == "Baltic Proper")]] <- tmp@polygons[[1]]
helcom_balsem <- helcom_balsem[helcom_balsem$Basin != "Gulf of Finland",]

# check
if (FALSE) {
  sp::plot(helcom_balsem, col = gplots::rich.colors(nrow(helcom_balsem), alpha=0.5))
  text(sp::coordinates(helcom_balsem), as.character(helcom_balsem$Basin), cex = 0.7)
}

# write
rgdal::writeOGR(helcom_balsem["Basin"], outputPath, "oxy_areas", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# add to zip
zip(file.path(outputPath, "oxy_areas.zip"), file.path(outputPath, dir(outputPath, pattern = "^oxy_areas*")))
