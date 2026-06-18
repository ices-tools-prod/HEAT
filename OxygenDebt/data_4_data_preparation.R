# ----------------------------
#
#   Prepare raw data for analysis
#
# ----------------------------

# load packages etc.
header("data")

# Define paths
inputPath <<- file.path("OxygenDebt/Input", assessmentPeriod)
outputPath <<- file.path("OxygenDebt/Output", assessmentPeriod)

# read in data ----------------------------------------------------------
library(data.table)
library(sf)

if (assessmentPeriod == '2011-2016') {
  ctd <- fread(input = file.path(inputPath, "StationSamples2011-2016CTD_2022-12-09.txt.gz"), sep = "\t", na.strings = "NULL", stringsAsFactors = FALSE, header = TRUE, check.names = TRUE)
  bot <- fread(input = file.path(inputPath, "StationSamples2011-2016BOT_2022-12-09.txt.gz"), sep = "\t", na.strings = "NULL", stringsAsFactors = FALSE, header = TRUE, check.names = TRUE)
  } else if (assessmentPeriod == "2016-2021") {
  ctd <- fread(input = file.path(inputPath, "StationSamples2016-2021CTD_2022-12-09.txt.gz"), sep = "\t", na.strings = "NULL", stringsAsFactors = FALSE, header = TRUE, check.names = TRUE)
  bot <- fread(input = file.path(inputPath, "StationSamples2016-2021BOT_2022-12-09.txt.gz"), sep = "\t", na.strings = "NULL", stringsAsFactors = FALSE, header = TRUE, check.names = TRUE)
  }

ctd <- ctd[QV.ODV.Depth..m. <= 1 & QV.ODV.Temperature..degC. <=1 & QV.ODV.Practical.Salinity..dmnless. <=1 & QV.ODV.Dissolved.Oxygen..ml.l. <= 1,.(Cruise, Year, Month, Day, Hour, Minute, Latitude = Latitude..degrees_north., Longitude = Longitude..degrees_east., Depth = Depth..m., Temperature = Temperature..degC., Salinity = Practical.Salinity..dmnless., Oxygen = Dissolved.Oxygen..ml.l.)]
bot <- bot[QV.ODV.Depth..m. <= 1 & QV.ODV.Temperature..degC. <=1 & QV.ODV.Practical.Salinity..dmnless. <=1 & QV.ODV.Dissolved.Oxygen..ml.l. <= 1,.(Cruise, Year, Month, Day, Hour, Minute, Latitude = Latitude..degrees_north., Longitude = Longitude..degrees_east., Depth = Depth..m., Temperature = Temperature..degC., Salinity = Practical.Salinity..dmnless., Oxygen = Dissolved.Oxygen..ml.l., Hydrogen_Sulphide = Hydrogen.Sulphide..H2S.S...umol.l.)]

oxy <- dplyr::full_join(ctd, bot, suffix = c(".ctd", ".bot"),
                        by = c("Year", "Month", "Day", "Hour", "Minute",
                               "Latitude", "Longitude",
                               "Cruise", "Depth"))

rm(ctd, bot)

# apply rules for which data to take:
#   keep ctd temperature over bottle, unless ctd is missing (i.e. is NA)
keep_x <- function(x, y) ifelse(is.na(x), y, x)
oxy$Oxygen <- keep_x(oxy$Oxygen.ctd, oxy$Oxygen.bot)
oxy$Temperature <- keep_x(oxy$Temperature.ctd, oxy$Temperature.bot)
oxy$Salinity <- keep_x(oxy$Salinity.ctd, oxy$Salinity.bot)
oxy$Type <- ifelse(is.na(oxy$Oxygen.ctd), "BOT", "CTD")
rm(keep_x)

# keep only data for the years given in the assessment
if (assessmentPeriod == '2011-2016') {
  oxy <- oxy[Year >= 2011 & Year <= 2016,]
} else if (assessmentPeriod == "2016-2021") {
  oxy <- oxy[Year >= 2016 & Year <= 2021,]
}

# Create profile ID
oxy[, ID := .GRP, by = .(Year, Month, Day, Hour, Minute, Latitude, Longitude, Cruise)]

# Sort by ID + Retain only columns we need
oxy <- oxy[order(ID), .(ID, Year, Month, Day, Latitude, Longitude, Depth, Type, Temperature, Salinity, Oxygen, Hydrogen_Sulphide)]

# ----------------------------
#
# Classify oxy stations into oxy areas
#
# ----------------------------

# Extract unique stations i.e. longitude/latitude pairs
stations <- unique(oxy[, .(Longitude, Latitude)])

# Make stations spatial keeping original latitude/longitude
stations <- st_as_sf(stations, coords = c("Longitude", "Latitude"), remove = FALSE, crs = 4326)

# Transform projection into UTM zone 34N
stations <- st_transform(stations, crs = 32634)

# Read oxy indicator modelling areas
oxy_areas <- st_read(file.path(outputPath, "oxy_areas.shp"))

# Classify stations into oxy areas
stations <- st_join(stations, oxy_areas, join = st_intersects)

# Delete stations not classified
stations <- na.omit(stations)

# Create x y columns with projected coordinates for later 
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

stations <- sfc_as_cols(stations)

# Remove spatial column and make into data table
stations <- st_set_geometry(stations, NULL) %>% as.data.table()

# Merge stations back into station samples - getting rid of station samples not classified into assessment units
oxy <- stations[oxy, on = .(Longitude, Latitude), nomatch = 0]

# ----------------------------
#
#  Merge auxiliary info table
#
# ----------------------------

# Read in auxilliary info
aux <- fread(file.path(inputPath, "Auxilliary.csv"))

# Merge
oxy <- merge(oxy, aux, by = "Basin", all.x = TRUE)
rm(aux)

# Drop regions that have no auxilliary info
oxy[!is.na(surfacedepth1),]

# ----------------------------
#
#  Write out data
#
# ----------------------------

fwrite(oxy, file.path(outputPath, "oxy.csv"))

