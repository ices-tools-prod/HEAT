
# load header function into top of search list
while("tools:oxydebt_funs" %in% search()) detach("tools:oxydebt_funs")
sys.source("scripts/OxygenDebt/zz_header.R", envir = attach(NULL, name = "tools:oxydebt_funs"))

# data
source("scripts/OxygenDebt/data_1_download.R")
source("scripts/OxygenDebt/data_2_make_assessment_area.R")
source("scripts/OxygenDebt/data_3_make_bathymetry_layer.R")
source("scripts/OxygenDebt/data_4_data_preparation.R")

# modelling prep
source("scripts/OxygenDebt/input_1_data_cleaning.R")

# modelling
source("scripts/OxygenDebt/model_1_profiles.R")
source("scripts/OxygenDebt/model_2_spatial_profiles.R")
source("scripts/OxygenDebt/model_3_spatial_predictions.R")

# output
source("scripts/OxygenDebt/output_1_tables.R")
source("scripts/OxygenDebt/output_2_plots.R")


if (FALSE) {
# some ideas for modelling oxygen debt profiles spatially

# load packages etc.
header("model")

# read in data
oxy <- read.csv("analysis/input/OxygenDebt/oxy_clean.csv")
bathy <- sf::st_read(dsn = "data/OxygenDebt/shapefiles", layer = "helcom_bathymetry", quiet = TRUE)
helcom <- sf::st_read(dsn = "data/OxygenDebt/shapefiles", layer = "helcom_areas", quiet = TRUE)

# fit salinity model everywhere


#g <- gam(Salinity ~ te(x, y, Depth, yday, d = c(2,1,1), m = c(2,1,2), k = c(30,5,6), bs = c("tp","tp","cc")), data = data)
g <- gam(Oxygen_deficit ~ #s(Depth, ID, bs = "fs", k = 5) + # a smooth random effect by station
                          ti(x, y, Depth, d = c(2,1), m = c(2,1), k = c(50,5), bs = c("tp","tp"), fx = c(FALSE, TRUE)) +
                          ti(yday, Depth, d = c(1,1), m = c(2,1), k = c(6,5), bs = c("cc","tp"))
         , data = data)
summary(g)
AIC(g)

# choose a location to plot
plot(sf::st_geometry(makeSpatial(data)), pch = 1)
xy <- locator(n=1)

# get depth
plot(sf::st_geometry(helcom))
loc <- bathy[which.min(colSums((t(sf::st_coordinates(sf::st_geometry(bathy))) - unlist(xy))^2)),]
points(xy$x, xy$y, col = "red", pch=16)
points(loc, col = "blue")

# plot model at location with nearby data
ids <- data$ID[which.min((data$x - xy$x)^2 + (data$y - xy$y)^2)]
pdata <- data[data$ID %in% ids, ]

pred <- data.frame(Depth = seq(g$var.summary$Depth[1], loc$depth, length = 100))
pred[c("x", "y")] <- xy
pred$yday <- pdata$yday[1]
pred$fit <- predict(g, newdata = pred)

plot(pred$Depth, pred$fit, type = "l", col = "red")
points(pdata$Depth, pdata$Oxygen_deficit, col = "red")

}
