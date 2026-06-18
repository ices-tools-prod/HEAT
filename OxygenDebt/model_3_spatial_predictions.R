# ----------------------------
#
#   Model spatial profiles
#
#     * load profile parameters and fit spatial temporal models
#
# ----------------------------

# load packages etc.
header("model")

# Define paths
inputPath <<- file.path("OxygenDebt/Input", assessmentPeriod)
outputPath <<- file.path("OxygenDebt/Output", assessmentPeriod)

# start timer
t0 <- proc.time()

# get assessment period
#config <- jsonlite::fromJSON(file.path(inputPath, "config.json"))
if (assessmentPeriod == "2011-2016"){
  years <- c(2011, 2012, 2013, 2014, 2015, 2016)
} else if (assessmentPeriod == "2016-2021") {
  years <- c(2016, 2017, 2018, 2019, 2020, 2021)
}

# load gam fits ('gams')
check <- load(file.path(outputPath, "oxy_gam_fits.RData"))
if (check != "gams") {
  stop("Error loading gam fits!\n\tTry rerunning model_2_spatial_profiles.R")
}
rm(check)


# read helcom assessment areas
helcom <- sf::st_read(dsn = file.path(outputPath), layer = "oxy_areas", quiet = TRUE)

# read depth layer (spatial points) for prediction
bathy <- sf::st_read(dsn = file.path(outputPath), layer = "oxy_bathymetry", quiet = TRUE)

# drop regions not in models!
if ("Basin" %in% names(gams[[1]]$var.summary)) {
  bathy <- bathy[bathy$Basin %in% levels(gams[[1]]$var.summary$Basin),]
}


# quick plot
if (FALSE) {
  plot(bathy["depth"], pch = ".", col = rev(viridis::magma(50, alpha=0.5))[cut(bathy$depth, 50)])
  plot(sf::st_geometry(helcom), add = TRUE, border = "red")
}

# NOTES:
# * predict over x, y (noting region) for each year:
# *     depth changepoint 2 (lower halocline depth)
# *     O2_deficit at lower halocline
# *     O2_slope slope
 
# create prediction data
bathy <- cbind(bathy, sf::st_coordinates(bathy))
surfaces <- dplyr::rename(data.frame(bathy), x = X, y = Y)
surfaces$yday <- 1 # predict at january 1st - this will be dropped later anyway

if (FALSE) {
  # reduce size of prediction grid for testing
  surfaces <- surfaces[sample(1:nrow(surfaces), 1000),]
}

# create prediction data for each year
surfaces <-
  do.call(rbind,
    lapply(years, function(y) cbind(surfaces, Year = y))
  )

# do predictions
surfaces <-
  do.call(cbind.data.frame,
    c(list(surfaces),
      sapply(gams,
        function(g)
        {
          # set seasonal to zero for annual variation only
          g$coefficients[grep("yday",names(coef(g)))] <- 0
          # make predictions
          unname(c(exp(predict(g, newdata = surfaces))))
        },
        simplify = FALSE)
      )
  )


# From these surfaces, compute volume specific oxygen deficit.
# *    CHECK WITH SAS CODE:  numerically integrate over depth, using bathymetry.
# *         calculate O2 decifit by 1m^3, from depth changepoint 2 to bathymetric depth.
# *         compute volume of prediction "transect"
# *         finally compute total O2 deficit / total volume
# *         this results in a surface of volume specific oxygen debt.


# create depth vectors
depths <-
  lapply(1:nrow(surfaces),
         function(i) {
           seq(0.5, surfaces$depth[i], by = 1)
         })

# add volumes
surfaces$surface <-
  sapply(1:nrow(surfaces), function(i) {
    sum(depths[[i]] <= surfaces$halocline[i] - surfaces$depth_gradient[i])
})

surfaces$uhalo <-
  sapply(1:nrow(surfaces), function(i) {
    sum(depths[[i]] > surfaces$halocline[i] - surfaces$depth_gradient[i] &
        depths[[i]] <= surfaces$halocline[i])
})

surfaces$lhalo <-
  sapply(1:nrow(surfaces), function(i) {
    sum(depths[[i]] > surfaces$halocline[i] &
        depths[[i]] <= surfaces$halocline[i] + surfaces$depth_gradient[i])
})

surfaces$bottom <-
  sapply(1:nrow(surfaces), function(i) {
    sum(depths[[i]] > surfaces$halocline[i] + surfaces$depth_gradient[i])
})


# create salinity vector predictions
salinities <-
  lapply(1:nrow(surfaces),
         function(i) {
           surfaces$sali_surf[i] + surfaces$sali_dif[i] * pnorm(depths[[i]], surfaces$halocline[i], surfaces$depth_gradient[i])
         })

# add salinities
surfaces$surface_salinity <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] <= surfaces$halocline[i] - surfaces$depth_gradient[i]
    sum(salinities[[i]][filt])
})

surfaces$uhalo_salinity <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] - surfaces$depth_gradient[i] &
            depths[[i]] <= surfaces$halocline[i]
    sum(salinities[[i]][filt])
})

surfaces$lhalo_salinity <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] &
            depths[[i]] <= surfaces$halocline[i] + surfaces$depth_gradient[i]
    sum(salinities[[i]][filt])
})

surfaces$bottom_salinity <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] + surfaces$depth_gradient[i]
    sum(salinities[[i]][filt])
})


# add oxygen debt
surfaces$surface_O2debt <- 0 # by definition

surfaces$uhalo_O2debt <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] - surfaces$depth_gradient[i] &
            depths[[i]] <= surfaces$halocline[i]
    if (sum(filt) == 0) return(0)
    O2debt <- (depths[[i]][filt] - (surfaces$halocline[i] - surfaces$depth_gradient[i])) *
               surfaces$O2def_below_halocline[i] / (2 * surfaces$depth_gradient[i])
    sum(O2debt)
})

surfaces$lhalo_O2debt <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] &
            depths[[i]] <= surfaces$halocline[i] + surfaces$depth_gradient[i]
    if (sum(filt) == 0) return(0)
    O2debt <- (depths[[i]][filt] - (surfaces$halocline[i] - surfaces$depth_gradient[i])) *
               surfaces$O2def_below_halocline[i] / (2 * surfaces$depth_gradient[i])
    sum(O2debt)
})

surfaces$bottom_O2debt <-
  sapply(1:nrow(surfaces), function(i) {
    filt <- depths[[i]] > surfaces$halocline[i] + surfaces$depth_gradient[i]
    if (sum(filt) == 0) return(0)
    O2debt <- (depths[[i]][filt] - (surfaces$halocline[i] + surfaces$depth_gradient[i])) *
               surfaces$O2def_slope_below_halocline[i] + surfaces$O2def_below_halocline[i]
    sum(O2debt)
})


# remaining interesting quantities are:
#    Lhalo_area = NA,
#    uhalo_area = NA,
#    bottom_area = NA,
#    surface_area = NA,
#    O2debt = NA,
#    hypoxic_volume = NA,
#    cod_rep_volume = NA,
#    hypoxic_area = NA

# check
if (FALSE) {
  summary(surfaces)
}


if (FALSE) {
# a previous incorrect implementation
surfaces$depth_change_point1 <- surfaces$halocline - surfaces$sali_dif
surfaces$depth_change_point2 <- surfaces$halocline + surfaces$sali_dif
surfaces$oxygendebt <-
  sapply(1:nrow(surfaces),
    function(i) {
    # dont predict if halocline is below max depth
    if (surfaces$halocline[i] >= surfaces$depth[i]) return (NA_real_)
    # predict oxygen debt at 1m intervals and take the average
    depths <- seq(surfaces$halocline[i], surfaces$depth[i], by = 1)
    mean(oxy_profile(depths, surfaces[i,])) # sum(o2) / volume
  })
}

# checks
if (FALSE) {
  summary(surfaces)

  # plot surfaces to check
  tmp <- sf::st_as_sf(surfaces[surfaces$Year == 2011, ], coords = c("x", "y"), crs = 32634)

  plot(tmp["halocline"], pch = 16, cex = 0.5, col = rev(viridis::magma(50, alpha=0.5))[cut(tmp$halocline, 50)])
  plot(sf::st_geometry(helcom), add = TRUE, border = "red")
  title(main = "Halocline depth")

  plot(tmp["O2def_below_halocline"], pch = 16, cex = 0.5, col = rev(viridis::magma(50, alpha=0.5))[cut(tmp$O2def_below_halocline, 50)])
  plot(sf::st_geometry(helcom), add = TRUE, border = "red")
  title(main = "O2def below halocline")

  plot(tmp["oxygendebt"], pch = 16, cex = 0.5, col = rev(viridis::magma(50, alpha=0.5))[cut(tmp$oxygendebt, 50)])
  plot(sf::st_geometry(helcom), add = TRUE, border = "red")
  title(main = "Volume specific oxygen debt")

}

# write out --------------

save(surfaces, file = file.path(outputPath, "oxy_gam_predictions.RData"))


# done -------------------

message(sprintf("time elapsed: %.2f seconds", (proc.time() - t0)["elapsed"]))

