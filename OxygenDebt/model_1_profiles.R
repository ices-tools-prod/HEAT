# ----------------------------
#
#   Model profiles
#
#     * load oxy data and fit salinity and oxygen profiles
#
# ----------------------------


# load packages etc.
header("model")

# Define paths
inputPath <<- file.path("OxygenDebt/Input", assessmentPeriod)
outputPath <<- file.path("OxygenDebt/Output", assessmentPeriod)

# start timer
t0 <- proc.time()

# read in data
oxy <- fread(file.path(outputPath, "oxy_clean.csv"))

# run all fits
profiles <-
  do.call(rbind,
          lapply(unique(oxy$ID),
                 function(id) doonefit_full(subset(oxy, ID == id), ID = id, debug = TRUE)))

# join all profile level variables onto profiles
out <- unique(dplyr::select(oxy, -Depth, -Type,
                                 -Temperature, -Salinity, -Oxygen, -Hydrogen_Sulphide,
                                 -Oxygen_ml, -Oxygen_deficit, -censor))

rownames(out) <- paste(out$ID)

# join this onto the profiles
profiles <- merge(profiles, out)


# select which data are reliable -----
drop_sal <- profiles$sali_dif.se > 2.5 | # only use profiles that are based on an estimate of the salinity difference estimate with +- 5 accuracy
            profiles$depth_change_point2.se > 10 | # only use salinity profiles that are based on an estimate of the lower halocline estimate with +- 20m accuracy
            profiles$halocline.se > 10 | # only use salinity profiles with halocline depth estimated to +- 20m accuracy
            profiles$sali_dif < 0 | profiles$sali_dif > 17 |   # big salinity difference estimates with resonable precision - dubious.
            profiles$halocline == profiles$salmod_halocline_lower | profiles$halocline > 100 |
            profiles$depth_gradient > 45 |
            profiles$depth_change_point1 < profiles$surfacedepth2 | profiles$depth_change_point1 > 90

drop_sal[is.na(drop_sal)] <- FALSE
profiles$sali_dif[drop_sal] <- NA
profiles$halocline[drop_sal] <- NA
profiles$depth_gradient[drop_sal] <- NA
profiles$depth_change_point1[drop_sal] <- NA
profiles$depth_change_point2[drop_sal] <- NA

profiles$O2def_below_halocline[drop_sal] <- NA
profiles$O2def_slope_below_halocline[drop_sal] <- NA

# drop off badly estmated O2 def slopes
drop_O2def <- profiles$O2def_slope_below_halocline > 1.5

profiles$O2def_slope_below_halocline[drop_O2def] <- NA

# write out
fwrite(profiles, file.path(outputPath, "oxy_profiles.csv"))

# done -------------------

message(sprintf("time elapsed: %.2f seconds", (proc.time() - t0)["elapsed"]))

# checks -----------------

if (FALSE) {
  par(mfrow = c(2,2))
  with(profiles[!is.na(profiles$sali_dif),],
    {
      plot(sali_dif, sali_dif.se)
      plot(sali_dif, halocline.se)
      plot(sali_dif, depth_gradient.se)
      plot(sali_dif, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$halocline),],
    {
      plot(halocline, sali_dif.se)
      plot(halocline, halocline.se)
      plot(halocline, depth_gradient.se)
      plot(halocline, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$depth_gradient),],
    {
      plot(depth_gradient, sali_dif.se)
      plot(depth_gradient, halocline.se)
      plot(depth_gradient, depth_gradient.se)
      plot(depth_gradient, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$depth_change_point1),],
    {
      plot(depth_change_point1, sali_dif.se)
      plot(depth_change_point1, halocline.se)
      plot(depth_change_point1, depth_gradient.se)
      plot(depth_change_point1, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$depth_change_point2),],
    {
      plot(depth_change_point2, sali_dif.se)
      plot(depth_change_point2, halocline.se)
      plot(depth_change_point2, depth_gradient.se)
      plot(depth_change_point2, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$O2def_below_halocline),],
    {
      plot(O2def_below_halocline, sali_dif.se)
      plot(O2def_below_halocline, halocline.se)
      plot(O2def_below_halocline, depth_gradient.se)
      plot(O2def_below_halocline, depth_change_point2.se)
    })

    with(profiles[!is.na(profiles$O2def_slope_below_halocline),],
    {
      plot(O2def_slope_below_halocline, sali_dif.se)
      plot(O2def_slope_below_halocline, halocline.se)
      plot(O2def_slope_below_halocline, depth_gradient.se)
      plot(O2def_slope_below_halocline, depth_change_point2.se)
    })

    plot(profiles[!is.na(profiles$halocline),
                  c("sali_dif", "halocline", "depth_gradient",
                    "sali_dif.se", "halocline.se", "depth_gradient.se",
                    "depth_change_point1",
                    "depth_change_point2", "depth_change_point2.se")])

    plot(profiles[!is.na(profiles$depth_gradient),
                  c("sali_dif", "halocline", "depth_gradient",
                    "depth_change_point1", "depth_change_point2")])

    plot(profiles[!is.na(profiles$depth_gradient),
                     c("O2def_below_halocline", "O2def_slope_below_halocline",
                       "O2def_slope_below_halocline.se")])


    summary(profiles[!is.na(profiles$depth_gradient),
                     c("sali_surf",
                       "sali_dif", "halocline", "depth_gradient",
                       "sali_dif.se", "halocline.se", "depth_gradient.se",
                       "depth_change_point1",
                       "depth_change_point2", "depth_change_point2.se",
                       "O2def_below_halocline", "O2def_slope_below_halocline")])

  by(profiles[!is.na(profiles$depth_gradient),],
     profiles$Basin[!is.na(profiles$depth_gradient)],
     function(x) summary(x[c("sali_surf", "sali_dif", "halocline", "depth_gradient",
                           "sali_dif.se", "depth_change_point1", "depth_change_point2",
                           "depth_change_point2.se", "O2def_below_halocline",
                           "O2def_slope_below_halocline")]))

  helcom <- sf::st_read(dsn = "data/OxygenDebt/shapefiles", layer = "helcom_areas", quiet = TRUE)
  plot(sf::st_geometry(helcom), col = gplots::rich.colors(nrow(helcom), alpha = 0.5))
  plot(sf::st_geometry(makeSpatial(profiles)), add = TRUE, pch = 16, cex = 0.5)
}
