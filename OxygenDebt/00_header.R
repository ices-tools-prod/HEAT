# ----------------------------
#
#   load packages and clear work space
#
# ----------------------------

header <-
function(phase = c("data", "input", "model", "output"))
{
  phase <- match.arg(phase)
  # phase specific packages -------------------
  pkg <-
    switch(phase,
      #data = c("sp", "rgdal", "rgeos", "dplyr"),
      data = c("sf", "dplyr"),
      input = c("lubridate"),
      model = c("stats", "survival", "dplyr", "mgcv", "sp", "rgdal"),
      output = c("tidyr", "dplyr")
    )

  # common packages -------------------
  pkg <- unique(c(pkg, "viridis", "gplots", "jsonlite"))

  # install dependencies used in the analysis
  newpkg <- setdiff(pkg, installed.packages())
  if (length(newpkg)) install.packages(newpkg, dependencies = TRUE, quiet = TRUE)

  # always load certain packages -------------------
  always_load <- c("mgcv")
  if (always_load %in% pkg) {
    tmp <- sapply(intersect(always_load, pkg), library, character.only = TRUE, quietly = TRUE, pos = 3)
  }

  # load utils into oxydebt_funs namespace
  tmp <-
    sapply(dir("OxygenDebt/Utils"),
           function(x) {
           sys.source(paste0("OxygenDebt/Utils/", x), envir = environment(header))
           })

  # clear workspace -------------------
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  
  # Define assessment period i.e. uncomment the period you want to run the assessment for!
  #assessmentPeriod <<- "2011-2016" # HOLAS II
  assessmentPeriod <<- "2016-2021" # HOLAS III
}




