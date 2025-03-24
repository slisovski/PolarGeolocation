# remotes::install_github("slisovski/PolarGeolocation", force = TRUE)
library(PolarGeolocation)

data("GRKN_calib")
data("GRKN_twlCalib")
data("GRKN_breed")

lon = 122.27
lat = -17.96

lux.calib <- GRKN_calib
twl       <- GRKN_twlCalib

tagdata   <- GRKN_breed

calib <- getTemplateCalib(tagdata = lux.calib,
                          lon = lon,
                          lat = lat,
                          zAdjust = 0,
                          bin = 5)

bboxRast <- initRegion(tagdata = GRKN_breed,
                       calibration = calib,
                       resolution = 100,
                       mask = "none",
                       plot = TRUE)

bbox     <- makeMask(init = bboxRast, quantile = 0.85, plot = TRUE)

fit <- templateEstimate(tagdata = GRKN_breed,
                        calibration = calib,
                        bbox = bbox,
                        resolution = 50,
                        adjust = 300,
                        contour = TRUE,
                        ncores = detectCores(),
                        plot = TRUE)


locationSummary(fit, quantile = 0.85, plot = T)
