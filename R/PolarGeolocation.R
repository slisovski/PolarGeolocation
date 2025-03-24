#' Creates a calibration table
#'
#'A central function for the location estimation process. This function creates a table
#'with estimated maximum light template as well as gamma distribution parameters for the
#'recorded light across intervals (bins) of zenith angles.
#'
#'Note: Plase refer to publication in Wader Study for detailed information (Lisovski 2018).
#'
#' The \code{adjust.interval} argument can be used to specify a
#' timing adjustment for tags that report the maximum light level
#' observed in the preceeding sampling interval. If this argument is
#' zero, no adjustment is made, otherwise the timestamps of sunset
#' intervals will be adjusted by this interval to compensate for the
#' recording behaviour of the tag.
#'
#' @title Calibration function
#'
#' @param tagdata a dataframe with columns \code{Date} and
#' \code{Light} that are the sequence of sample times (as POSIXct)
#' and light levels recorded by the tag at known location.
#' @param twl a dataframe with columns \code{Twilight} and
#' \code{Rise} that are the sequence of defined sunrise and sunset times.
#' @param lon calibration longitude.
#' @param lat calibration latitude.
#' @param bin the grouping factor of light recordings across zenith angles.
#' @param adjust  timing adjustment for sunset intervals.
#' @param zAdjust adjust maximum template curve to the right.
#' @param plot logical, if \code{TRUE} a plot will be shown.
#' @return a dataframe with one row for each bin of zenith angles. Columns provide minimum and maximum zenith
#' angle (range), maximum light value (the template) and the rate and scale paramters of a fitted gamma distribution.
#' @importFrom MASS fitdistr
#' @importFrom zoo na.approx
#' @importFrom graphics par plot lines mtext polygon abline
#' @export
getTemplateCalib <- function(tagdata,
                             lon,
                             lat,
                             bin = 5,
                             adjust = 300,
                             zAdjust = 1,
                             plot = T) {


  sun <- solar(tagdata$Date)
  z   <- zenith(sun, lon, lat)

  if(!is.null(adjust)) {
    tagdata$Date[c(0,diff(z))>0] <- tagdata$Date[c(0,diff(z))>0]-adjust
    sun <- solar(tagdata$Date)
    z   <- zenith(sun, lon, lat)
  }

  top    <- max(tagdata$Light)+0.05
  bottom <- min(tagdata$Light)+0.05

  z_seq  <- seq(min(z), max(z), by = 0.1)
  dat    <- tibble(zenith = z, light = tagdata$Light)
  maxDat <- tibble(zenith = z_seq, exp = sapply(z_seq, function(x) {
    dat %>% filter(zenith>=x, zenith<x+5) %>% pull(light) %>% max()
  }))

  mid <- (top - bottom)/2
  c   <- max(maxDat[maxDat$exp>(mid-mid*0.1) &
                      maxDat$exp<(mid+mid*0.1),1]) + zAdjust

  slTab <- sapply(seq(0, 5, length = 100), function(s) {
    exp <- (top) + ((bottom - top)/(1+exp((c-z_seq)/s)))
    sum((exp - maxDat %>% pull(exp))<0)
  })
  slope <- seq(0, 5, length = 100)[which.min(slTab)]

  if(plot) {
    plot(z, tagdata$Light, pch = 16, cex = 0.5,
         col = adjustcolor('grey10', alpha.f = 0.5))
    lines(z_seq, (top) + ((bottom - top)/(1+exp((c-z_seq)/slope))),
          lwd = 2, col = "darkorange")
  }

  mL <- approxfun(x = z_seq,
                  y = (top) + ((bottom - top)/(1+exp((c-z_seq)/slope))),
                  rule = 2)


  z_seq <- round(seq(min(z), max(z)+bin, by = bin),2)
  gr    <- cut(z, z_seq, include.lowest = TRUE, labels = FALSE)

  calibTab <- cbind(z_seq[-length(z_seq)], do.call("rbind", lapply(sort(unique(gr)), function(g) {
    x <- abs(mL(z[gr==g] - tagdata$Light[gr==g]))
    if(all(x<0.1)) {
      out = cbind(1, 5)
    } else {
      fit1 <- suppressWarnings(fitdistr(x, "log-normal"))
      cbind(matrix(fit1$estimate[1:2], nrow = 1))
    }
  })))

  list(MaxL = mL,
       calibTab = calibTab[order(calibTab[,1]),])
}



#' Initialise the region
#
#' @title Initialisation
#'
#' @param tagdata a dataframe with columns \code{Date} and
#' \code{Light} that are the sequence of sample times (as POSIXct)
#' and light levels recorded by the tag at known location.
#' @param calibration ...
#' @param adjust  timing adjustment for sunset intervals.
#' @param plot logical, if \code{TRUE} a plot will be shown.
#' @return ...
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_geometry st_union st_transform st_intersection st_centroid st_as_sf st_bbox
#' @importFrom parallel makeCluster clusterEvalQ stopCluster
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom stars st_as_stars st_rasterize
#' @importFrom scales rescale
#' @export
initRegion <- function(tagdata,
                       calibration,
                       resolution = 50,
                       adjust = 300,
                       arctic = TRUE,
                       buffer = 3000,
                       mask = "land",
                       ncores = detectCores(),
                       plot = FALSE) {


  sf_use_s2(FALSE)
  polarBuffer <- 0

  if(arctic) {
     proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea", 0, 90)
  } else proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea", 0, -90)

  bbox <- st_point(c(0,0)) %>% st_sfc() %>% st_set_crs(proj) %>% st_buffer(buffer*1000)

  map <- ne_countries(scale = 50) %>%
    st_geometry() %>% st_union() %>%
    st_transform(proj) %>% st_intersection(bbox) %>% suppressMessages()

  mapRast <- st_rasterize(map %>% st_as_sf(),
                          st_as_stars(st_bbox(map), dx = resolution*1000, dy = resolution*1000, values = NA_real_))

  if(mask=="ocean") {
    mapRast <- mapRast %>% mutate(ID = ifelse(is.na(ID), 1, NA))
  }

  if(mask=="none") {
    mapRast <- mapRast %>% mutate(ID = 1) %>% st_crop(bbox)
  }

  crds0    <- st_as_sf(mapRast) %>% st_centroid() %>%
    mutate(dist = st_distance(., st_point(c(0,0)) %>% st_sfc(crs = proj))) %>% st_transform(4326) %>%
    suppressMessages() %>% suppressWarnings()

  crds     <- crds0[as.numeric(crds0$dist)>=(polarBuffer*1000),] %>%
    st_coordinates() %>% suppressMessages() %>% suppressWarnings()

  dat <- tagdata
  clb <- calibration
  adj <- adjust
  crd <- crds

  if(ncores>1) {

    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl, c('dat', 'clb', 'adj', 'crd', 'zenith', 'solar', 'refracted'), envir = environment())

    parOut <- parLapply(cl, 1:nrow(crd), function(x) {

      date <- dat$Date
      ss   <- solar(date)
      zX   <- refracted(zenith(ss, crd[x, 1], crd[x, 2]))
      ct   <- c(zX[-length(zX)]>zX[-1])
      ct   <- c(ct, ct[length(ct)])

      if(!is.null(adj)) {
        date[!ct] <- date[!ct]-adj
        ss        <- solar(date)
        zX        <- refracted(zenith(ss, crd[x, 1], crd[x, 2]))
      }

      ExpL  <- clb$MaxL(zX)
      diffL <- (ExpL -  tagdata$Light)+1e-5

      ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
      sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],
                                                        calibration$calibTab[x, 2],
                                                        calibration$calibTab[x,3], log = TRUE))))
    }) %>% Reduce("c", .)

    parallel::stopCluster(cl)

  } else {

    parOut <- lapply(1:nrow(crds), function(x) {

      date <- tagdata$Date
      ss   <- solar(date)
      zX   <- refracted(zenith(ss, crds[x, 1], crds[x, 2]))
      ct   <- c(zX[-length(zX)]>zX[-1])
      ct   <- c(ct, ct[length(ct)])

      if(!is.null(adjust)) {
        date[!ct] <- date[!ct]-adjust
        ss        <- solar(date)
        zX        <- refracted(zenith(ss, crds[x, 1], crds[x, 2]))
      }

      ExpL  <- calibration$MaxL(zX)
      diffL <- (ExpL -  tagdata$Light)+1e-5

      ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
      sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],
                                                        calibration$calibTab[x, 2],
                                                        calibration$calibTab[x,3], log = TRUE))))
    }) %>% Reduce("c", .)

  }

  rastOut <- tibble(lon = crds[,1], lat = crds[,2], p = rescale(parOut, c(1,0))) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(proj) %>%
    st_rasterize(st_as_stars(st_bbox(map), dx = resolution*1000, dy = resolution*1000, values = NA_real_))

  if(plot) {
    pl <- ggplot() +
      geom_stars(data = rastOut, show.legend = F) +
      scale_fill_binned(type = "viridis", na.value = "transparent") +
      geom_sf(data = map, fill = NA, linewidth = 0.5) +
      theme_void()
    print(pl)
  }

  list(rastOut, map, mask = mask, polarBuffer = polarBuffer)

}

#' Extract bounding box from initialisation
#
#' @title Bounding box
#'
#' @param init output from \code{initRegion}
#' @param quantile values/region to include in bounding box
#' @param plot logical, if \code{TRUE} a plot will be shown.
#' @return ...
#' @importFrom sf st_union st_transform st_intersection st_centroid st_as_sfc st_bbox
#' @importFrom dplyr mutate
#' @export
makeMask <- function(init, quantile = 0.9, buffer = 100, plot = TRUE) {

  t <- init[[1]] %>% st_as_sf() %>% st_centroid() %>% filter(!is.infinite(p)) %>% st_transform(4326) %>%
    mutate(lon      = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%  st_drop_geometry() %>%
    arrange(lon) %>% mutate(r  = rollapply(p, 5, max, fill=NA, align="right")) %>%
    arrange(desc(r)) %>% suppressWarnings()

  proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea",
                  t$lon[1], min(abs(t$lat))) %>% suppressWarnings()

  quants <- init[[1]] %>% st_as_sf() %>% filter(!is.infinite(p)) %>%
    filter(p > quantile(p, probs = quantile)) %>% st_union() %>%
    st_buffer(buffer*1000) %>% st_transform(proj)

  box <- quants %>% st_transform(proj) %>% st_bbox() %>% st_as_sfc(crs = proj)
  map <- init[[2]] %>% st_transform(proj) %>% st_intersection(box)
  col <- init[[1]] %>% st_as_sf() %>% filter(!is.infinite(p)) %>% st_transform(proj) %>% st_intersection(box) %>% suppressWarnings()

  if(plot) {
    pl <- ggplot() +
      geom_sf(data = map, fill = "grey90", color = "grey10", show.legend = F) +
      geom_sf(data = col, mapping = aes(fill = p), show.legend = F, alpha = 0.7) +
      scale_fill_binned(type = "viridis", na.value = "transparent") +
      geom_vline(xintercept = t$lon[1], linewidth = 1, color = "orange", linetype = 2) +
      geom_sf(data = box, fill = NA, color = "black") +
      theme_void()
    print(pl)
  }

  list(bbox = box, map = map, mask = init$mask, polarBuffer = init$polarBuffer)
}

#' Template fit
#
#' @title Template fit
#'
#' @export
templateEstimate <- function(tagdata,
                             calibration,
                             bbox,
                             resolution = 15,
                             adjust = 300,
                             ncores = detectCores(),
                             plot = TRUE) {

  sf_use_s2(FALSE)

  mapRast <- st_rasterize(bbox$map %>% st_as_sf(), st_as_stars(st_bbox(bbox$map),
                                                               dx = resolution*1000, dy = resolution*1000, values = NA_real_))

  crds  <- mapRast %>%
    st_as_sf(na.rm = ifelse(bbox$mask=='none', FALSE, TRUE)) %>% st_centroid() %>% st_transform(4326) %>%
    mutate(dist = st_distance(., st_point(c(0,90)) %>% st_sfc(crs = 4326))) %>%
    filter(as.numeric(dist) >= as.numeric((bbox$polarBuffer*1000))) %>%
    st_coordinates() %>%
    suppressMessages() %>% suppressWarnings()

  cl <- makeCluster(getOption("cl.cores", ncores))
  clusterExport(cl, c('tagdata', 'calibration', 'crds', 'zenith', 'solar', 'refracted', 'adjust'), envir = environment())

  pOut <- unlist(parLapply(cl, 1:nrow(crds), function(p) {

    date <- tagdata$Date
    ss   <- solar(date)
    zX   <- refracted(zenith(ss, crds[p,1], crds[p,2]))
    ct   <- c(zX[-length(zX)]>zX[-1])
    ct   <- c(ct, ct[length(ct)])

    if(!is.null(adjust)) {
      date[!ct] <- date[!ct]-adjust
      ss        <- solar(date)
      zX        <- refracted(zenith(ss, crds[p,1],  crds[p,2]))
    }

    ExpL  <- calibration$MaxL(zX)
    diffL <- ExpL - tagdata$Light

    ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)

    sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],
                                                      calibration$calibTab[x,2],
                                                      calibration$calibTab[x,3], log = TRUE))))
  }))

  stopCluster(cl)

  rastOut <- tibble(lon = crds[,1], lat = crds[,2], p = pOut) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(mapRast)) %>%
    st_rasterize(st_as_stars(st_bbox(bbox$map), dx = (resolution+2)*1000, dy = (resolution+2)*1000, values = NA_real_)) %>%
    suppressWarnings() %>% suppressMessages()

  cont    <- st_contour(rastOut, breaks = quantile(rastOut[[1]], probs = seq(0, 0.4, length = 10), na.rm = T))

  if(plot) {
    pl <- ggplot() +
      geom_sf(data = bbox$map, fill = "grey90", color = "grey10", show.legend = F) +
      geom_stars(data = rastOut, alpha = 0.7, show.legend = F) +
      scale_fill_viridis_c(na.value = NA, direction = -1) +
      theme_light()
    if(contour) pl <- pl + geom_sf(data = cont, fill = NA, color = "orange") %>% suppressMessages()
    print(pl) %>% suppressWarnings()
  }


  list(fit = rastOut, map = bbox$map)
}



#' Extract Location
#
#' @title Extract breeding site
#'
#' @export
locationSummary <- function(fit, quantile = 0.8, plot = T) {

  tmp <- fit$fit %>% st_as_sf() %>% st_centroid() %>% st_transform(4326) %>%
    mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% suppressWarnings()

  lon  <- tmp[tmp$p < quantile(p, probs = quantile),] %>% pull(lon) %>% quantile(., probs = c(0,0.5,1))
  lat <- tmp[tmp$p < quantile(p, probs = quantile),] %>% pull(lat) %>% quantile(., probs = c(0,1))

  t1 <- tibble(mask = 'none', lon = lon[2], lat = lat[1])

  out <- t1 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(fit$map))

  ### land
  pts <- st_as_sf(tibble(lon = seq(-180, 180, length = 360*2),
                         lat = as.numeric(t1$lat)), coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(st_crs(fit$map)) %>% st_intersection(fit$map) %>%
    mutate(dist = st_distance(., out)) %>% arrange(dist) %>% st_transform(4326) %>% st_coordinates()

  t2 <- tibble(mask = 'mask', lon = pts[1,1], lat = pts[1,2])

  out <- rbind(out, t2 %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(fit$map)))

  if(plot) {
    pl <- ggplot() +
      geom_sf(data = fit$map, fill = "grey90", color = "grey10", show.legend = F) +
      geom_stars(data = fit$fit, alpha = 0.6, show.legend = F) +
      scale_fill_viridis_c(na.value = NA, direction = -1) +
      geom_sf(data = out, mapping = aes(color = mask), shape = 21, fill = "white", size = 5, linewidth = 10) +
      theme_light()
    print(pl) %>% suppressWarnings()
  }


  rbind(t1, t2)

}

## Functions copied from SGAT and temporariliy included into PolarGeolocation:
## 1) solar
## 2) zenith
## 3) refracted
## 4) twilight
## 5) geolight.convert



##' Calculate solar time, the equation of time and solar declination
##'
##' The solar time, the equation of time and the sine and cosine of
##' the solar declination are calculted for the times specified by
##' \code{tm} using the same methods as
##' \url{https://gml.noaa.gov/grad/solcalc/}.
##' @title Solar Time and Declination
##' @param tm a vector of POSIXct times.
##' @return A list containing the following vectors.
##' \item{\code{solarTime}}{the solar time (degrees)}
##' \item{\code{eqnTime}}{the equation of time (minutes of time)}
##' \item{\code{sinSolarDec}}{sine of the solar declination}
##' \item{\code{cosSolarDec}}{cosine of the solar declination}
##' @seealso \code{\link{zenith}}
##' @examples
##' ## Current solar time
##' solar(Sys.time())
##' @export
solar <- function(tm) {

  rad <- pi/180

  ## Time as Julian day (R form)
  Jd <- as.numeric(tm)/86400.0+2440587.5

  ## Time as Julian century [G]
  Jc <- (Jd-2451545)/36525

  ## The geometric mean sun longitude (degrees) [I]
  L0 <- (280.46646+Jc*(36000.76983+0.0003032*Jc))%%360


  ## Geometric mean anomaly for the sun (degrees) [J]
  M <- 357.52911+Jc*(35999.05029-0.0001537*Jc)

  ## The eccentricity of earth's orbit [K]
  e <- 0.016708634-Jc*(0.000042037+0.0000001267*Jc)

  ## Equation of centre for the sun (degrees) [L]
  eqctr <- sin(rad*M)*(1.914602-Jc*(0.004817+0.000014*Jc))+
    sin(rad*2*M)*(0.019993-0.000101*Jc)+
    sin(rad*3*M)*0.000289

  ## The true longitude of the sun (degrees) [M]
  lambda0 <- L0 + eqctr

  ## The apparent longitude of the sun (degrees) [P]
  omega <- 125.04-1934.136*Jc
  lambda <- lambda0-0.00569-0.00478*sin(rad*omega)


  ## The mean obliquity of the ecliptic (degrees) [Q]
  seconds <- 21.448-Jc*(46.815+Jc*(0.00059-Jc*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60

  ## The corrected obliquity of the ecliptic (degrees) [R]
  omega <- 125.04-1934.136*Jc
  obliq <- obliq0 + 0.00256*cos(rad*omega)

  ## The equation of time (minutes of time) [U,V]
  y <- tan(rad*obliq/2)^2
  eqnTime <- 4/rad*(y*sin(rad*2*L0) -
                      2*e*sin(rad*M) +
                      4*e*y*sin(rad*M)*cos(rad*2*L0) -
                      0.5*y^2*sin(rad*4*L0) -
                      1.25*e^2*sin(rad*2*M))

  ## The sun's declination (radians) [T]
  solarDec <- asin(sin(rad*obliq)*sin(rad*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)

  ## Solar time unadjusted for longitude (degrees) [AB!!]
  ## Am missing a mod 360 here, but is only used within cosine.
  solarTime <- ((Jd-0.5)%%1*1440+eqnTime)/4
  #solarTime <- ((Jd-2440587.5)*1440+eqnTime)/4

  ## Return solar constants
  list(solarTime=solarTime,
       eqnTime=eqnTime,
       sinSolarDec=sinSolarDec,
       cosSolarDec=cosSolarDec)
}



##' Calculate the solar zenith angle for given times and locations
##'
##' \code{zenith} uses the solar time and declination calculated by
##' \code{solar} to compute the solar zenith angle for given times and
##' locations, using the same methods as
##' \url{https://gml.noaa.gov/grad/solcalc/}.  This function does not
##' adjust for atmospheric refraction see \code{\link{refracted}}.
##' @title Solar Zenith Angle
##' @param sun list of solar time and declination computed by \code{solar}.
##' @param lon vector of longitudes.
##' @param lat vector latitudes.
##' @return A vector of solar zenith angles (degrees) for the given
##' locations and times.
##' @seealso \code{\link{solar}}
##' @examples
##' ## Approx location of Sydney Harbour Bridge
##' lon <- 151.211
##' lat <- -33.852
##' ## Solar zenith angle for noon on the first of May 2000
##' ## at the Sydney Harbour Bridge
##' s <- solar(as.POSIXct("2000-05-01 12:00:00","EST"))
##' zenith(s,lon,lat)
##' @export
zenith <- function(sun,lon,lat) {

  rad <- pi/180

  ## Suns hour angle (degrees) [AC!!]
  hourAngle <- sun$solarTime+lon-180
  #hourAngle <- sun$solarTime%%360+lon-180

  ## Cosine of sun's zenith [AD]
  cosZenith <- (sin(rad*lat)*sun$sinSolarDec+
                  cos(rad*lat)*sun$cosSolarDec*cos(rad*hourAngle))

  ## Limit to [-1,1] [!!]
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1

  ## Ignore refraction correction
  acos(cosZenith)/rad
}


##' Adjust the solar zenith angle for atmospheric refraction.
##'
##' Given a vector of solar zeniths computed by \code{\link{zenith}},
##' \code{refracted} calculates the solar zeniths adjusted for the
##' effect of atmospheric refraction.
##'
##' \code{unrefracted} is the inverse of \code{refracted}. Given a
##' (single) solar zenith adjusted for the effect of atmospheric
##' refraction, \code{unrefracted} calculates the solar zenith as
##' computed by \code{\link{zenith}}.
##'
##' @title Atmospheric Refraction
##' @param zenith zenith angle (degrees) to adjust.
##' @return vector of zenith angles (degrees) adjusted for atmospheric
##' refraction.
##' @examples
##' ## Refraction causes the sun to appears higher on the horizon
##' refracted(85:92)
##' ## unrefracted gives unadjusted zenith (see SGAT)
##'
##' @export
refracted <- function(zenith) {
  rad <- pi/180
  elev <- 90-zenith
  te <- tan((rad)*elev)
  ## Atmospheric Refraction [AF]
  r <- ifelse(elev>85,0,
              ifelse(elev>5,58.1/te-0.07/te^3+0.000086/te^5,
                     ifelse(elev>-0.575,
                            1735+elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711))),-20.772/te)))
  ## Corrected Zenith [90-AG]
  zenith-r/3600
}
