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
#' @param plot logical, if \code{TRUE} a plot will be shown.
#' @return a dataframe with one row for each bin of zenith angles. Columns provide minimum and maximum zenith
#' angle (range), maximum light value (the template) and the rate and scale paramters of a fitted gamma distribution.
#' @importFrom MASS fitdistr
#' @importFrom zoo na.approx
#' @importFrom graphics par plot lines mtext polygon abline
#' @export
getTemplateCalib <- function(tagdata,
                             twl,
                             lon,
                             lat,
                             bin = 1,
                             adjust = 300,
                             max.light.fit = NULL,
                             plot = TRUE) {

  ### first split into segments
  sun <- solar(tagdata$Date)
  z   <- zenith(sun, lon, lat)
  ct  <- c(z[-length(z)]>z[-1])
  ct <- c(ct, ct[length(ct)])

  ### second split with adjusted twilights
  if(!is.null(adjust)) {
    tagdata$Date[!ct] <- tagdata$Date[!ct]-adjust
    sun <- solar(tagdata$Date)
    z   <- zenith(sun, lon, lat)
    ct  <- c(ct, ct[length(ct)])
  }

  ind <- c(1, which(ct[-length(ct)]!=ct[-1]), length(z))

  segs <- NA
  segs[ind+1] <- 1:length(ind)
  segs <- na.approx(segs, method = "constant", rule = 2)

  tagdata$Elev <- z
  tagdata$Rise <- ifelse(ct[-1], TRUE, FALSE)
  tagdata$Segment <- na.approx(segs[-length(segs)], method = "constant", rule = 2)
  tagdata <- subset(tagdata, !is.na(Rise))


  tmp03 <- split(cbind(ID = 1:nrow(tagdata), tagdata), f = tagdata$Segment)

  bin.seq0 <- seq(min(tagdata$Elev), max(tagdata$Elev), by = bin)
  bin.seq  <- cbind(bin.seq0[-length(bin.seq0)], bin.seq0[-1])

  out1 <- apply(bin.seq, 1, function(x1) {
    as.numeric(unlist(lapply(tmp03, function(x2) subset(x2, Elev>=x1[1] & Elev<=x1[2], select = ID))))
  })


  out2_1 <- lapply(out1, function(x) max(tagdata[x,"Light"]))



  ###################
  #### Inser new ####
  ###################
  bin0 <- seq(bin-round((bin*0.25),1), bin+round((bin*0.25),1), by = 0.1)


  maxTab <- cbind(z = NA, l = NA)
  for(i in 1:length(bin0)) {
    bin.seq0 <- seq(min(tagdata$Elev), max(tagdata$Elev), by = bin0[i])
    bin.seq  <- cbind(bin.seq0[-length(bin.seq0)], bin.seq0[-1])

    out1 <- apply(bin.seq, 1, function(x1) {
      as.numeric(unlist(lapply(tmp03, function(x2) subset(x2, Elev>=x1[1] & Elev<=x1[2], select = ID))))
    })

    maxTmp   <- lapply(out1, function(x) max(tagdata[x,"Light"]))
    maxTab   <- rbind(maxTab, cbind(z = apply(bin.seq, 1, mean), l = unlist(maxTmp)))
  }

  if(is.null(max.light.fit)) max.light.fit <- max(maxTab[,2], na.rm = T)

  t1 <- maxTab[maxTab[,1]>max(maxTab[which(maxTab[,2]==max(tagdata$Light)),1]) &
               maxTab[,2]>min(tagdata$Light) &
               maxTab[,1]<quantile(tagdata$Elev[tagdata$Light==min(tagdata$Light)], probs = 0.1) &
               maxTab[,2]<max.light.fit  , ]
  mod <- lm(l~z, data = as.data.frame(t1))
     # points(as.data.frame(t1), pch = 16)
     # abline(mod)

  zNew <- seq((max(tagdata$Light) -coefficients(mod)[1])/coefficients(mod)[2],
              (min(tagdata$Light) -coefficients(mod)[1])/coefficients(mod)[2], length = 50)

  mL   <- approxfun(x = zNew, y = predict(mod, newdata = data.frame(z = zNew)), rule = 3)



  out2_2 <- apply(cbind(1:length(out1)), 1, function(x) apply(tagdata[out1[[x]],c("Light", "Elev")], 1, function(y) abs((mL(y[2])+0.01)-y[1])))
  med    <- apply(cbind(1:nrow(bin.seq)), 1, function(x) median(tagdata[out1[[x]],"Light"]))

  out3 <- lapply(out2_2, function(x) {
    if(all(x==0.01)) {
      out = cbind(1, 5)
    } else {
      fit1 <- suppressWarnings(fitdistr(x, "log-normal"))
      cbind(matrix(fit1$estimate[1:2], nrow = 1))
    }
  })

  slp <- cbind(apply(bin.seq, 1, mean), c(NA, NA, as.numeric(rollapply(med, 5, function(x) coefficients(lm(x~c(1:length(x))))[2])), NA, NA))
  cut <- max(slp[which(slp[-nrow(slp),2]> -0.25 & slp[-1,2]<= -0.25),1])

  if(plot) {
    colSet  <- colorRampPalette(c("darkblue", "aliceblue"))
    colRise <- colorRampPalette(c("darkred", "rosybrown1"))

    opar <- par(mar = c(6,6,1,1), cex.lab = 1.2, las = 1)
    plot(NA, xlim = c(min(z), 110), ylim = c(0, max(tagdata$Light, na.rm = T)), xlab = "Zenith", ylab = "Light")

    x1 <- seq(min(z), 110, length = 500)
    y1 <- mL(x1)
    polygon(c(x1, rev(x1)),
            c(y1, rep(max(tagdata$Light)+0.3, length(x1))),
            col = adjustcolor("orange", alpha.f = 0.4), border = "grey60", lty = 2, lwd = 2)

    abline(v = bin.seq0, col = adjustcolor("grey30", alpha.f = 0.4), lty = 3)

    t <- lapply(tmp03, function(x) {
      if(x[1,5]) {
        cl <- colRise(length(unique(tagdata$Segment[tagdata$Rise])))[unique(x[,"Segment"])]
      } else {
        cl <- colSet(length(unique(tagdata$Segment[!tagdata$Rise])))[unique(x[,"Segment"])]
      }
      lines(x[x[,4]<110,4], x[x[,4]<110,3], type = "o", pch = 16, col = cl, cex = 0.75)
    })

    lines(c(cut, cut), c(0, max(tagdata$Light, na.rm = T)+0.3), lwd = 3, lty = 4)
    par(opar)
  }


  list(MaxL = mL, calibTab = cbind(apply(bin.seq, 1, mean), do.call("rbind", out3)), cut = cut)

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
