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
#' @importFrom GeoLight solar zenith refracted
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


#' Initialise mask
#'
#'
initMask <- function(tagdata,
                     calibration,
                     res = 50,
                     adjust = 300,
                     quantile = 0.85,
                     hemisphere = "north",
                     mask = "land",
                     plot = FALSE) {

  require(sf)
  sf_use_s2(FALSE)
  require(stars)
  require(dplyr)
  require(ggplot2)

  if(is.null(map)) {
    if(hemisphere=="north") {
        proj_start <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea", 0, 90)
    } else {
        proj_start <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea", 0, -90)
    }
    bbox <- st_point(c(0,0)) %>% st_sfc() %>% st_set_crs(proj_start) %>% st_buffer(30e5)

    map <- rnaturalearth::ne_countries(scale = 50) %>%
      st_geometry() %>% st_union() %>%
      st_transform(proj_start) %>% st_intersection(bbox) %>% suppressMessages()
  }

  mapRast <- st_rasterize(map %>% st_as_sf(), st_as_stars(st_bbox(map), dx = res*1000, dy = res*1000, values = NA_real_))

  if(mask!="land") {
    mapRast <- mapRast %>% mutate(ID = ifelse(is.na(ID), 1, NA))
  }

  crds    <- st_as_sf(mapRast, na.rm = FALSE) %>% st_centroid() %>% st_transform(4326) %>%
    st_coordinates() %>% suppressMessages() %>% suppressWarnings()

  mapRast <- mapRast %>% mutate(lon = crds[,1], lat = crds[,2]) %>% merge(names = 'layers')

  testOut <- st_apply(mapRast, 1:2, function(x) {

    if(!is.na(x[1])) {

    date <- tagdata$Date
    ss   <- solar(date)
    zX   <- refracted(zenith(ss, x[2], x[3]))
    ct   <- c(zX[-length(zX)]>zX[-1])
    ct   <- c(ct, ct[length(ct)])

    if(!is.null(adjust)) {
      date[!ct] <- date[!ct]-adjust
      ss        <- solar(date)
      zX        <- refracted(zenith(ss, x[2], x[3]))
    }

    ExpL  <- calibration$MaxL(zX)
    diffL <- (ExpL -  tagdata$Light)+1e-5

    ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
    sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],  calibration$calibTab[x, 2], calibration$calibTab[x,3], log = FALSE))))

    } else NA

  }, FUTURE = TRUE)

  cont <-  testOut %>% setNames("probs") %>% mutate(probs = probs >= quantile(probs, probs = quantile, na.rm = T)) %>%
    st_as_sf() %>% filter(probs) %>% st_union()

  if(plot) {
  pl <- ggplot() +
    geom_stars(data = testOut, show.legend = F) +
    scale_fill_binned(type = "viridis", na.value = "transparent") +
    geom_sf(data = map, fill = NA, linewidth = 0.5) +
    geom_sf(data = cont, fill = NA, color = "red", linewidth = 0.85) +
    theme_void()
  print(pl)
  }

  cont

}


#' Location Estimation
#'
templateEstimate <- function(tagdata,
                             calibration,
                             bbox = bbox,
                             res = 15,
                             window = 24,
                             adjust = 300,
                             quantile = c(0.4, 0.6),
                             mask = "land",
                             cutoff = 3,
                             parallel = TRUE,
                             ncores = 4,
                             path = FALSE,
                             plot = TRUE) {

    require(sf)
    sf_use_s2(FALSE)
    require(stars)
    require(dplyr)
    require(ggplot2)

    bbox_out <- bbox %>% st_bbox %>% st_as_sfc(st_set_crs(bbox))
    proj     <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea",
                        (bbox_out %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,1],
                        (bbox_out %>% st_centroid() %>% st_transform(4326) %>% st_coordinates())[,2])

    map <- rnaturalearth::ne_countries(scale = 50) %>%
        st_geometry() %>% st_union() %>%
        st_transform(proj) %>% st_intersection(bbox_out %>% st_transform(proj)) %>% suppressMessages()

    mapRast <- st_rasterize(map %>% st_as_sf(), st_as_stars(st_bbox(map), dx = res*1000, dy = res*1000, values = NA_real_))

    if(mask!="land") {
      mapRast <- mapRast %>% mutate(ID = ifelse(is.na(ID), 1, NA))
    }

    crds <- mapRast %>% st_as_sf() %>% st_centroid() %>% st_transform(4326) %>% st_coordinates() %>% suppressWarnings()

    dat  <- tagdata %>% as_tibble()
    ss   <- solar(dat$Date)
    zX   <- refracted(zenith(ss,
                (bbox_out %>% st_transform(proj) %>% st_centroid() %>% st_centroid() %>% st_coordinates())[,1],
                (bbox_out %>% st_transform(proj) %>% st_centroid() %>% st_centroid() %>% st_coordinates())[,2]))
    dat  <- dat %>% mutate(ct = c(c(zX[-length(zX)]>zX[-1]), TRUE),
                           Date = as.POSIXct(ifelse(!ct, Date - adjust, Date), origin = "1970-01-01")) %>%
            group_split(Window = cut(as.numeric((Date - Date[1])/60/60), window*c(0:100), labels = FALSE, include.lowest = T))


    if(parallel) {

        require(parallel)
        cl <- makeCluster(getOption("cl.cores", ncores))
        clusterExport(cl, c('crds', "calibration"))
        clusterEvalQ(cl, {
          library(SGAT)

        })

        test <- parLapply(cl, dat, function(x) {
            p_log <- apply(crds, 1, function(p) {
                          zX    <- refracted(zenith(solar(x$Date), p[1], p[2]))
                          ExpL  <- calibration$MaxL(zX)
                          diffL <- (ExpL -  x$Light)+1e-5

                          ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
                          sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],  calibration$calibTab[x, 2], calibration$calibTab[x,3], log = FALSE))))
                        })
            p_log
        }) %>% Reduce("cbind", .)
    } else {
      test <- lapply(dat, function(x) {
        p_log <- apply(crds, 1, function(p) {
          zX    <- refracted(zenith(solar(x$Date), p[1], p[2]))
          ExpL  <- calibration$MaxL(zX)
          diffL <- (ExpL -  x$Light)+1e-5

          ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
          sum(unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],  calibration$calibTab[x, 2], calibration$calibTab[x,3], log = FALSE))))
        })
        p_log
      }) %>% Reduce("cbind", .)
    }

    if(parallel) stopCluster(cl)

    maxPts <- as_tibble(cbind(crds, p = apply(test[,-c(1:cutoff, (ncol(test)-cutoff+1):ncol(test))], 1, max))) %>%
      setNames(c("X", "Y", "p"))

    if(plot) {
      pl <- ggplot() +
        geom_point(data = tibble(X = crds[,1],
                                 Y = crds[,2],
                                 p_log = apply(test[,-c(1:cutoff, (ncol(test)-cutoff+1):ncol(test))], 1, sum)), aes(x = X, y = Y, fill = p_log),
                   shape = 21, color = "transparent", show.legend = FALSE, size = 2) +
        scale_fill_viridis_c() +
        if(path) geom_line(data = maxPts[apply(test, 2, which.max),], mapping = aes(x = X, y = Y)) +
        if(path) geom_point(data = maxPts[apply(test, 2, which.max),], mapping = aes(x = X, y = Y, size = p), fill = "transparent") +
        theme_light()
      print(pl)
    }


    outList <- maxPts %>% filter(p>=quantile(p, probs = quantile)) %>% dplyr::select(X, Y) %>%
      apply(2, function(x) list(lon = median(x), lower = quantile(x, probs = quantile[1]), upper = quantile(x, probs = quantile[2])))

    tibble(lon = outList$X[[1]], lat = outList$Y[[1]],
           lon_lower = outList$X[[2]], lon_upper = outList$X[[3]],
           lat_lower = outList$Y[[2]], lat_upper = outList$Y[[3]])
}






#'
#' #' This function uses the calibration table to evaluate the likelihood that the light has been recorded at any location of the grid.
#' #'
#' #' @title Calculation of spatial likelihood surface
#' #'
#' #' @param tagdata a dataframe with columns \code{Date} and
#' #' \code{Light} that are the sequence of sample times (as POSIXct)
#' #' and light levels recorded by the tag at known location.
#' #' @param calibration calibration table created by \code{getTemplateCalib}.
#' #' @param window number of days (24h) to put together for likelihood estimation (see details).
#' #' @param mask set mask for areas animals are assumed to be restricted to, "sea", "land", or "none"
#' #' @param useMask logical, if \code{TRUE}, the mask is not only used to define the spatial extent but the pre-defined likelihoods (e.g. land vs. ocean) will be taken into account.
#' #' @param cores parallel computing if >1.
#' #' @param message logical, if TRUE messages reporting the progress will be printed in the console.
#' #' @importFrom GeoLight solar refracted zenith
#' #' @importFrom stats dlnorm
#' #' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ parApply stopCluster
#' #'
#' #' @export
#' templateEstimate <- function(tagdata,
#'                              calibration,
#'                              window = 1,
#'                              adjust = 300,
#'                              mask,
#'                              useMask = FALSE,
#'                              cores = detectCores()-1,
#'                              message = TRUE) {
#'
#'   # Define segment by date
#'   seg  <- floor((as.numeric(tagdata$Date)- as.numeric(min(tagdata$Date)))/(24*60*60))
#'   nseg <- length(unique(seg))
#'
#'   ind <- cbind(c(1:nseg)[-c((nseg-(window-1)):nseg)], c(1:nseg)[-c(1:window)])
#'
#'   # Split into `slices`
#'   slice <- apply(ind, 1, function(x) tagdata[seg%in%c(x[1]:x[2]),])
#'
#'   n   <- length(slice)
#'   if(useMask) pts <- mask[[1]][mask[[1]][,3]==1,1:2] else pts <- mask[[1]]
#'
#'
#'   logp <- function(c, w, adjust) {
#'
#'     date <- slice[[w]]$Date
#'
#'     ### first split into segments
#'     ss  <- solar(date)
#'     zX  <- refracted(zenith(ss, c[1L], c[2L]))
#'     ct  <- c(zX[-length(zX)]>zX[-1])
#'     ct  <- c(ct, ct[length(ct)])
#'
#'     ### second split with adjusted twilights
#'     if(!is.null(adjust)) {
#'       date[!ct] <- date[!ct]-adjust
#'       ss        <- solar(date)
#'       zX        <- refracted(zenith(ss, c[1L], c[2L]))
#'     }
#'
#'     ExpL  <- calibration$MaxL(zX)
#'     diffL <- (ExpL -  slice[[w]]$Light)+1e-5
#'
#'     # plot(slice[[w]]$Date, slice[[w]]$Light, type = "o", pch = 16, col = adjustcolor("grey", alpha.f = 0.5))
#'     # points(slice[[w]]$Date, ExpL, type = "o", cex = 0.5)
#'     # par(new = TRUE)
#'     # plot(date, diffL, type = "l", xaxt = "n", xaxt = "n", lwd = 2, col = "orange")
#'     # axis(4)
#'
#'     ind <- cut(zX, breaks = calibration$calibTab[,1], labels = FALSE)
#'     out1  <- unlist(sapply(unique(ind), function(x) dlnorm(diffL[ind==x],  calibration$calibTab[x, 2], calibration$calibTab[x,3], log = FALSE)))
#'
#'     cbind(sum(out1, na.rm = T), ifelse(all(zX<=calibration$cut), 100, sum(diffL<0)))
#'   }
#'
#'
#'   if(message) cat("making cluster\n")
#'   mycl <- makeCluster(cores)
#'   tmp  <- clusterSetRNGStream(mycl)
#'   tmp  <- clusterExport(mycl,c("slice", "calibration", "logp"), envir=environment())
#'   tmp  <- clusterEvalQ(mycl, library("GeoLight"))
#'
#'   ## Compute likelihood
#'   ll <- array(dim = c(nrow(pts), n, 2))
#'
#'   for(i in 1:n) {
#'
#'     if(message) {
#'       cat("\r", "window", i, " of ", n)
#'       flush.console()
#'     }
#'
#'     ll[,i,] <- t(parApply(mycl, pts, 1, FUN = logp, w = i, adjust = adjust))
#'
#'   }
#'
#'   end <- stopCluster(mycl)
#'
#'   list(crds = pts, probTab = ll)
#'
#' }
#'
#'
#'
#'
#'
#'
#' templateSummary <- function(tempEst,
#'                             mask,
#'                             probs = c(0.05, 0.1),
#'                             cutoff = 1) {
#'
#'
#'   maskP   <- project(mask[[1]][,-3], proj = proj4string(mask$map))
#'     maskR <- rasterize(maskP, mask$grid, field = mask[[1]][,3])
#'
#'   crdsLL <- project(tempEst$crds[,1:2], proj = proj4string(mask$map))
#'
#'     mLik   <- apply(tempEst$probTab[,,1], 1, function(x) sum(x))
#'     ind    <- apply(tempEst$probTab[,,2], 1, function(x) sum(x)>cutoff)
#'       tt     <- rasterize(crdsLL[!ind,], mask$grid, field = mLik[!ind])
#'       tt[]   <- values(tt)/max(values(tt), na.rm = T)
#'
#'     ttU <- tt
#'       ttU[!maskR[]] <- NA
#'     ttM <- tt
#'       ttM[!is.na(ttU)] <- NA
#'
#'
#'
#'   crd0 <- coordinates(tt)[which.max(ttU[]),]
#'
#'   crd1 <- coordinates(tt)[which(tt[]>=(1-probs[2])),]
#'   crd2 <- coordinates(tt)[which(tt[]>=(1-probs[1])),]
#'
#'
#'   invCrd <- project(matrix(c(crd0, apply(crd1,2,min), apply(crd1,2,max), apply(crd2,2,min), apply(crd2,2,max)), ncol = 2, byrow = T),
#'                     proj =  proj4string(mask$map), inv = T)
#'
#'   out <- data.frame(t(as.vector(t(invCrd))[c(1,2,3,7,9,5,4,6,8,10)]))
#'     names(out) <- c("Lon", "Lat", "Lon.lower1", "Lon.lower2", "Lon.upper2", "Lon.upper1", "Lat.lower1", "Lat.lower2", "Lat.upper2", "Lat.upper1")
#'
#'
#'     opar <- par(mar = c(0,2,0,0), bty = "n")
#'     brks <- seq(min(tt[], na.rm = T), max(tt[], na.rm = T), length = 100)
#'
#'     plot(ttU, legend = FALSE, breaks = brks, col = rev(rainbow(100, start = 0, end = 0.7)),
#'          xaxt = "n", yaxt = "n")
#'     plot(ttM, legend = FALSE, breaks = brks, col = rev(rainbow(100, start = 0, end = 0.7, s = 0.2)),
#'          xaxt = "n", yaxt = "n", add = T)
#'
#'     plot(mask$map, add = T)
#'     arrows(crd0[1], min(crd1[,2]), crd0[1], max(crd1[,2]), length = 0, lwd = 0.7)
#'     arrows(crd0[1], min(crd2[,2]), crd0[1], max(crd2[,2]), length = 0, lwd = 3)
#'
#'     arrows(min(crd1[,1]), crd0[2], max(crd1[,1]), crd0[2], length = 0, lwd = 0.7)
#'     arrows(min(crd2[,1]), crd0[2], max(crd2[,1]), crd0[2], length = 0, lwd = 3)
#'
#'     points(crd0[1],  crd0[2], pch = 21, cex = 3, lwd = 2, bg = "white")
#'     par(opar)
#'
#'
#'   out
#' }
#'
#'
