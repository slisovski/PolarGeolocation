% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolarGeolocation.R
\name{templateEstimate}
\alias{templateEstimate}
\title{Calculation of spatial likelihood surface}
\usage{
templateEstimate(tagdata, calibration, window = 1, adjust = 300, mask,
  useMask = FALSE, cores = detectCores() - 1, message = TRUE)
}
\arguments{
\item{tagdata}{a dataframe with columns \code{Date} and
\code{Light} that are the sequence of sample times (as POSIXct)
and light levels recorded by the tag at known location.}

\item{calibration}{calibration table created by \code{getTemplateCalib}.}

\item{window}{number of days (24h) to put together for likelihood estimation (see details).}

\item{mask}{set mask for areas animals are assumed to be restricted to, "sea", "land", or "none"}

\item{useMask}{logical, if \code{TRUE}, the mask is not only used to define the spatial extent but the pre-defined likelihoods (e.g. land vs. ocean) will be taken into account.}

\item{cores}{parallel computing if >1.}

\item{message}{logical, if TRUE messages reporting the progress will be printed in the console.}
}
\description{
Location Estimation
}
\details{
This function uses the calibration table to evaluate the likelihood that the light has been recorded at any location of the grid.
}
