% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvaluateSpectrum.R
\name{EvaluateSpectrum}
\alias{EvaluateSpectrum}
\title{EvaluateSpectrum.}
\usage{
EvaluateSpectrum(res = NULL, dat = NULL, drt = 2, dmz = 0.005,
  ionization = c("APCI", "ESI")[1], smooth = 0)
}
\arguments{
\item{res}{res-object.}

\item{dat}{List of raw data files.}

\item{drt}{Time deviation. Passed to \link{DeconvoluteSpectrum}.}

\item{dmz}{Mass deviation. Passed to \link{DeconvoluteSpectrum}.}

\item{ionization}{Ion source. Passed to \link{DeconvoluteSpectrum}.}

\item{smooth}{Smoothing parameter. Passed to \link{DeconvoluteSpectrum}.}
}
\value{
Extended res-object.
}
\description{
\code{EvaluateSpectrum} will extend a res-object (results) by information regarding mass spectra extracted from raw data files.
}
\details{
If still present, rows containing NA in mz2 column will be filtered additionally.
}
\keyword{internal}
