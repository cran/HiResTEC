% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMultipleBPC.R
\name{getMultipleBPC}
\alias{getMultipleBPC}
\title{getMultipleBPC.}
\usage{
getMultipleBPC(
  x,
  mz = NULL,
  mz_dev = 0.005,
  rt = NULL,
  rt_dev = 2,
  zeroVal = NA,
  smooth = 0,
  returnEIC = FALSE
)
}
\arguments{
\item{x}{xcmsRaw object.}

\item{mz}{mass vector or NULL (default) to return the TIC.}

\item{mz_dev}{allowed deviations (can be a single numeric, a vector, a matrix with one row (lower bound, upper bound) or a matrix with \code{length(mz)} rows giving lower and upper bound for each mz).}

\item{rt}{target time point or NULL (default) to use full scan time.}

\item{rt_dev}{allowed window.}

\item{zeroVal}{Set values <=0 to NA or keep as is with NULL.}

\item{smooth}{Window size for moving average smoother, 0 = no smoothing.}

\item{returnEIC}{Return EIC instead of BPC?}
}
\value{
A matrix with scan wise (rows) intensities for all requested masses (columns)
    as either EIC or BPC.
}
\description{
\code{getMultipleBPC} will extract multiple BPCs from an xcmsRaw
    object for a vector of mz within the limits given by rt, rt_dev and mz_dev.
}
\details{
While there are other functions to extract BPC information from raw data,
    this one is particularly useful to get all traces belonging to a isotopologue
    group. It will attach several derived values to the results object,
    i.e. describing the observed mass shift (deviation from expected value) which
    is helpful in QC for non-targeted tracer analyses.
}
\examples{
# see \link{plotBPC} for an example
}
\references{
Uses C code modified from XCMS (see \code{citation("xcms")}).
}
