% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenerateQCPlots.R
\name{GenerateQCPlots}
\alias{GenerateQCPlots}
\title{GenerateQCPlots.}
\usage{
GenerateQCPlots(
  res_list = NULL,
  pdf_file = NULL,
  mfrow = NULL,
  skip_plots = NULL
)
}
\arguments{
\item{res_list}{A list of result objects (each testing an individual mz pair).}

\item{pdf_file}{Either APCI or ESI. Choice will modify some internal parameters and checks performed.}

\item{mfrow}{If NULL automatically determined, otherwise useful to specify a layout.}

\item{skip_plots}{NULL or numeric vector in which case plots with numbers in skip_plots will be empty.}
}
\value{
NULL.
}
\description{
\code{GenerateQCPlots} will produce QC plots for a list containing test results objects.
}
\details{
For individual candidates screen output is reasonable, otherwise a target PDF file should be specified.
}
\examples{
# load evaluation result of example data
data(res_list)
# generate Figures on screen (use PDF output for mass evaluation)
GenerateQCPlots(res_list[1])

}
