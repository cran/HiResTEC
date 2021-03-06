% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMID.R
\name{plotMID}
\alias{plotMID}
\title{plotMID.}
\usage{
plotMID(mid = NULL, gr = NULL, name = "unknown", contr = NULL,
  stackedbars = FALSE, subplot_ylim = 100, ...)
}
\arguments{
\item{mid}{Matrix of measured ion intensities corrected using CalcMID.}

\item{gr}{Groups, a factor.}

\item{name}{Name of metabolite.}

\item{contr}{Contrasts. Not yet clear if useful.}

\item{stackedbars}{Alternative plotting layout using stacked bar plot.}

\item{subplot_ylim}{Calculate ylim individually per subplot if 0, show full range in all subplots if 100 and limit to the minimal specified number otherwise.}

\item{...}{Further arguments to 'boxplot'.}
}
\value{
NULL.
}
\description{
\code{plotMID} will plot a Mass Isotopomer Distribution (MID) as calculated by CalcMID.
}
\details{
Not yet.
}
\examples{
mid <- matrix(c(seq(0,0.3,0.1), seq(1,0.7,-0.1)), byrow=TRUE, nrow=2)
gr <- gl(2,2,labels=letters[1:2])
plotMID(mid=mid, gr=gr, name="Metabolite X")
plotMID(mid=mid, gr=gr, stackedbars=TRUE, las=1, col=2:3, xlab="MID")

}
