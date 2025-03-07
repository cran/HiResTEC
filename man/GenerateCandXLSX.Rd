% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenerateCandXLSX.R
\name{GenerateCandXLSX}
\alias{GenerateCandXLSX}
\title{Generate a table for the candidates obtained in \link{EvaluateCandidateListAgainstRawData}.}
\usage{
GenerateCandXLSX(res_list = NULL, xlsx_file = NULL, rejected = FALSE)
}
\arguments{
\item{res_list}{A list of result objects (each testing an individual mz pair).}

\item{xlsx_file}{File name.}

\item{rejected}{Logical. Prepare table of rejected candidates if TRUE.}
}
\value{
Candidate table as data.frame.
}
\description{
\code{GenerateCandXLSX} will produce a XLSX of a list containing
    test results objects.
}
\details{
Just a wrapper, to get the important information in a tabular layout.
}
\examples{
# load evaluation result of example data and
# generate table within R (use parameter xlsx_file to write to file)
x <- GenerateCandXLSX(HiResTEC::res_list)
str(x)
x[,1:5]

}
