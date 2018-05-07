#'@title CalcTheoreticalMDV.
#'
#'@description
#'\code{CalcTheoreticalMDV} will compute the Mass Distribution Vectors of isotopologues as it is used for correction matrix in \code{\link{CalcMID}} computations.
#'
#'@details
#'\code{CalcTheoreticalMDV} basically is a convenience function using Rdisop to generate the isotopologue distribution at natural abundance of 13C for a given formula.
#'It will break this down into a matrix where the components of the MID constitute the rows and the expected relative ion intensities are within the columns.
#'The number of exported ion intensities and MID components can be limited if numeric values for "nmz" and/or "nbio" provided as attributes with the formula.
#'
#'@param fml The chemical formula of the compound.
#'
#'@return
#'A matrix with dimensions according to the attributes of \code{fml} or the number of carbons respectively.
#'
#'@export
#'
#'@importFrom Rdisop getMolecule
#'
#'@examples
#'# standard distribution matrix
#'fml <- "C5H6Si1"
#'CalcTheoreticalMDV(fml=fml)
#'attr(fml,"nmz") <- 4
#'CalcTheoreticalMDV(fml=fml)
#'attr(fml,"nbio") <- 2
#'CalcTheoreticalMDV(fml=fml)

CalcTheoreticalMDV <- function(fml=NULL) {
  # establish number of ions measured (estimate from formula if not provided)
  if (is.null(attr(fml, "nmz"))) {
    cce <- InterpretMSSpectrum::CountChemicalElements(x=fml, ele=c("C","Si"))
    #cce <- CountChemicalElements(x=fml, ele=c("C","Si"))
    maxisotopes <- cce["C"]-3*cce["Si"]
  } else {
    maxisotopes <- attr(fml, "nmz")
  }
  # establish number of biological carbons (estimate from formula if not provided)
  if (is.null(attr(fml, "nbio")) || attr(fml, "nbio")>maxisotopes) {
    maxmid <- maxisotopes
  } else {
    maxmid <- attr(fml, "nbio")
  }
  # get base distribution and convert into matrix including normalization
  bd <- Rdisop::getMolecule(fml, maxisotopes = maxisotopes+1)$isotopes[[1]][2,]
  td <- matrix(NA, ncol=length(bd), nrow=maxmid+1, dimnames=list(paste0("M",0:maxmid),paste0("M+",0:maxisotopes)))
  for (i in 1:nrow(td)) td[i,] <- c(rep(0,i-1), bd[1:(length(bd)-(i-1))])
  td <- t(apply(td, 1, function(x) {x/sum(x)}))
  # return matrix
  return(td)
}
