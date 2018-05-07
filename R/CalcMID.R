#'@title CalcMID.
#'
#'@description
#'\code{CalcMID} will compute a MID (Mass Isotopomer Distribution) based on measured ion intensities in GC-APCI-MS.
#'
#'@details
#'Let's assume we measured the ion intensities of all 3 isotopes of an individual compound containing 2 carbons and observe a vector of {978,22,0}.
#'We may calculate the enrichment (E) out of this data, i.e. the relative proportion of 13C vs total carbon which will amount to about 1.1% (the natural 13C abundance) under standard conditions.
#'The equvivalent MID vector would be {1,0,0}, indicating that the non-labeled isotopologue (where non-labeled means non-labled above the natural 1.1%) is the only component observed.
#'During a labelling experiment we may change the measurement values in different ways (either labelling only one carbon or both), which potentially can translate into similar values for E being larger 1.1%.
#'The MIDs will provide additional information about the isotopolouge fraction which gave rise to the observed E's (cf. examples).
#'
#'@param int Vector of measured ion intensities of a fragment.
#'@param fml Chemical formula of fragment.
#'@param ratio If NULL M+H/M+ ratio will be determined from the data if necessary. Can be specified explicitly here.
#'@param nmz Attached as attr to fml for \link{CalcTheoreticalMDV}.
#'@param nbio Attached as attr to fml for \link{CalcTheoreticalMDV}.
#'
#'@return
#'Percent representation of each isotopologue measured (=MIDs).
#'
#'@export
#'
#'@examples
#'#tbd

CalcMID <- function(int=NULL, fml="", ratio=NULL, nmz=NULL, nbio=NULL) {

  length.out <- 3
  inperc <- TRUE

  # Helper functions
  FitMID <- function(md=NULL, td=NULL, r=0, length.out=3) {
    # md normalized (i.e. sum=1) raw intensity vector
    # td theoretical intensity distribution (using function 'CalcTheoreticalMDV')
    # r ratio of M+/M+H
    # Helper Function
    poss_local <- function(mid=NULL, d=NULL) {
      lst <- lapply(mid, function(x) {seq(x-d, x+d, length.out = length.out)})
      lst <- lapply(lst, function(x) {x[x>=0 & x<=1]})
      tmp <- expand.grid(lst)
      flt <- apply(tmp,1,sum)==1
      return(tmp[flt,])
    }
    # set starting conditions
    mid_start <- rep(0.5,length(md)-ifelse(r!=-1,1,0))
    names(mid_start) <- paste0("M",0:(length(mid_start)-1))
    # approximate best MID itteratively
    for (d in c(0.5,0.25,0.125,0.0625,0.03125,0.01,0.005,0.001)) {
      brute_force_local <- poss_local(mid=mid_start, d=d)
      test <- apply(brute_force_local,1,function(x){
        # [JL] shouldnt the next line compare against nrow??
        if (ncol(td)!=length(x)) {
          warning("Something went wrong. Check your 'nmz' parameter.")
          stop()
          #browser()
        }
        midc <- apply(td*x,2,sum)
        if (r!=-1) midc <- (r*c(midc,0)+(1-r)*c(0,midc))
        return(sqrt(sum((midc/sum(midc)-md/sum(md))^2)))
      })
      mid_start <- brute_force_local[which.min(test),]
    }
    return(unlist(mid_start))
  }

  # set non-finite values to zero
  if (any(!is.finite(int))) int[!is.finite(int)] <- 0

  # compute theoretic distribution:
  attr(fml, "nbio") <- nbio
  attr(fml, "nmz") <- nmz
  theo_dist <- CalcTheoreticalMDV(fml=fml)

  # compute starting parameters for optimization
  if (sum(int)==0) return(rep(NA, ncol(theo_dist)))
  md <- int/sum(int)
  if (length(int)>ncol(theo_dist)) {
    if (is.null(ratio)) ratio <- int[1]/(int[1]+int[2]-int[2]*theo_dist[1,2]/theo_dist[1,1])
  } else {
    ratio <- -1
  }
  if (is.finite(ratio)) {
    out <- FitMID(md=md, td=theo_dist, r=ratio, length.out = length.out)
  } else {
    out <- unlist(sapply(theo_dist[,1],function(x){return(NA)}))
  }

  # return relative representation of each isotopologue measured (= corrected MIDs in per cent)
  if (inperc==TRUE) out <- round(100*out,2)
  return(out)
}
