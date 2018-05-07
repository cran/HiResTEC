#'@title plotBPC.
#'
#'@description
#'\code{plotBPC} will plot for each item of a list of result-ojects from \link{getMultipleBPC} the BPC traces and the spectrum at the scan where the summed intensity of all ions is max.
#'
#'@details
#'not yet
#'
#'@param bpc A bpc object (list of intensity matrixes, rt x mz, including several attributes as attached by \link{getMultipleBPC}).
#'@param mfrow Specify mfrow explicitely (is optimized internally if NULL to cover n=length(bpc)).
#'@param skip_plots Allows to block certain subplots in the mfrow matrix to bettern align replicates.
#'@param ylim Can be specified specifically, will be adjusted to overall min/max otherwise.
#'@param col Specific color vector for masses may be provided.
#'@param ids Specific plot ids may be explicitely provided.
#'
#'@return
#'A plot to the graphics device and NULL as invisible.
#'
#'@examples
#'#load example raw data
#'data(res_list)
#'plotBPC(bpc = res_list[[1]][["bpc"]][c(1:2,13:14)])
#'
#'@export
#'
#'@importFrom graphics lines
#'@importFrom graphics axis
#'@importFrom graphics abline
#'@importFrom graphics text
#'@importFrom grDevices grey
#'@importFrom graphics box

plotBPC <- function(bpc=NULL, mfrow=NULL, skip_plots=NULL, ylim=NULL, col=NULL, ids=NULL) {

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # number of valid samples
  n <- ifelse(is.list(bpc), length(bpc), 1)

  # remove empty samples
  if (any(sapply(bpc, is.null))) {
    idx <- which(sapply(bpc, is.null))
    if (length(idx)==n) stop(print("All samples appear empty."))
    bpc <- bpc[-idx]
    ids <- ids[-idx]
    n <- n-length(idx)
    for (i in idx) {
      skip_plots <- c(skip_plots, i+sum(skip_plots<i))
    }
  }

  # check for differing number of ion traces in samples
  num_ion_traces <- unique(sapply(bpc,ncol))
  if (!(length(num_ion_traces)==1)) {
    warning("Different numbers of ion traces exported in bpc. Possibly one file is empty")
    stop()
  }

  # set color vector
  if (is.null(col)) {
    if (num_ion_traces>=2) {
      col <- 1:num_ion_traces
    } else {
      col <- grDevices::grey(0.9)
    }
  } else {
    if (length(col)!=num_ion_traces) col <- rep(col, length.out=num_ion_traces)
  }

  # set ylim
  if (is.null(ylim)) {
    ylim <- c(
      min(sapply(bpc, function(x) { ifelse(any(is.finite(unlist(x))), min(x, na.rm=T), 100) })),
      max(sapply(bpc, function(x) { ifelse(any(is.finite(unlist(x))), max(x, na.rm=T), 100) }))
    )
  }

  # switch between BPC only and spectra aside plot
  if (num_ion_traces>=2) {
    if (is.null(mfrow)) mfrow <- grDevices::n2mfrow(n)*c(1,2) else mfrow <- mfrow*c(1,2)

    par(mfrow=mfrow)

    cor_val <- 0
    for (i in 1:(n+length(skip_plots))) {
      if (i %in% skip_plots) {
        cor_val <- cor_val+1
        for (i in 1:2) plot(1,1,axes=F,ann=F,type="n")
      } else {
        j <- i-cor_val
        if (is.null(bpc[[j]])) {
          for (i in 1:2) plot(1,1,axes=F,ann=F,type="n")
        } else {
          int <- bpc[[j]]
          rt <- attr(bpc[[j]],"rt")
          # plot isotope BPCs of sample
          par(mar=c(2,2.4,3,0))
          plot(y=int[,1], x=rt, main=ifelse(is.null(ids),j,ids[j]), type="l", col=col[1], ylim=ylim, xlab="RT", ylab="", log="y", panel.first=abline(v=rt[attr(bpc[[j]],"maxBPC")], lty=1, lwd=3, col=grDevices::grey(0.9)))
          points(int[,1]~rt, pch=21, bg=col[1])
          for (k in 2:ncol(int)) {
            lines(int[,k]~rt,col=col[k])
            points(int[,k]~rt, pch=rep(c(21,22,24,25),3)[k], bg=col[k])
          }
          # plot spectrum of sample
          spec <- bpc[[j]][attr(bpc[[j]], "maxBPC"),]
          spec <- spec/ifelse(any(is.finite(unlist(spec))), max(spec, na.rm=T), 1)
          mz <- attr(bpc[[j]], "mz")
          par(mar=c(2,0.5,3,0.1))
          plot(y=spec, x=mz, type="h", lwd=3, ylim=c(0,1), xlim=c(floor(min(mz)),ceiling(max(mz))), ylab="", xlab="", main="", col=col, panel.first=abline(h=seq(0,1,0.1), lty=2, lwd=1, col=grDevices::grey(0.9)), axes=FALSE)
          axis(1); box()
          # annotate with mass defect
          md <- attr(bpc[[j]], "mass_defect")
          for (k in 1:length(md)) text(y=spec[k], x=mz[k], labels=md[k], adj=c(ifelse(spec[k]>0.8,1,0), 1.25), srt=90, col=grDevices::grey(0.4), cex=1)
        }
      }
    }

  } else {
    # simplified version (only BPC, no spectrum)
    if (is.null(mfrow)) mfrow <- grDevices::n2mfrow(n) else mfrow <- mfrow
    par(mfrow=mfrow)
    cor_val <- 0
    for (i in 1:(n+length(skip_plots))) {
      if (i %in% skip_plots) {
        cor_val <- cor_val+1
        plot(1,1,axes=F,ann=F,type="n")
      } else {
        j <- i-cor_val
        int <- bpc[[j]]
        rt <- attr(bpc[[j]],"rt")
        # plot BPC of sample j
        par(mar=c(2,2.4,3,0))
        plot(y=int[,1], x=rt, main=ifelse(is.null(ids),j,ids[j]), type="l", ylim=ylim, xlab="RT", ylab="", log="y", panel.first=abline(v=rt[attr(bpc[[j]],"maxBPC")], lty=1, lwd=3, col=col))
      }
    }
  }
  invisible(NULL)
}
