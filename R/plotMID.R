#'@title plotMID.
#'
#'@description
#'\code{plotMID} will plot a Mass Isotopomer Distribution (MID) as calculated by CalcMID.
#'
#'@details
#'Not yet.
#'
#'@param mid Matrix of measured ion intensities corrected using CalcMID.
#'@param gr Groups, a factor.
#'@param name Name of metabolite.
#'@param contr Contrasts. Not yet clear if useful.
#'@param stackedbars Alternative plotting layout using stacked bar plot.
#'@param subplot_ylim Calculate ylim individually per subplot if 0, show full range in all subplots if 100 and limit to the minimal specified number otherwise.
#'@param ... Further arguments to 'boxplot'.
#'
#'@return
#'NULL.
#'
#'@importFrom graphics barplot
#'@importFrom grDevices grey
#'
#'@export
#'
#'@examples
#'mid <- matrix(c(seq(0,0.3,0.1), seq(1,0.7,-0.1)), byrow=TRUE, nrow=2)
#'gr <- gl(2,2,labels=letters[1:2])
#'plotMID(mid=mid, gr=gr, name="Metabolite X")
#'plotMID(mid=mid, gr=gr, stackedbars=TRUE, las=1, col=2:3, xlab="MID")
#'
plotMID <- function(mid=NULL, gr=NULL, name="unknown", contr=NULL, stackedbars=FALSE, subplot_ylim=100, ...) {
  if (is.null(gr)) gr <- gl(n = 1, k = ncol(mid))
  if (stackedbars) {
    # get group medians
    tmp <- t(apply(mid, 1, function(x) { sapply(split(x, gr), median) }))
    # readjust to sum=100
    tmp <- apply(tmp, 2, function(x) {100*x/sum(x)})
    graphics::barplot(tmp, ylab=name, ...)
  } else {
    tmp <- apply(mid, 1, function(x) { split(x, gr) })
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow=c(1,length(tmp)))
    par(mar=c(3,4,1,0)+0.5)
    ylim <- range(mid, na.rm=T)
    for (k in 1:length(tmp)) {
      if (subplot_ylim<100) {
        ylim <- range(tmp[[k]], na.rm=T)
        while(diff(ylim)<subplot_ylim) {
          ylim <- round(ylim+c(-1,1))
          if (ylim[1]<0) ylim[1] <- 0
          if (ylim[2]>100) ylim[2] <- 100
        }
      }
      graphics::boxplot(tmp[[k]], main="", ylab="", ylim=ylim, ...)
      graphics::mtext(text = paste0("M",k-1), side = 3, adj = 1)
      if (k==1) {
        graphics::title(ylab=name, cex.lab=2)
      }
      if (!is.null(contr) && all(contr %in% 1:length(tmp[[1]]))) {
        if (length(contr)==1) {
          for (l in (1:length(tmp[[1]]))[-contr]) {
            p <- try(stats::t.test(tmp[[k]][[l]], tmp[[k]][[contr]])$p.value)
            if (is.numeric(p)) {
              pos <- ifelse(mean(tmp[[k]][[l]])<50,3,1)
              graphics::text(x = l+0.5, y = mean(tmp[[k]][[l]]), labels = formatC(p, format = "e", digits = 1), pos = pos, col=ifelse(p<0.01, 2, grDevices::grey(0.9)))
            }
          }
        } else {
          for (l in contr) {
            p <- try(stats::t.test(tmp[[k]][[l]], tmp[[k]][[l+1]])$p.value)
            if (is.numeric(p)) {
              pos <- ifelse(mean(tmp[[k]][[l]])<50,3,1)
              graphics::text(x = l+0.5, y = mean(tmp[[k]][[l]]), labels = formatC(p, format = "e", digits = 1), pos = pos, col=ifelse(p<0.01, 2, grDevices::grey(0.9)))
            }
          }
        }
      }
    }
    par(mfrow=c(1,1))
  }
  invisible(NULL)
}
