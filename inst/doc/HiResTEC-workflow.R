## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_raw_data------------------------------------------------------------
raw <- HiResTEC::raw
sam <- data.frame(
  "gr" = rep(rep(c("A","B"), each=3), 3), 
  "tp" = rep(0:2, each=6), 
  "rep" = rep(1:3,6)
)
sam[,"ids"] <- apply(sam, 1, paste0, collapse="_")
head(sam)

## ----show_raw_bpc, fig.dim=c(8,6), out.width="98%"----------------------------
length(raw)
class(raw[[1]])
raw[[1]]@mzrange
bpc <- lapply(raw, HiResTEC::getMultipleBPC)
HiResTEC::plotBPC(bpc, type = "bpc", mfrow = c(3,6), ids = sam[,"ids"])

## ----show_spectrum------------------------------------------------------------
rt <- 1026.5
s <- HiResTEC::DeconvoluteSpectrum(dat = raw, rt = rt)
InterpretMSSpectrum::PlotSpec(x = s)

## ----show_mid, fig.dim=c(8,6), out.width="98%"--------------------------------
mz <- 556.263
mz_dev <- 0.4
bpc <- lapply(raw, function(x) { HiResTEC::getMultipleBPC(x = x, mz = mz + c(0:6)*1.003355, mz_dev = mz_dev) })
flt <- sam[,"tp"] %in% c(0, 2) & sam[,"gr"] == "A"
HiResTEC::plotBPC(bpc[flt], mfrow = c(2,3), ids = sam[flt,"ids"])

## ----get_peak_list------------------------------------------------------------
int <- sapply(bpc, function(x) { x[attr(x, "maxBPC"),] })
colnames(int) <- sam[,"ids"]
xg <- data.frame("mz" = as.numeric(rownames(int)), "rt" = rep(rt,7), int, row.names = NULL)
xg[,1:8]

## ----EvaluatePairsFromXCMSSet, results="hide"---------------------------------
preCL <- HiResTEC::EvaluatePairsFromXCMSSet(xg = xg, tp = sam$tp, gr = sam$gr, dmz = 0.04)

## ----EvaluatePairsFromXCMSSet_result------------------------------------------
head(preCL[order(preCL[,"P"]),3:7])

## ----EvaluateCandidateListAgainstRawData, results="hide"----------------------
finCL <- HiResTEC::EvaluateCandidateListAgainstRawData(
  x = preCL, tp = sam$tp, gr = sam$gr, dat = raw, dmz = 0.04, rolp = "all"
)

## ----QCplot1, fig.dim=c(9,5), out.width="98%", echo=FALSE---------------------
HiResTEC:::CandidateBoxplot(res = finCL[[1]])

## ----QCplot2, fig.dim=c(8,6), out.width="98%", echo=FALSE---------------------
flt <- sam[,"rep"] == 1
bpc <- finCL[[1]][["bpc"]]
HiResTEC::plotBPC(bpc[flt], mfrow = c(3,2), ids = sam[flt,"ids"])

## ----QCplot3, fig.dim=c(8,6), out.width="98%", echo=FALSE---------------------
x <- finCL[[1]]
mz1 <- x[["mz1"]]
mz2 <- x[["mz2"]]
tmp.txt <- data.frame("pos" = c(mz1, mz2), "val" = c("mz1", "mz2"))
ylim <- c(0, max(c(x[["s"]][abs(x[["s"]][, "mz"] - (mz1 + 4)) < 10, "int"], x[["s2"]][abs(x[["s2"]][, "mz"] - (mz1 + 4)) < 10, "int"]), na.rm = T))
graphics::layout(matrix(c(1, 2), ncol = 1))
InterpretMSSpectrum::PlotSpec(x[["s"]], rellab = x[["s"]][which.min(abs(x[["s"]][, 1] - mz1)), 1], masslab = 0, xlim = mz1 + c(-6, 15), txt = tmp.txt, ylim = ylim)
InterpretMSSpectrum::PlotSpec(x[["s2"]], rellab = x[["s2"]][which.min(abs(x[["s2"]][, 1] - mz2)), 1], masslab = 0, xlim = mz + c(-6, 15), txt = tmp.txt, ylim = ylim)

## ----CorMID-------------------------------------------------------------------
fml <- "C22H46N5O4Si4"
attr(fml, "nbio") <- 5

# re-extract the BPCs including [M+]
mz <- finCL[[1]]$mz1
rt <- finCL[[1]]$rt
bpc <- lapply(raw, function(x) { 
  HiResTEC::getMultipleBPC(x = x, mz = mz + c(-1:6)*1.003355, mz_dev = 0.04, rt = rt) 
})

# BPCs show about 2.5% [M+] intensity, define r accordingly
r <- setNames(c(0.975, 0.025), nm = c("M+H", "M+"))
int <- sapply(bpc, function(x) { x[attr(x, "maxBPC"),] })
rownames(int) <- paste0("M",-1:6)
colnames(int) <- sam$ids
mid <- apply(int, 2, function(x) { 
  CorMID::CorMID(int = x, fml = fml, r = r) 
})

# show mean group corrected MID
sapply(split(as.data.frame(t(mid)), interaction(sam[,"gr"], sam[,"tp"])), function(x) {
  round(apply(x,2,mean),1)
})

