#' @title GenerateQCPlots.
#'
#' @description
#' \code{GenerateQCPlots} will produce QC plots for a list containing test results objects.
#'
#' @details
#' For individual candidates screen output is reasonable, otherwise a target PDF file should be specified.
#'
#' @param res_list A list of result objects (each testing an individual mz pair).
#' @param pdf_file Either APCI or ESI. Choice will modify some internal parameters and checks performed.
#' @param mfrow If NULL automatically determined, otherwise useful to specify a layout.
#' @param skip_plots NULL or numeric vector in which case plots with numbers in skip_plots will be empty.
#'
#' @return
#' NULL.
#'
#' @examples
#' # load evaluation result of example data
#' data(res_list)
#' # generate Figures on screen (use PDF output for mass evaluation)
#' GenerateQCPlots(res_list[1])
#'
#' @importFrom InterpretMSSpectrum PlotSpec
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom graphics layout
#'
#' @export
#'
GenerateQCPlots <- function(res_list = NULL, pdf_file = NULL, mfrow = NULL, skip_plots = NULL) {
  # open connection to pdf if file name is given
  if (!is.null(pdf_file)) {
    grDevices::pdf(file = pdf_file, width = 14)
    on.exit(grDevices::dev.off())
  } else {
    # ensure that par is restored
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
  }

  # convert back to list in case user selected single item
  if (!is.null(names(res_list)) && names(res_list)[1] == "mz1") res_list <- list(res_list)

  for (i in 1:length(res_list)) {
    res <- res_list[[i]]

    # boxplots with stats
    CandidateBoxplot(res = res)

    # start asking after first plot
    if (is.null(pdf_file)) par(ask = TRUE)

    # plot BPC traces
    # specify mfrow if not provided by user
    n <- length(res[["tp"]])
    if (is.null(mfrow)) {
      if (n > 12) {
        mfrow_local <- grDevices::n2mfrow(length(levels(res[["inter"]])))
        # take one sample from each gr*tp with intermediate enrichment
        tmp <- attr(res[["enr"]], "Enrichment")
        if (length(tmp) == n) {
          names(tmp) <- 1:n
        } else {
          browser()
        }
        mfrow_select <- sapply(split(tmp, res[["inter"]]), function(x) {
          if (any(is.finite(x))) x <- x[is.finite(x)] else x <- x[1]
          if (length(x) == 1) {
            return(as.numeric(names(x)))
          } else {
            return(as.numeric(names(sort(x)[ceiling(length(x) / 2)])))
          }
        })
      } else {
        mfrow_local <- c(4, 3)
        mfrow_select <- 1:n
      }
    } else {
      mfrow_local <- mfrow
      mfrow_select <- 1:n
    }
    plotBPC(res[["bpc"]][mfrow_select], mfrow = mfrow_local, ids = paste0(as.character(res[["inter"]])[mfrow_select], " [", mfrow_select, "]"), skip_plots = skip_plots)

    # plot spectra
    mz <- res[["mz1"]]
    mz2 <- res[["mz2"]]
    tmp.col <- rep(1, nrow(res[["s"]]))

    olmz <- as.numeric(strsplit(res[["OverlappingMasses"]], ", ")[[1]])
    for (m in olmz) if (any(abs(res[["s"]] - m) < 0.1)) tmp.col[which.min(abs(res[["s"]] - m))] <- 3
    for (m in c(mz, mz2)) if (any(abs(res[["s"]] - m) < 0.1)) tmp.col[which.min(abs(res[["s"]] - m))] <- 2
    if (nrow(res[["s"]]) >= 1) {
      tmp.txt <- data.frame("pos" = c(mz, mz2), "val" = c("mz1", "mz2"))
      graphics::layout(matrix(1, ncol = 1))
      InterpretMSSpectrum::PlotSpec(res[["s"]], masslab = 1, cutoff = 0, txt = tmp.txt, cols = tmp.col, ylim = c(0, max(res[["s"]][, 2])))
      ylim <- c(0, max(c(res[["s"]][abs(res[["s"]][, "mz"] - (mz + 4)) < 10, "int"], res[["s2"]][abs(res[["s2"]][, "mz"] - (mz + 4)) < 10, "int"]), na.rm = T))
      if (nrow(res[["s2"]]) >= 1) {
        graphics::layout(matrix(c(1, 2), ncol = 1))
        InterpretMSSpectrum::PlotSpec(res[["s"]], rellab = res[["s"]][which.min(abs(res[["s"]][, 1] - mz)), 1], masslab = 0, xlim = mz + c(-6, 15), txt = tmp.txt, ylim = ylim)
        InterpretMSSpectrum::PlotSpec(res[["s2"]], rellab = res[["s2"]][which.min(abs(res[["s2"]][, 1] - mz2)), 1], masslab = 0, xlim = mz + c(-6, 15), txt = tmp.txt, ylim = ylim)
      } else {
        InterpretMSSpectrum::PlotSpec(res[["s"]], rellab = TRUE, masslab = 0.1, xlim = mz + c(-1, 1) * 20, txt = tmp.txt, cols = tmp.col, ylim = ylim)
      }
    } else {
      graphics::layout(matrix(1, ncol = 1))
      for (i in 1:2) plot(1, 1, axes = F, ann = F, type = "n")
    }
  }

  invisible(NULL)
}
