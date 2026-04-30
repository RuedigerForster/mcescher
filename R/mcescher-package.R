#' @description
#' \if{html}{\figure{mcescher.png}{options: style='float: right;' alt='mcescher logo' width='120'}}
#'
#' Automated chromatographic peak integration with uncertainty quantification
#' for target compound analysis.
#' @import Matrix
#' @import methods
#' @importFrom grDevices adjustcolor dev.off png
#' @importFrom graphics abline legend lines par polygon text
#' @importFrom stats approx coef convolve cor gaussian glm median nls
#'   nls.control predict quantile reshape sd
#' @importFrom utils capture.output head
"_PACKAGE"

#' Grob mix chromatographic dataset
#'
#' A numeric matrix of 18 GC-FID chromatograms of a Grob test mixture
#' at varying concentration levels, used as the package example dataset.
#'
#' @format A numeric matrix with 53998 rows (time points) and 18 columns
#'   (chromatograms). Column names are acquisition timestamps. The
#'   \code{sample_name} attribute holds human-readable concentration labels
#'   (e.g. \code{"OGE 75ug/mL_1"}).
#'
#' @source UDE/OGE collaboration dataset (Grob mix, 2022).
#' @usage data(Grob_mix)
"Grob_mix"
