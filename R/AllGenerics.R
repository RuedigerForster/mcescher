#' Extract the peaks table from a ChromResult.
#'
#' @param object A \code{ChromResult} object.
#' @return A data.frame of detected peaks.
#' @export
setGeneric("peaks", function(object) standardGeneric("peaks"))

#' Extract run metadata from a ChromResult.
#'
#' @param object A \code{ChromResult} object.
#' @return A named list of run metadata and method settings.
#' @export
setGeneric("meta", function(object) standardGeneric("meta"))

#' Extract per-chromatogram signal/baseline traces from a ChromResult.
#'
#' @param object A \code{ChromResult} object.
#' @return A list of trace data.frames (one per chromatogram).
#' @export
setGeneric("traces", function(object) standardGeneric("traces"))
