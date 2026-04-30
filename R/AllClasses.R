#' ChromResult S4 class — unified container for chromatography pipeline output.
#'
#' @slot meta          Named list: run metadata and method settings.
#' @slot peaks         data.frame: one row per detected peak (chrom, RT, name,
#'   height, width, area_*, ensemble statistics, uncertainty columns).
#' @slot lod_loq       data.frame: per-chromatogram LOD/LOQ estimates.
#' @slot alignment     data.frame or NULL: per-segment ATSA alignment shifts.
#' @slot baseline_rmse data.frame or NULL: per-chromatogram baseline RMSE.
#' @slot report        character or NULL: full text pipeline report captured
#'   from \code{pipeline_summary()}.
#' @slot traces        list or NULL: per-chromatogram signal/baseline trace
#'   data.frames for use by \code{plot()}.
#' @exportClass ChromResult
setClass("ChromResult",
  slots = c(
    meta          = "list",
    peaks         = "data.frame",
    lod_loq       = "data.frame",
    alignment     = "ANY",
    baseline_rmse = "ANY",
    report        = "ANY",
    traces        = "ANY"
  ),
  prototype = list(
    meta          = list(),
    peaks         = data.frame(),
    lod_loq       = data.frame(),
    alignment     = NULL,
    baseline_rmse = NULL,
    report        = NULL,
    traces        = NULL
  )
)
