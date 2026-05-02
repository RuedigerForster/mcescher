# loq_filter.R
# Per-row limit-of-quantification filter for mcescher peak tables.
#
# Filters out individual peak areas whose mcescher integration uncertainty
# (area_cv_pct) exceeds a threshold.  This is complementary to
# replicate_quality_gate(), which detects cross-run outliers:
#
#   replicate_quality_gate()  -- flags runs that are outliers within a sequence
#   loq_filter()              -- flags runs where the signal itself is too weak
#                                to yield a reliable area estimate
#
# Typical pipeline:
#   peaks_ok <- peaks_raw |>
#     replicate_quality_gate() |>
#     loq_filter()             |>
#     (\(df) df[df$quality_ok, ])()

#' Filter peak areas below the limit of quantification.
#'
#' Marks rows in a mcescher peak table where the per-peak integration
#' uncertainty \code{area_cv_pct} exceeds \code{max_area_cv_pct}.  These rows
#' are unreliable as training labels for machine-learning models and should be
#' excluded from downstream use.
#'
#' The filter operates row-wise on the theoretical mcescher uncertainty
#' (baseline ensemble spread + detector noise + SASS denoising artifact,
#' combined in quadrature).  It is independent of replicate context and can
#' therefore be applied to single-run tables as well as sequence batches.
#'
#' If \code{quality_ok} already exists in \code{peaks_df} (e.g. from
#' \code{\link{replicate_quality_gate}}), the LOQ verdict is ANDed with it so
#' that a single column reflects all quality criteria.
#'
#' @param peaks_df       Data frame with at least columns \code{area_cv_pct}
#'                       and \code{name}.
#' @param max_area_cv_pct  Maximum acceptable integration CV\% (k=1 relative
#'                       standard uncertainty, in percent).  Rows with
#'                       \code{area_cv_pct > max_area_cv_pct} receive
#'                       \code{quality_ok = FALSE}.  Default 5.0.
#'                       \cr\cr
#'                       Rationale: empirical analysis of the OGE-UDE Grob mix
#'                       dataset shows a sharp instability zone for co-eluting
#'                       peak pairs (e.g. Decane/Thiophene) in the 5-20 mV·s
#'                       area range, where mcescher reports CV\% > 100\%.
#'                       A 5\% threshold eliminates this zone while retaining
#'                       > 98\% of well-resolved, single-peak entries.
#'
#' @return \code{peaks_df} with two additional or updated columns:
#'   \describe{
#'     \item{loq_ok}{\code{logical}. \code{TRUE} if \code{area_cv_pct <=
#'       max_area_cv_pct} and \code{area_cv_pct} is not \code{NA}.}
#'     \item{quality_ok}{\code{logical}. \code{TRUE} if the row passes all
#'       applied quality criteria.  Created if absent; ANDed with existing
#'       values if already present.}
#'   }
#'
#' @seealso \code{\link{replicate_quality_gate}}
#'
#' @examples
#' \dontrun{
#' peaks_qc <- peaks_raw |>
#'   replicate_quality_gate() |>
#'   loq_filter()
#' peaks_ok <- peaks_qc[peaks_qc$quality_ok, ]
#' }
#'
#' @export
loq_filter <- function(peaks_df, max_area_cv_pct = 5.0) {
  cv <- peaks_df$area_cv_pct
  loq_ok <- !is.na(cv) & cv <= max_area_cv_pct

  peaks_df$loq_ok <- loq_ok

  if ("quality_ok" %in% names(peaks_df)) {
    peaks_df$quality_ok <- peaks_df$quality_ok & loq_ok
  } else {
    peaks_df$quality_ok <- loq_ok
  }

  peaks_df
}
