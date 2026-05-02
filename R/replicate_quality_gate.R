# replicate_quality_gate.R
# Post-processing quality gate for mcescher peak tables.
#
# Within each sequence processed by mcescher, all runs are replicates of the
# same sample.  For each compound the within-sequence distribution of
# area_cons is inspected with two complementary criteria:
#
#   1. Modified Z-score (Iglewicz-Hoaglin 1993, criterion |M_i| > 3.5):
#      flags individual runs whose area deviates more than z_threshold
#      scaled MADs from the within-sequence median.
#
#   2. Sequence CV% guard: if the within-sequence coefficient of variation
#      exceeds max_seq_cv_pct, ALL runs for that compound in the sequence
#      are flagged -- the sequence is too noisy regardless of which run is
#      the apparent outlier.  Typical causes: partially resolved co-elution
#      below the LOQ, baseline instability, column bleed episode.
#
# Sequences with fewer than min_n replicates for a compound pass through
# unchanged (quality_ok = TRUE, diagnostic columns = NA).

#' Flag unreliable peak areas using within-sequence replicate statistics.
#'
#' For each compound in a mcescher peak table, the within-sequence distribution
#' of \code{area_cons} is tested with a modified Z-score (Iglewicz-Hoaglin) and
#' a sequence-level CV% guard.  Rows that fail either criterion receive
#' \code{quality_ok = FALSE} and should be excluded from downstream use
#' (CNN training, APPAC calibration, publication tables).
#'
#' @param peaks_df     Data frame produced by \code{peaks(pipeline_summary(...))},
#'                     or a combined table from multiple sequences.  Must contain
#'                     columns \code{name} and \code{area_cons}.
#' @param min_n        Minimum number of replicates required before any criterion
#'                     is applied.  Compounds with fewer replicates pass through
#'                     with \code{quality_ok = TRUE} and \code{NA} diagnostics.
#'                     Default 5.
#' @param z_threshold  Iglewicz-Hoaglin modified Z-score threshold.  The standard
#'                     recommendation is 3.5.
#' @param max_seq_cv_pct  Sequence-level CV% threshold.  If the within-sequence
#'                     CV% of \code{area_cons} for a compound exceeds this value,
#'                     all replicates for that compound are flagged.  Set to
#'                     \code{Inf} to disable.  Default 3.0.
#'
#' @return \code{peaks_df} with four additional columns:
#'   \describe{
#'     \item{seq_n}{Integer. Number of replicates of this compound in the sequence.}
#'     \item{seq_cv_pct}{Numeric. Within-sequence CV\% of \code{area_cons}.}
#'     \item{modified_z}{Numeric. Iglewicz-Hoaglin modified Z-score.}
#'     \item{quality_ok}{Logical. \code{TRUE} if the row passes both criteria.}
#'   }
#'
#' @details
#' The modified Z-score is defined as \eqn{M_i = 0.6745 (x_i - \tilde{x}) /
#' \mathrm{MAD}}, where \eqn{\tilde{x}} is the within-group median and MAD uses
#' the consistency factor 1.4826 (normal-distribution equivalent).  Rows with
#' \eqn{|M_i| > } \code{z_threshold} are flagged.
#'
#' When \code{MAD = 0} (all replicates identical) only the CV\% criterion applies.
#'
#' The function is designed for use with one sequence at a time (as produced by
#' a single \code{compare_empower.R} or \code{run_ude.R} call).  It can also be
#' applied to a combined multi-sequence table if a \code{sequence_id} grouping
#' variable is added first and the function is called per group via
#' \code{split} / \code{lapply}.
#'
#' @references
#' Iglewicz, B. and Hoaglin, D. C. (1993).
#' \emph{How to Detect and Handle Outliers.}
#' ASQ Quality Press, Milwaukee, WI.
#'
#' @examples
#' \dontrun{
#' peaks_qc <- replicate_quality_gate(peaks_raw)
#' peaks_ok  <- peaks_qc[peaks_qc$quality_ok, ]
#' }
#'
#' @export
replicate_quality_gate <- function(
  peaks_df,
  min_n          = 5L,
  z_threshold    = 3.5,
  max_seq_cv_pct = 3.0
) {
  peaks_df$seq_n      <- NA_integer_
  peaks_df$seq_cv_pct <- NA_real_
  peaks_df$modified_z <- NA_real_
  peaks_df$quality_ok <- TRUE

  for (nm in unique(peaks_df$name)) {
    idx <- which(peaks_df$name == nm)
    a   <- peaks_df$area_cons[idx]
    n   <- length(a)
    peaks_df$seq_n[idx] <- n

    if (n < min_n) next

    med  <- stats::median(a, na.rm = TRUE)
    mad_ <- stats::mad(a, center = med, constant = 1.4826, na.rm = TRUE)
    mn   <- mean(a, na.rm = TRUE)
    cv   <- if (!is.na(mn) && mn > 0) 100 * stats::sd(a, na.rm = TRUE) / mn
            else NA_real_
    peaks_df$seq_cv_pct[idx] <- cv

    # Criterion 2: high sequence CV — flag entire compound for this sequence
    if (!is.na(cv) && is.finite(max_seq_cv_pct) && cv > max_seq_cv_pct) {
      peaks_df$quality_ok[idx] <- FALSE
      next
    }

    # Criterion 1: per-run modified Z-score
    if (!is.na(mad_) && mad_ > 0) {
      mz <- 0.6745 * (a - med) / mad_
      peaks_df$modified_z[idx] <- mz
      peaks_df$quality_ok[idx] <- abs(mz) <= z_threshold
    }
  }

  peaks_df
}
