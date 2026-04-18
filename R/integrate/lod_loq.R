# lod_loq.R
# Per-chromatogram LOD and LOQ estimation.
#
# Approach (ICH Q2(R1) signal-to-noise):
#   σ_noise  = robust noise estimate from peak-free baseline regions
#   LOD      = 3.3 * σ_noise   (S/N = 3.3)
#   LOQ      = 10  * σ_noise   (S/N = 10)
#
# Values are in detector signal units.  If a calibration slope is supplied
# (signal units / concentration), LOD_conc and LOQ_conc are also returned.
#
# Functions:
#   estimate_noise(signal, RT, peaks, method = "mad")
#   lod_loq(signal, RT, peaks,
#           method_lod = NULL, method_loq = NULL,
#           cal_slope  = NULL,
#           noise_method = "mad")
#   lod_loq_batch(mat, RT, peaks_list, ...)
#
# Noise estimation methods:
#   "mad"  – median absolute deviation of the baseline-region signal
#             (robust, preferred for asymmetric distributions)
#   "sd"   – standard deviation of baseline-region signal
#   "rms"  – root mean square (appropriate when mean ≈ 0 after baseline removal)


# =============================================================================
# Internal: identify baseline (peak-free) regions
# =============================================================================

# Returns a logical vector: TRUE = this sample is far enough from any peak
# to be considered baseline.
.peak_free_mask <- function(n, peaks, guard_factor = 2.0) {
  mask <- rep(TRUE, n)
  if (is.null(peaks) || nrow(peaks) == 0L) return(mask)
  for (i in seq_len(nrow(peaks))) {
    lo <- max(1L, round(peaks$lo_idx[i] - guard_factor * (peaks$idx[i] - peaks$lo_idx[i])))
    hi <- min(n,  round(peaks$hi_idx[i] + guard_factor * (peaks$hi_idx[i] - peaks$idx[i])))
    mask[lo:hi] <- FALSE
  }
  mask
}


# =============================================================================
# Noise estimation
# =============================================================================

#' Estimate signal noise from peak-free regions.
#'
#' @param signal       Baseline-corrected signal (residual after baseline removal).
#' @param RT           Retention time vector (minutes).
#' @param peaks        Data frame from find_peaks() used to mask peak regions.
#'                     If NULL, uses the entire signal.
#' @param method       "mad" | "sd" | "rms"
#' @param guard_factor Multiplier on half-width to extend the exclusion zone
#'                     around each peak.
#'
#' @return Named list: noise_sd, n_baseline_pts, method, baseline_RT_ranges
estimate_noise <- function(signal, RT, peaks = NULL,
                           method       = c("mad", "sd", "rms"),
                           guard_factor = 2.0) {
  method <- match.arg(method)
  mask   <- .peak_free_mask(length(signal), peaks, guard_factor)
  bl_sig <- signal[mask]

  if (length(bl_sig) < 5L) {
    warning("Fewer than 5 baseline points available; noise estimate unreliable.")
    return(list(noise_sd = NA_real_, n_baseline_pts = sum(mask),
                method = method, baseline_RT_ranges = NULL))
  }

  noise_sd <- switch(method,
    mad = stats::mad(bl_sig, constant = 1.4826, na.rm = TRUE),
    sd  = sd(bl_sig, na.rm = TRUE),
    rms = sqrt(mean(bl_sig^2, na.rm = TRUE))
  )

  # identify contiguous baseline ranges for reporting
  runs <- rle(mask)
  ends   <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  bl_ranges <- lapply(seq_along(runs$values), function(k) {
    if (runs$values[k])
      c(RT[starts[k]], RT[ends[k]])
    else NULL
  })
  bl_ranges <- Filter(Negate(is.null), bl_ranges)

  list(noise_sd        = noise_sd,
       n_baseline_pts  = sum(mask),
       method          = method,
       baseline_RT_ranges = bl_ranges)
}


# =============================================================================
# LOD / LOQ per chromatogram
# =============================================================================

#' Compute LOD and LOQ for a single chromatogram.
#'
#' @param signal       Signal vector (ideally baseline-corrected).
#' @param RT           Retention time vector (minutes).
#' @param peaks        Detected peaks data frame (find_peaks output).
#' @param method_lod   Method LOD in the same signal units (for comparison).
#'                     NULL = no comparison.
#' @param method_loq   Method LOQ in signal units. NULL = no comparison.
#' @param cal_slope    Calibration slope (signal / concentration unit).
#'                     If given, LOD_conc and LOQ_conc are returned.
#' @param noise_method "mad" | "sd" | "rms"
#'
#' @return Data frame with one row containing:
#'   noise_sd, LOD, LOQ,
#'   LOD_conc (if cal_slope given), LOQ_conc (if cal_slope given),
#'   LOD_ok (TRUE if LOD ≤ method_lod), LOQ_ok,
#'   n_peaks_above_LOD, n_peaks_above_LOQ
lod_loq <- function(signal, RT, peaks = NULL,
                    method_lod   = NULL,
                    method_loq   = NULL,
                    cal_slope    = NULL,
                    noise_method = c("mad", "sd", "rms")) {
  noise_method <- match.arg(noise_method)
  ne  <- estimate_noise(signal, RT, peaks, method = noise_method)
  sig <- ne$noise_sd

  LOD <- 3.3 * sig
  LOQ <- 10.0 * sig

  row <- data.frame(
    noise_sd          = sig,
    noise_n_pts       = ne$n_baseline_pts,
    noise_method      = noise_method,
    LOD               = LOD,
    LOQ               = LOQ,
    LOD_conc          = NA_real_,
    LOQ_conc          = NA_real_,
    method_LOD        = if (!is.null(method_lod)) method_lod else NA_real_,
    method_LOQ        = if (!is.null(method_loq)) method_loq else NA_real_,
    LOD_ok            = if (!is.null(method_lod)) LOD <= method_lod else NA,
    LOQ_ok            = if (!is.null(method_loq)) LOQ <= method_loq else NA,
    n_peaks_above_LOD = NA_integer_,
    n_peaks_above_LOQ = NA_integer_,
    stringsAsFactors  = FALSE
  )

  if (!is.null(cal_slope) && !is.na(cal_slope) && cal_slope > 0) {
    row$LOD_conc <- LOD / cal_slope
    row$LOQ_conc <- LOQ / cal_slope
  }

  if (!is.null(peaks) && nrow(peaks) > 0L) {
    row$n_peaks_above_LOD <- sum(peaks$height >= LOD, na.rm = TRUE)
    row$n_peaks_above_LOQ <- sum(peaks$height >= LOQ, na.rm = TRUE)
  }

  row
}


# =============================================================================
# Batch: one row per chromatogram column
# =============================================================================

#' Apply lod_loq() to every column of a matrix.
#'
#' @param mat        n_samples × n_chromatograms matrix.
#' @param RT         Retention time vector (minutes).
#' @param peaks_list List of peaks data frames, one per column (or NULL for all).
#' @param ...        Additional arguments passed to lod_loq().
#'
#' @return Data frame: one row per chromatogram, with a column "chrom_idx".
lod_loq_batch <- function(mat, RT, peaks_list = NULL, ...) {
  nc  <- ncol(mat)
  rows <- vector("list", nc)
  for (j in seq_len(nc)) {
    pk  <- if (!is.null(peaks_list)) peaks_list[[j]] else NULL
    row <- lod_loq(mat[, j], RT, peaks = pk, ...)
    row$chrom_idx <- j
    rows[[j]] <- row
  }
  do.call(rbind, rows)
}
