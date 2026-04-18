# atsa.R
# Automatic Time-Shift Alignment (ATSA) for chromatographic data.
#
# Reference:
#   Zheng et al. (2017) Automatic time-shift alignment method for
#   chromatographic data analysis. Sci. Rep. 7:256.
#   doi:10.1038/s41598-017-00390-7
#
# The algorithm has three stages:
#   1. Baseline correction (LMV-RSA) and peak detection.
#   2. Preliminary alignment: segment-based TPC-optimised shift + warping.
#   3. Precise alignment: peak-to-peak sub-segment alignment + warping.
#
# Usage:
#   source("./atsa/atsa.R")   # also sources lmv.R and peak_detect.R
#   result <- atsa(mat, RT)
#   aligned_mat <- result$aligned

source("./lmv/lmv.R")
source("./atsa/peak_detect.R")


# =============================================================================
# 1. Internal helpers
# =============================================================================

# Linear warp: stretch / compress a vector to n_out samples.
.warp <- function(y, n_out) {
  if (length(y) == n_out) return(y)
  if (length(y) < 2)       return(rep(y[1], n_out))
  approx(seq_along(y), y, xout = seq(1, length(y), length.out = n_out))$y
}

# Fill NA values in a vector using the nearest non-NA neighbour.
.fill_na <- function(x) {
  bad  <- which(is.na(x))
  if (!length(bad)) return(x)
  good <- which(!is.na(x))
  if (!length(good)) return(x)
  x[bad] <- x[good[pmin(findInterval(bad, good), length(good))]]
  x
}


# =============================================================================
# 2. TPC (Total Peak Correlation) criterion  — eq. (2-3) in the paper
# =============================================================================

# Compute TPC between a reference segment and a test segment.
#   ref_sub   : reference signal vector (already trimmed to segment)
#   test_sub  : test signal vector, already warped to same length as ref_sub
#   seg_peaks : data frame of peaks within ref_sub (idx re-based to 1)
#   dt        : time step (minutes per sample)
.tpc <- function(ref_sub, test_sub, seg_peaks, dt) {
  N <- nrow(seg_peaks)
  if (N == 0L) return(0)

  total_w <- 0; total_wc <- 0; I_matched <- 0L
  n_ref <- length(ref_sub)

  for (i in seq_len(N)) {
    pk <- seg_peaks$idx[i]
    hw <- max(1L, round(seg_peaks$width[i] / dt / 2))
    lo <- max(1L, pk - hw); hi <- min(n_ref, pk + hw)
    r_sub <- ref_sub[lo:hi]; t_sub <- test_sub[lo:hi]

    # Check whether a corresponding test peak is present (sd > 0)
    if (length(r_sub) < 2 || sd(r_sub) < 1e-12) {
      # reference peak with no counterpart in test → c_i = -1  (eq. note)
      c_i <- -1
    } else if (sd(t_sub) < 1e-12) {
      c_i <- -1
    } else {
      c_i <- cor(r_sub, t_sub)
      I_matched <- I_matched + 1L
    }
    w_i      <- seg_peaks$area[i] / max(1, seg_peaks$width[i] / dt)
    total_w  <- total_w  + w_i
    total_wc <- total_wc + w_i * c_i
  }

  if (total_w == 0) return(0)
  (total_wc / total_w) * (I_matched / N)
}


# =============================================================================
# 3. Segment initialisation (Stage 2a)
# =============================================================================

# Build reference segments from detected peaks.
#   ref_peaks    : output of find_peaks() for the reference chromatogram
#   RT           : full retention time vector
#   segment_size : initial segment size (minutes)
#   min_peaks    : segments with fewer peaks are merged into a neighbour
#
# Returns a data frame: begin, end (sample indices), begin_RT, end_RT
.build_segments <- function(ref_peaks, RT, segment_size = 3, min_peaks = 3) {
  n_total <- length(RT)

  if (nrow(ref_peaks) == 0L)
    return(data.frame(begin=1L, end=n_total,
                      begin_RT=RT[1], end_RT=RT[n_total]))

  pks <- ref_peaks[order(ref_peaks$RT), ]

  # --- assign each peak to a segment based on distance to segment-first peak
  seg_id      <- integer(nrow(pks))
  seg_id[1]   <- 1L
  seg_first_RT <- pks$RT[1]
  cur_seg     <- 1L

  for (i in 2:nrow(pks)) {
    if ((pks$RT[i] - seg_first_RT) >= segment_size) {
      cur_seg      <- cur_seg + 1L
      seg_first_RT <- pks$RT[i]
    }
    seg_id[i] <- cur_seg
  }
  n_seg <- cur_seg

  # --- merge segments with < min_peaks into the smaller neighbour
  repeat {
    counts  <- tabulate(seg_id, nbins = max(seg_id))
    small   <- which(counts > 0L & counts < min_peaks)
    if (!length(small)) break
    s <- small[1]
    if (s == 1L) {
      seg_id[seg_id == s] <- 2L
    } else {
      # merge with whichever neighbour has fewer peaks
      left  <- if (s > 1L)            counts[s - 1L] else Inf
      right <- if (s < length(counts)) counts[s + 1L] else Inf
      target <- if (left <= right) s - 1L else s + 1L
      seg_id[seg_id == s] <- target
    }
    # re-number contiguously
    unique_ids <- sort(unique(seg_id))
    new_ids    <- seq_along(unique_ids)
    for (k in seq_along(unique_ids))
      seg_id[seg_id == unique_ids[k]] <- new_ids[k]
  }
  n_seg <- max(seg_id)

  # --- compute segment boundaries from first/last peak positions
  rt_to_idx <- function(t) which.min(abs(RT - t))

  start_RT <- numeric(n_seg); end_RT <- numeric(n_seg)
  for (s in seq_len(n_seg)) {
    sp <- pks[seg_id == s, ]
    start_RT[s] <- min(sp$RT)
    end_RT[s]   <- max(sp$RT)
  }

  # --- modify boundaries to midpoints between adjacent segments
  for (s in seq_len(n_seg - 1L)) {
    mid         <- (end_RT[s] + start_RT[s + 1L]) / 2
    end_RT[s]   <- mid
    start_RT[s + 1L] <- mid
  }
  start_RT[1]    <- RT[1]
  end_RT[n_seg]  <- RT[n_total]

  data.frame(
    begin    = vapply(start_RT, rt_to_idx, 1L),
    end      = vapply(end_RT,   rt_to_idx, 1L),
    begin_RT = start_RT,
    end_RT   = end_RT
  )
}


# =============================================================================
# 4. Shift search per segment (Stage 2b)
# =============================================================================

.find_segment_shifts <- function(ref_sig, test_sig, ref_peaks,
                                  segments, shift_pts, dt) {
  n_seg  <- nrow(segments)
  shifts <- integer(n_seg)
  n_test <- length(test_sig)

  for (s in seq_len(n_seg)) {
    r_lo    <- segments$begin[s]; r_hi <- segments$end[s]
    ref_sub <- ref_sig[r_lo:r_hi]
    seg_len <- length(ref_sub)

    # reference peaks within segment, re-indexed to 1-based within segment
    sp <- ref_peaks[ref_peaks$idx >= r_lo & ref_peaks$idx <= r_hi, ]
    if (nrow(sp) > 0) sp$idx <- sp$idx - r_lo + 1L

    best_tpc <- -Inf; best_sh <- 0L

    for (sh in seq(-shift_pts, shift_pts)) {
      # test window is reference window displaced by -sh
      t_lo <- max(1L, r_lo - sh); t_hi <- min(n_test, r_hi - sh)
      test_raw <- test_sig[t_lo:t_hi]
      test_sub <- .warp(test_raw, seg_len)
      sc <- .tpc(ref_sub, test_sub, sp, dt)
      if (sc > best_tpc) { best_tpc <- sc; best_sh <- sh }
    }
    shifts[s] <- best_sh
  }
  shifts
}


# =============================================================================
# 5. Outlier shift correction (Stage 2c, eq. 4-5)
# =============================================================================

.correct_outlier_shifts <- function(shifts, ref_sig, test_sig,
                                     segments, shift_pts, dt,
                                     thresh = 2.5) {
  if (length(shifts) < 3L) return(shifts)

  med_s <- median(shifts)
  sigma <- 1.483 * median(abs(shifts - med_s))
  if (sigma < 1e-10) return(shifts)

  d_k      <- abs(shifts - med_s) / sigma
  outliers <- which(d_k > thresh)
  n_test   <- length(test_sig)

  for (s in outliers) {
    r_lo    <- segments$begin[s]; r_hi <- segments$end[s]
    ref_sub <- ref_sig[r_lo:r_hi]
    seg_len <- length(ref_sub)

    # re-align within ±2.5σ using plain Pearson correlation
    lo_sh <- round(med_s - 2.5 * sigma)
    hi_sh <- round(med_s + 2.5 * sigma)
    best_cor <- -Inf; best_sh <- round(med_s)

    for (sh in lo_sh:hi_sh) {
      t_lo <- max(1L, r_lo - sh); t_hi <- min(n_test, r_hi - sh)
      test_sub <- .warp(test_sig[t_lo:t_hi], seg_len)
      if (sd(test_sub) < 1e-12) next
      cc <- cor(ref_sub, test_sub)
      if (!is.na(cc) && cc > best_cor) { best_cor <- cc; best_sh <- sh }
    }
    shifts[s] <- best_sh
  }
  shifts
}


# =============================================================================
# 6. Warp a full chromatogram according to per-segment shifts (Stage 2d)
# =============================================================================

.warp_chromatogram <- function(test_sig, segments, shifts, n_out) {
  out    <- rep(NA_real_, n_out)
  n_test <- length(test_sig)

  for (s in seq_len(nrow(segments))) {
    r_lo    <- segments$begin[s]; r_hi <- segments$end[s]
    sh      <- shifts[s]
    t_lo    <- max(1L, r_lo - sh); t_hi <- min(n_test, r_hi - sh)
    seg_len <- r_hi - r_lo + 1L
    out[r_lo:r_hi] <- .warp(test_sig[t_lo:t_hi], seg_len)
  }
  .fill_na(out)
}


# =============================================================================
# 7. Precise alignment — Stage 3
# =============================================================================

# Sub-segment boundaries: midpoints between adjacent reference peak positions.
.subseg_bounds <- function(ref_peaks, n) {
  if (nrow(ref_peaks) == 0L) return(c(1L, n))
  pks    <- ref_peaks[order(ref_peaks$RT), ]
  bounds <- c(1L)
  for (i in seq_len(nrow(pks) - 1L)) {
    mid    <- round((pks$idx[i] + pks$idx[i + 1L]) / 2)
    bounds <- c(bounds, mid)
  }
  unique(sort(c(bounds, n)))
}

.precise_align <- function(ref_sig, test_sig, ref_peaks, test_peaks,
                            RT, dt, n) {
  if (nrow(ref_peaks) == 0L || nrow(test_peaks) == 0L) return(test_sig)
  ref_peaks  <- ref_peaks[order(ref_peaks$RT), ]
  test_peaks <- test_peaks[order(test_peaks$RT), ]
  bounds     <- .subseg_bounds(ref_peaks, n)
  n_sub      <- length(bounds) - 1L
  out        <- rep(NA_real_, n)
  n_test     <- length(test_sig)

  for (s in seq_len(n_sub)) {
    r_lo <- bounds[s]; r_hi <- bounds[s + 1L]
    seg_len <- r_hi - r_lo + 1L

    # reference peak(s) in this sub-segment
    sp <- ref_peaks[ref_peaks$idx >= r_lo & ref_peaks$idx <= r_hi, ]
    if (nrow(sp) == 0L) {
      out[r_lo:r_hi] <- .warp(test_sig[r_lo:min(n_test, r_hi)], seg_len)
      next
    }
    ref_pk_RT <- sp$RT[which.max(sp$height)]

    # nearest test peak
    nearest   <- which.min(abs(test_peaks$RT - ref_pk_RT))
    sh        <- round((ref_pk_RT - test_peaks$RT[nearest]) / dt)

    t_lo <- max(1L, r_lo - sh); t_hi <- min(n_test, r_hi - sh)
    out[r_lo:r_hi] <- .warp(test_sig[t_lo:t_hi], seg_len)
  }
  .fill_na(out)
}


# =============================================================================
# 8. Main entry point
# =============================================================================

#' Automatic Time-Shift Alignment for chromatographic matrices.
#'
#' @param mat            Numeric matrix, rows = time points, columns = samples.
#' @param RT             Retention time vector in minutes (length = nrow(mat)).
#' @param segment_size   Initial segment size for preliminary alignment (min).
#' @param shift_max      Maximum time-shift search range (min).
#' @param smooth_width   Smoothing window (samples) for peak detection.
#' @param amp_thresh     Minimum peak height for detection.
#' @param slope_thresh   Minimum derivative slope for peak detection.
#' @param baseline_correct  Apply LMV-RSA baseline correction before alignment.
#' @param lmv_span       Window parameter for LMV-RSA.
#' @param ref_peaks_RT   Optional numeric vector of known retention times (min)
#'                       to use as alignment anchors instead of auto-detection.
#'                       Typically derived from the method config peak RT_ref
#'                       values.  When supplied, find_peaks() is skipped for the
#'                       reference chromatogram; the signal maximum within
#'                       ±ref_peaks_width of each anchor is used as the peak
#'                       apex.  Test-column peaks (Stage 3) are still
#'                       auto-detected.
#' @param ref_peaks_width Half-width (min) of the search window around each
#'                       anchor when ref_peaks_RT is supplied.  Default 0.1 min.
#' @param verbose        Print progress messages.
#'
#' @return A list with:
#'   \item{aligned}{Aligned matrix (same dimensions as mat).}
#'   \item{baseline_corrected}{Baseline-corrected matrix used for alignment.}
#'   \item{shifts}{Integer matrix (n_segments x n_columns) of preliminary shifts.}
#'   \item{reference_idx}{Column index of the chosen reference chromatogram.}
#'   \item{ref_peaks}{Peaks data frame for the reference chromatogram.}
#'   \item{segments}{Segment table used for preliminary alignment.}
atsa <- function(mat, RT,
                 segment_size     = 3.0,
                 shift_max        = 0.5,
                 smooth_width     = 11L,
                 amp_thresh       = 0,
                 slope_thresh     = 0,
                 baseline_correct = TRUE,
                 lmv_span         = 50L,
                 ref_peaks_RT     = NULL,
                 ref_peaks_width  = 0.1,
                 verbose          = TRUE) {

  n  <- nrow(mat)
  nc <- ncol(mat)
  dt <- RT[2] - RT[1]
  shift_pts <- max(1L, round(shift_max / dt))

  # ---- Stage 1: baseline correction ----------------------------------------
  if (verbose) message(sprintf("[ATSA] Stage 1: baseline correction (%s)",
                                if (baseline_correct) "LMV-RSA" else "skipped"))
  mat_bc <- mat
  if (baseline_correct) {
    for (j in seq_len(nc)) {
      bl <- tryCatch(
        { b <- lmv_rsa(mat[, j], span = lmv_span); b[is.na(b)] <- 0; b },
        error = function(e) rep(0, n))
      mat_bc[, j] <- pmax(0, mat[, j] - bl)
    }
  }

  # ---- Stage 1: peak detection ---------------------------------------------
  if (verbose) message("[ATSA] Stage 1: peak detection")
  peaks_list <- lapply(seq_len(nc), function(j)
    find_peaks(mat_bc[, j], RT,
               smooth_width = smooth_width,
               amp_thresh   = amp_thresh,
               slope_thresh = slope_thresh))

  n_peaks <- vapply(peaks_list, nrow, 1L)
  if (verbose) message(sprintf("[ATSA]   peaks detected per column: min=%d, max=%d, median=%.0f",
                                min(n_peaks), max(n_peaks), median(n_peaks)))

  # ---- Reference selection: max mean pairwise Pearson correlation ----------
  if (verbose) message("[ATSA] Selecting reference chromatogram")
  cor_mat <- cor(mat_bc)
  ref_idx <- which.max(colMeans(cor_mat))
  ref_sig <- mat_bc[, ref_idx]
  if (verbose) message(sprintf("[ATSA]   reference = column %d", ref_idx))

  # ---- Reference peaks: from config anchors or auto-detection --------------
  if (!is.null(ref_peaks_RT)) {
    w_samp <- max(3L, round(ref_peaks_width / dt))
    ref_rows <- lapply(ref_peaks_RT, function(rt) {
      # centre on the signal maximum within ±ref_peaks_width of the anchor
      centre <- which.min(abs(RT - rt))
      lo     <- max(1L, centre - w_samp)
      hi     <- min(n,  centre + w_samp)
      pk     <- lo - 1L + which.max(ref_sig[lo:hi])
      h      <- ref_sig[pk]
      data.frame(idx    = pk,
                 RT     = RT[pk],
                 height = h,
                 width  = ref_peaks_width,
                 area   = h * ref_peaks_width,
                 lo_idx = max(1L, pk - w_samp),
                 hi_idx = min(n,  pk + w_samp))
    })
    ref_peaks <- do.call(rbind, ref_rows)
    ref_peaks <- ref_peaks[order(ref_peaks$RT), ]
    if (verbose)
      message(sprintf("[ATSA]   reference peaks from config (%d anchors)",
                      nrow(ref_peaks)))
  } else {
    ref_peaks <- peaks_list[[ref_idx]]
    if (verbose)
      message(sprintf("[ATSA]   reference peaks auto-detected (%d peaks)",
                      nrow(ref_peaks)))
  }

  # ---- Build reference segments --------------------------------------------
  segments <- .build_segments(ref_peaks, RT, segment_size)
  if (verbose) message(sprintf("[ATSA]   %d segments", nrow(segments)))

  # ---- Stage 2: preliminary alignment + warping ----------------------------
  if (verbose) message("[ATSA] Stage 2: preliminary alignment")
  aligned   <- mat_bc
  shift_log <- matrix(0L, nrow = nrow(segments), ncol = nc)

  for (j in seq_len(nc)) {
    if (j == ref_idx) next
    if (verbose) message(sprintf("[ATSA]   column %d / %d", j, nc))

    sh <- .find_segment_shifts(ref_sig, mat_bc[, j], ref_peaks,
                                segments, shift_pts, dt)
    sh <- .correct_outlier_shifts(sh, ref_sig, mat_bc[, j],
                                   segments, shift_pts, dt)
    shift_log[, j] <- sh
    aligned[, j]   <- .warp_chromatogram(mat_bc[, j], segments, sh, n)
  }

  # ---- Stage 3: precise alignment ------------------------------------------
  if (verbose) message("[ATSA] Stage 3: precise alignment")
  for (j in seq_len(nc)) {
    if (j == ref_idx) next
    new_peaks  <- find_peaks(aligned[, j], RT,
                              smooth_width = smooth_width,
                              amp_thresh   = amp_thresh,
                              slope_thresh = slope_thresh)
    aligned[, j] <- .precise_align(ref_sig, aligned[, j],
                                    ref_peaks, new_peaks, RT, dt, n)
  }

  if (verbose) message("[ATSA] Done.")

  list(aligned            = aligned,
       baseline_corrected = mat_bc,
       shifts             = shift_log,
       reference_idx      = ref_idx,
       ref_peaks          = ref_peaks,
       segments           = segments)
}
