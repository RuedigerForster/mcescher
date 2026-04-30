# peak_detect.R
# Peak detection for chromatographic data.
# Ported from findpeaksplot.m (T.C. O'Haver, 1995-2016).
#
# Strategy: pseudo-Gaussian smoothing of the first derivative,
# then detect downward zero-crossings that exceed slope/amplitude thresholds.
# Peak width is estimated at half-maximum.  Area is approximated as the
# Gaussian-peak formula: area = 1.0646 * height * width.


# ---- internal smoothing helpers ---------------------------------------------

.sliding_avg <- function(y, w) {
  # Single-pass boxcar (sliding average) of width w.
  w <- max(1L, round(w))
  n <- length(y)
  out <- numeric(n)
  cs  <- cumsum(c(0, y))
  hw  <- floor(w / 2)
  for (i in seq_len(n)) {
    lo <- max(1L, i - hw)
    hi <- min(n,  i + hw)
    out[i] <- (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
  }
  out
}

.pseudo_gaussian <- function(y, w) {
  # Three passes of sliding average  Gaussian convolution.
  if (w < 2) return(y)
  .sliding_avg(.sliding_avg(.sliding_avg(y, w), w), w)
}

.central_diff <- function(y) {
  # First derivative via two-point central difference.
  n <- length(y)
  d <- numeric(n)
  d[1]   <- y[2] - y[1]
  d[n]   <- y[n] - y[n - 1L]
  if (n > 2)
    d[2:(n-1)] <- (y[3:n] - y[1:(n-2)]) / 2
  d
}


# ---- log-parabola Gaussian fit (O'Haver gaussfit) ---------------------------

.gaussfit <- function(x, y) {
  # Returns peak position by fitting a parabola to log(y).
  # Falls back to the x at max(y) on failure.
  maxy <- max(y)
  y[y < maxy / 100] <- maxy / 100
  logy <- log(y)
  coef <- tryCatch(
    suppressWarnings(.polyfit_3(x, logy)),
    error = function(e) NULL)
  if (is.null(coef) || any(is.na(coef)) || coef[1] >= 0)
    return(x[which.max(y)])
  pos <- -coef[2] / (2 * coef[1])
  # If the vertex lies outside the fitting window the parabola is extrapolating
  # (happens when the smoothed-derivative zero-crossing is offset from the raw
  # apex, leaving the apex at the window edge). Fall back to the raw maximum.
  if (pos < x[1L] || pos > x[length(x)]) return(x[which.max(y)])
  pos
}

.polyfit_3 <- function(x, y) {
  # Least-squares fit of degree-2 polynomial; returns coefficients [a, b, c].
  X <- cbind(x^2, x, 1)
  qr.solve(X, y)
}


# ---- public interface -------------------------------------------------------

.find_peaks <- function(y, RT,
                       smooth_width = 11,
                       slope_thresh = 0,
                       amp_thresh   = 0,
                       peak_group   = 5) {
  n  <- length(y)
  dt <- if (n > 1) RT[2] - RT[1] else 1

  # smooth derivative
  d <- .central_diff(y)
  if (smooth_width > 1) d <- .pseudo_gaussian(d, smooth_width)

  pg    <- max(3L, round(peak_group))
  half  <- round(pg / 2)
  start <- max(2L, half + 1L)
  stop  <- min(n - 1L, n - half)

  peaks <- vector("list", 200L)
  k     <- 0L

  for (j in start:stop) {
    # downward zero-crossing
    if (sign(d[j]) <= sign(d[j + 1L])) next
    if ((d[j] - d[j + 1L]) <= slope_thresh) next
    if (y[j] <= amp_thresh) next

    # sub-group for position refinement
    lo  <- max(1L, j - half)
    hi  <- min(n,  j + half)
    pos <- tryCatch(.gaussfit(RT[lo:hi], y[lo:hi]),
                    error = function(e) RT[j])
    h   <- y[j]

    # half-maximum width
    lo_w <- j; while (lo_w > 1L && y[lo_w] > h / 2) lo_w <- lo_w - 1L
    hi_w <- j; while (hi_w < n  && y[hi_w] > h / 2) hi_w <- hi_w + 1L
    w <- max(dt, RT[hi_w] - RT[lo_w])

    k <- k + 1L
    peaks[[k]] <- list(idx    = j,
                       RT     = pos,
                       height = h,
                       width  = w,
                       area   = GAUSS_AREA_FACTOR * h * w,
                       lo_idx = lo_w,
                       hi_idx = hi_w)
  }

  if (k == 0L)
    return(data.frame(idx=integer(), RT=numeric(), height=numeric(),
                      width=numeric(), area=numeric(),
                      lo_idx=integer(), hi_idx=integer()))

  result <- do.call(rbind, lapply(peaks[1:k], as.data.frame))
  result[order(result$RT), ]
}
