# peak_integrate.R
library(parallel)
# Peak area integration via four methods:
#   "PD"    – Perpendicular Drop  (O'Haver / measurepeaks.m)
#   "TS"    – Tangent Skim / valley–valley (O'Haver)
#   "gauss" – Gaussian nonlinear fit
#   "EGH"   – Exponentially Modified Gaussian Hybrid (Lan et al. 2001,
#              J. Chromatogr. A 915:1-13; MATLAB port of ExponentialGaussian.m)
#
# Main function:
#   integrate_peaks(signal, RT, peaks, baseline = NULL,
#                   methods = c("PD","TS","gauss","EGH"))
#
# Args:
#   signal   : numeric vector, raw or baseline-corrected chromatogram
#   RT       : retention time vector (minutes), same length as signal
#   peaks    : data frame from find_peaks() with columns
#              idx, RT, height, width, lo_idx, hi_idx
#   baseline : optional numeric vector (same length as signal).
#              If supplied, signal – baseline is used for PD/TS area,
#              but the raw signal is used for fitting boundaries.
#   methods  : character vector, subset of c("PD","TS","gauss","EGH")
#
# Returns: data frame, one row per peak, columns:
#   peak_no, RT, height, width,
#   area_PD, area_TS, area_gauss, area_EGH,
#   fit_rmse_gauss, fit_rmse_EGH,
#   left_valley_RT, right_valley_RT


# =============================================================================
# Internal helpers
# =============================================================================

# Find left/right valley indices around peak at position pk_idx.
# Searches within ±search_width samples from the peak.
.valleys <- function(y, pk_idx, width_samples, search_factor = 1.3) {
  n  <- length(y)
  hw <- max(3L, round(width_samples * search_factor))
  lo <- max(1L, pk_idx - hw)
  hi <- min(n,  pk_idx + hw)

  left_val  <- lo + which.min(y[lo:pk_idx]) - 1L
  right_val <- pk_idx + which.min(y[pk_idx:hi]) - 1L
  c(left = left_val, right = right_val)
}

# Numerical integration (trapezoidal) of y over RT[lo:hi].
.trapz <- function(RT, y, lo, hi) {
  if (lo >= hi) return(0)
  idx <- lo:hi
  sum(diff(RT[idx]) * (y[idx[-length(idx)]] + y[idx[-1]]) / 2)
}

# Gaussian curve: H * exp(-(t - mu)^2 / (2 * sigma^2))
.gauss_curve <- function(t, H, mu, sigma) {
  H * exp(-(t - mu)^2 / (2 * sigma^2))
}

# EGH curve (Lan et al. 2001 eq.1):
# f(t) = H * exp(-(t-tR)^2 / (2*sigma^2 + tau*(t-tR)))  where denominator > 0, else 0
.egh_curve <- function(t, H, tR, sigma, tau) {
  denom <- 2 * sigma^2 + tau * (t - tR)
  out   <- ifelse(denom > 0, H * exp(-(t - tR)^2 / denom), 0)
  out
}

# EGH area (Lan et al. 2001, eq. A4-A5 via proportionality factor e0)
.egh_area <- function(H, sigma, tau) {
  t_angle <- atan(abs(tau) / sigma)
  # 6th-order polynomial for proportionality constant (ExponentialGaussian.m)
  factors <- c(4, -6.293724, 9.232834, -11.34291, 9.123978, -4.173753, 0.827797)
  e0 <- sum(factors * t_angle^(0:6))
  H * (sigma * sqrt(pi / 8) + abs(tau)) * e0
}


# =============================================================================
# Per-peak integration functions
# =============================================================================

.integrate_one_PD <- function(signal_bc, RT, pk_idx, width_samp) {
  v  <- .valleys(signal_bc, pk_idx, width_samp)
  lo <- v["left"]; hi <- v["right"]
  area <- .trapz(RT, signal_bc, lo, hi)
  list(area = area, left_RT = RT[lo], right_RT = RT[hi])
}

.integrate_one_TS <- function(signal_bc, RT, pk_idx, width_samp) {
  v   <- .valleys(signal_bc, pk_idx, width_samp)
  lo  <- v["left"]; hi <- v["right"]
  # subtract linear baseline between the two valley points
  bl  <- approx(c(RT[lo], RT[hi]),
                c(signal_bc[lo], signal_bc[hi]),
                xout = RT[lo:hi])$y
  seg <- signal_bc[lo:hi] - bl
  seg[seg < 0] <- 0
  area <- .trapz(RT, c(rep(0, lo - 1L), seg, rep(0, length(signal_bc) - hi)),
                 lo, hi)
  list(area = area, left_RT = RT[lo], right_RT = RT[hi])
}

.integrate_one_gauss <- function(signal_bc, RT, pk_idx, width_samp) {
  v   <- .valleys(signal_bc, pk_idx, width_samp)
  lo  <- v["left"]; hi <- v["right"]
  seg_x <- RT[lo:hi]; seg_y <- signal_bc[lo:hi]

  H0     <- signal_bc[pk_idx]
  mu0    <- RT[pk_idx]
  sigma0 <- max(1e-6, (RT[hi] - RT[lo]) / 4)

  fit <- tryCatch(
    nls(seg_y ~ .gauss_curve(seg_x, H, mu, sigma),
        start    = list(H = H0, mu = mu0, sigma = sigma0),
        lower    = list(H = 0,  mu = RT[lo], sigma = 1e-6),
        algorithm = "port",
        control  = nls.control(maxiter = 50, warnOnly = TRUE)),
    error = function(e) NULL)

  if (is.null(fit)) {
    # fall back: analytical Gaussian from width at half-maximum
    sigma_hw <- max(1e-6, (RT[hi] - RT[lo]) / (2 * sqrt(2 * log(2))))
    area     <- H0 * sigma_hw * sqrt(2 * pi)
    return(list(area = area, rmse = NA_real_, left_RT = RT[lo], right_RT = RT[hi]))
  }
  cf    <- coef(fit)
  area  <- abs(cf["H"]) * abs(cf["sigma"]) * sqrt(2 * pi)
  yhat  <- predict(fit)
  rmse  <- sqrt(mean((seg_y - yhat)^2))
  list(area = area, rmse = rmse, left_RT = RT[lo], right_RT = RT[hi])
}

.integrate_one_EGH <- function(signal_bc, RT, pk_idx, width_samp) {
  v   <- .valleys(signal_bc, pk_idx, width_samp)
  lo  <- v["left"]; hi <- v["right"]
  seg_x <- RT[lo:hi]; seg_y <- signal_bc[lo:hi]

  H0     <- signal_bc[pk_idx]
  tR0    <- RT[pk_idx]
  sigma0 <- max(1e-6, (RT[hi] - RT[lo]) / 4)
  tau0   <- 0  # start symmetric; let optimizer find asymmetry

  fit <- tryCatch(
    nls(seg_y ~ .egh_curve(seg_x, H, tR, sigma, tau),
        start     = list(H = H0, tR = tR0, sigma = sigma0, tau = tau0),
        lower     = list(H = 0,  tR = RT[lo], sigma = 1e-6, tau = -10 * sigma0),
        upper     = list(H = 10 * H0, tR = RT[hi], sigma = 10 * sigma0, tau = 10 * sigma0),
        algorithm = "port",
        control   = nls.control(maxiter = 50, warnOnly = TRUE)),
    error = function(e) NULL)

  if (is.null(fit)) {
    # fall back to Gaussian
    res <- .integrate_one_gauss(signal_bc, RT, pk_idx, width_samp)
    res$rmse <- NA_real_
    return(res)
  }
  cf   <- coef(fit)
  area <- tryCatch(.egh_area(cf["H"], cf["sigma"], cf["tau"]),
                   error = function(e) NA_real_)
  if (is.na(area) || area <= 0) {
    area <- abs(cf["H"]) * abs(cf["sigma"]) * sqrt(2 * pi)  # Gaussian fallback
  }
  yhat <- predict(fit)
  rmse <- sqrt(mean((seg_y - yhat)^2))
  list(area = area, rmse = rmse, left_RT = RT[lo], right_RT = RT[hi])
}


# =============================================================================
# Public interface
# =============================================================================

#' Integrate chromatographic peaks by one or more methods.
#'
#' @param signal   Raw or smoothed signal vector.
#' @param RT       Retention time vector (minutes).
#' @param peaks    Data frame from find_peaks() (cols: idx, RT, height, width, lo_idx, hi_idx).
#' @param baseline Optional baseline vector. If given, signal–baseline is used for
#'                 PD/TS; fitting methods operate on the baseline-subtracted signal.
#' @param methods  Character vector: any of "PD", "TS", "gauss", "EGH".
#'
#' @return Data frame, one row per peak.
integrate_peaks <- function(signal, RT, peaks,
                            baseline = NULL,
                            methods  = c("PD", "TS", "gauss", "EGH")) {
  methods <- match.arg(methods, c("PD", "TS", "gauss", "EGH"), several.ok = TRUE)
  dt      <- RT[2] - RT[1]

  signal_bc <- if (!is.null(baseline)) signal - baseline else signal

  n_peaks <- nrow(peaks)
  if (n_peaks == 0L) return(data.frame())

  out <- data.frame(
    peak_no        = seq_len(n_peaks),
    RT             = peaks$RT,
    height         = peaks$height,
    width          = peaks$width,
    area_PD        = NA_real_,
    area_TS        = NA_real_,
    area_gauss     = NA_real_,
    area_EGH       = NA_real_,
    fit_rmse_gauss = NA_real_,
    fit_rmse_EGH   = NA_real_,
    left_valley_RT = NA_real_,
    right_valley_RT= NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_peaks)) {
    pk  <- peaks$idx[i]
    if (!is.null(peaks$lo_idx) && !is.na(peaks$lo_idx[i])) {
      w_samp <- max(3L, peaks$hi_idx[i] - peaks$lo_idx[i])
    } else {
      w_samp <- max(3L, round(peaks$width[i] / dt))
    }

    if ("PD" %in% methods) {
      r <- tryCatch(.integrate_one_PD(signal_bc, RT, pk, w_samp),
                    error = function(e) list(area = NA_real_, left_RT = NA_real_, right_RT = NA_real_))
      out$area_PD[i]         <- r$area
      out$left_valley_RT[i]  <- r$left_RT
      out$right_valley_RT[i] <- r$right_RT
    }
    if ("TS" %in% methods) {
      r <- tryCatch(.integrate_one_TS(signal_bc, RT, pk, w_samp),
                    error = function(e) list(area = NA_real_, left_RT = NA_real_, right_RT = NA_real_))
      out$area_TS[i] <- r$area
      if (is.na(out$left_valley_RT[i])) {
        out$left_valley_RT[i]  <- r$left_RT
        out$right_valley_RT[i] <- r$right_RT
      }
    }
    if ("gauss" %in% methods) {
      r <- tryCatch(.integrate_one_gauss(signal_bc, RT, pk, w_samp),
                    error = function(e) list(area = NA_real_, rmse = NA_real_,
                                             left_RT = NA_real_, right_RT = NA_real_))
      out$area_gauss[i]     <- r$area
      out$fit_rmse_gauss[i] <- r$rmse
    }
    if ("EGH" %in% methods) {
      r <- tryCatch(.integrate_one_EGH(signal_bc, RT, pk, w_samp),
                    error = function(e) list(area = NA_real_, rmse = NA_real_,
                                             left_RT = NA_real_, right_RT = NA_real_))
      out$area_EGH[i]     <- r$area
      out$fit_rmse_EGH[i] <- r$rmse
    }
  }

  # drop columns for methods not requested
  keep_cols <- c("peak_no", "RT", "height", "width",
                 if ("PD"    %in% methods) "area_PD",
                 if ("TS"    %in% methods) "area_TS",
                 if ("gauss" %in% methods) c("area_gauss", "fit_rmse_gauss"),
                 if ("EGH"   %in% methods) c("area_EGH",   "fit_rmse_EGH"),
                 "left_valley_RT", "right_valley_RT")
  out[, keep_cols]
}
