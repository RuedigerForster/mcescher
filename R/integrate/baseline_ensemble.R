# baseline_ensemble.R
library(parallel)
# Estimate peak-area uncertainty due to baseline algorithm and parameter choice.
#
# Strategy: run a defined set of (algorithm, parameters) combinations on the
# same chromatogram, integrate each resulting baseline-corrected signal with a
# chosen method, then summarise the spread across combinations per peak.
#
# Algorithms wrapped here (all must be sourced before calling):
#   airPLS  – source("./airpls/airPLS.R")
#   arPLS   – source("./arpls/doArPLS.R")
#   LMV-RSA – source("./lmv/lmv.R")
#   SASS    – source("./sass/sass.R")
#
# Main function:
#   baseline_ensemble(signal, RT, peaks,
#                     integrate_method = "TS",
#                     configs          = default_ensemble_configs())
#
# Returns a list:
#   $areas        – matrix n_peaks × n_configs with individual areas
#   $summary      – data frame per peak: RT, height, area_mean, area_sd,
#                   area_cv_pct, area_min, area_max, area_range
#   $baselines    – matrix n_samples × n_configs with fitted baselines
#   $config_names – character vector labelling each configuration

source("./airpls/airPLS.R")
source("./arpls/doArPLS.R")
source("./lmv/lmv.R")
source("./flatfit/flatfit.R")
source("./integrate/peak_integrate.R")

# Compile and load the BEADS C++ baseline (compiled once per session).
# sourceCpp caches the shared library; subsequent loads are fast.
if (!exists("beads_baseline", mode = "function")) {
  Rcpp::sourceCpp("./beads/BEADS/beads_rcpp.cpp")
}


# =============================================================================
# Algorithm wrappers  (signal → baseline vector)
# =============================================================================

.bl_airpls <- function(signal, lambda = 1e7, order = 2, wep = 0.1, p = 0.05) {
  tryCatch({
    mat <- matrix(signal, nrow = 1)
    res <- airPLS(mat, lambda = lambda, order = order, wep = wep, p = p)
    as.numeric(res$Z[1, ])   # airPLS returns list(Xc, Z); Z[1,] is the baseline row
  }, error = function(e) rep(0, length(signal)))
}

.bl_arpls <- function(signal, lambda = 1e7, ratio = 1e-6) {
  tryCatch({
    res <- doArPLS(signal, lambda = lambda, ratio = ratio)
    as.numeric(res$z)   # doArPLS returns list(z, bslPts); z is the baseline
  }, error = function(e) rep(0, length(signal)))
}

.bl_lmv <- function(signal, span = 50) {
  tryCatch({
    bl <- lmv_rsa(signal, span = span)
    bl[is.na(bl)] <- 0
    bl
  }, error = function(e) rep(0, length(signal)))
}

.bl_flatfit <- function(signal, smoothness = 1e-3, p = 0.05) {
  tryCatch({
    flatfit_baseline(signal, smoothness = smoothness, p = p)
  }, error = function(e) rep(0, length(signal)))
}

# BEADS – Baseline Estimation and Denoising with Sparsity
# (Ning, Selesnick, Duval 2014; C++ port via Rcpp/RcppEigen)
# fc: filter cut-off frequency (cycles/sample); smaller = slower baseline.
# Default regularisation: lam0=0.4, lam1=4.0, lam2=3.2 (amp=0.8 scaling).
.bl_beads <- function(signal, d = 1, fc = 0.01,
                       r = 6.0, lam0 = 0.4, lam1 = 4.0, lam2 = 3.2) {
  tryCatch({
    as.numeric(beads_baseline(signal, d = d, fc = fc, r = r,
                               lam0 = lam0, lam1 = lam1, lam2 = lam2))
  }, error = function(e) rep(0, length(signal)))
}


# =============================================================================
# Default configuration set
# =============================================================================

#' Return the default list of (algorithm, parameter) configurations.
#' Each element is a named list: name, fun (function(signal) → baseline).
default_ensemble_configs <- function() {
  list(
    list(name = "airPLS_lam1e6",  fun = function(s) .bl_airpls(s, lambda = 1e6)),
    list(name = "airPLS_lam1e7",  fun = function(s) .bl_airpls(s, lambda = 1e7)),
    list(name = "airPLS_lam1e8",  fun = function(s) .bl_airpls(s, lambda = 1e8)),
    list(name = "arPLS_lam1e6",   fun = function(s) .bl_arpls(s,  lambda = 1e6)),
    list(name = "arPLS_lam1e7",   fun = function(s) .bl_arpls(s,  lambda = 1e7)),
    list(name = "arPLS_lam1e8",   fun = function(s) .bl_arpls(s,  lambda = 1e8)),
    list(name = "LMV_span30",     fun = function(s) .bl_lmv(s,    span = 30)),
    list(name = "LMV_span50",     fun = function(s) .bl_lmv(s,    span = 50)),
    list(name = "LMV_span100",    fun = function(s) .bl_lmv(s,    span = 100)),
    list(name = "flatfit_p003",    fun = function(s) .bl_flatfit(s, p = 0.03)),
    list(name = "flatfit_p005",    fun = function(s) .bl_flatfit(s, p = 0.05)),
    list(name = "flatfit_p010",    fun = function(s) .bl_flatfit(s, p = 0.10)),
    list(name = "BEADS_fc0005",    fun = function(s) .bl_beads(s, fc = 0.005)),
    list(name = "BEADS_fc001",     fun = function(s) .bl_beads(s, fc = 0.01)),
    list(name = "BEADS_fc005",     fun = function(s) .bl_beads(s, fc = 0.05))
  )
}


# =============================================================================
# Main function
# =============================================================================

#' Compute peak area uncertainty across a set of baseline configurations.
#'
#' @param signal           Raw signal vector.
#' @param RT               Retention time vector (minutes).
#' @param peaks            Data frame from find_peaks().
#' @param integrate_method Integration method passed to integrate_peaks():
#'                         "PD", "TS", "gauss", or "EGH".
#' @param configs          List of configuration objects from
#'                         default_ensemble_configs() or custom list.
#' @param verbose          Print progress.
#'
#' @return List with elements $areas, $summary, $baselines, $config_names.
baseline_ensemble <- function(signal, RT, peaks,
                               integrate_method = "TS",
                               configs          = default_ensemble_configs(),
                               verbose          = FALSE,
                               signal_full      = NULL,
                               RT_full          = NULL) {
  n_cfg   <- length(configs)
  n_peaks <- nrow(peaks)
  area_col <- paste0("area_", integrate_method)

  # If full-resolution signal provided, integrate there; otherwise integrate on signal.
  use_full <- !is.null(signal_full) && !is.null(RT_full)

  areas     <- matrix(NA_real_, nrow = n_peaks,  ncol = n_cfg)
  baselines <- matrix(NA_real_, nrow = length(signal), ncol = n_cfg)  # stored compressed
  cfg_names <- vapply(configs, `[[`, "", "name")

  cfg_results <- mclapply(seq_len(n_cfg), function(k) {
    if (verbose) message(sprintf("  [ensemble] %d/%d  %s", k, n_cfg, cfg_names[k]))
    bl_comp <- configs[[k]]$fun(signal)
    ar <- rep(NA_real_, n_peaks)
    if (n_peaks > 0L) {
      if (use_full) {
        # Decompress baseline to full resolution via RT interpolation
        bl_int  <- approx(RT, bl_comp, xout = RT_full, rule = 2)$y
        int_sig <- signal_full
        int_RT  <- RT_full
      } else {
        bl_int  <- bl_comp
        int_sig <- signal
        int_RT  <- RT
      }
      res <- tryCatch(
        integrate_peaks(int_sig, int_RT, peaks,
                        baseline = bl_int,
                        methods  = integrate_method),
        error = function(e) NULL)
      if (!is.null(res) && area_col %in% names(res))
        ar <- res[[area_col]]
    }
    list(bl = bl_comp, ar = ar)
  }, mc.cores = getOption("mc.cores", 4L), mc.preschedule = FALSE)

  for (k in seq_len(n_cfg)) {
    baselines[, k] <- cfg_results[[k]]$bl
    areas[, k]     <- cfg_results[[k]]$ar
  }

  # ---- summary per peak ----------------------------------------------------
  if (n_peaks > 0L) {
    area_mean  <- rowMeans(areas, na.rm = TRUE)
    area_sd    <- apply(areas, 1, sd,  na.rm = TRUE)
    area_min   <- apply(areas, 1, min, na.rm = TRUE)
    area_max   <- apply(areas, 1, max, na.rm = TRUE)
    area_cv    <- ifelse(area_mean > 0, 100 * area_sd / area_mean, NA_real_)

    cons       <- two_step_consensus(areas)

    summary_df <- data.frame(
      peak_no    = seq_len(n_peaks),
      RT         = peaks$RT,
      height     = peaks$height,
      area_mean  = area_mean,
      area_sd    = area_sd,
      area_cv_pct= area_cv,
      area_min   = area_min,
      area_max   = area_max,
      area_range = area_max - area_min,
      area_cons  = cons$area_cons,
      u_step1    = cons$u_step1,
      Z_spread   = cons$Z_spread,
      penalty    = cons$penalty,
      u_cons     = cons$u_cons,
      U_cons     = cons$U_cons
    )
  } else {
    summary_df <- data.frame()
  }

  colnames(areas) <- cfg_names
  colnames(baselines) <- cfg_names

  list(areas        = areas,
       summary      = summary_df,
       baselines    = baselines,
       config_names = cfg_names)
}


# =============================================================================
# Two-step consensus combiner
# =============================================================================

#' Combine K baseline-config area estimates using a two-step approach:
#'
#' Step 1 — correlation-based config weights:
#'   Each config's weight = mean pairwise correlation with all other configs
#'   across the peak-area vector.  Globally rogue baselines (distort all peaks
#'   together) receive low weight.  Produces the consensus point estimate â_p
#'   and a first-pass uncertainty u1_p (weighted SD).
#'
#' Step 2 — Z-score penalty for locally unstable peaks:
#'   Z_{p,k} = (a_{p,k} − â_p) / u1_p anchors deviations to the step-1
#'   result (removes circularity vs. anchoring to the unweighted mean).
#'   The weighted mean |Z| per peak (using step-1 weights) measures local
#'   integration difficulty: peaks on sloping baselines or fused pairs show
#'   high spread even among the trusted configs.  The step-1 uncertainty is
#'   inflated by penalty = max(1, weighted_mean|Z|), which equals 1 for
#'   well-resolved peaks and grows proportionally for difficult ones.
#'
#' @param areas  n_peaks × K numeric matrix (one column per baseline config).
#' @return Data frame: peak_no, area_cons, u_step1, Z_spread, penalty,
#'         u_cons (penalised), U_cons (k = 2).
two_step_consensus <- function(areas) {
  n_peaks <- nrow(areas)
  K       <- ncol(areas)

  # fallback: unweighted mean when there are too few peaks or configs
  if (K < 2L || n_peaks < 2L) {
    a_mean <- rowMeans(areas, na.rm = TRUE)
    a_sd   <- if (K > 1L) apply(areas, 1, sd, na.rm = TRUE) else rep(NA_real_, n_peaks)
    return(data.frame(peak_no   = seq_len(n_peaks),
                      area_cons = a_mean,
                      u_step1   = a_sd,
                      Z_spread  = NA_real_,
                      penalty   = 1,
                      u_cons    = a_sd,
                      U_cons    = 2 * a_sd))
  }

  # --- Step 1: correlation-based config weights ----------------------------
  C <- cor(areas, use = "pairwise.complete.obs")
  diag(C) <- NA_real_
  w <- colMeans(C, na.rm = TRUE)
  w[is.na(w) | w < 0] <- 0
  if (sum(w) == 0) w <- rep(1, K)
  w <- w / sum(w)

  area_cons <- as.numeric(areas %*% w)

  u1 <- sqrt(vapply(seq_len(n_peaks), function(p)
    sum(w * (areas[p, ] - area_cons[p])^2, na.rm = TRUE), 0.0))

  # --- Step 2: Z-score penalty for locally difficult peaks -----------------
  # Denominator: step-1 weighted SD (floored to avoid division by zero)
  denom <- pmax(u1, 1e-10)

  # Z_{p,k}: deviation of config k from step-1 consensus, in units of u1_p
  Z <- sweep(areas, 1, area_cons, "-") / denom   # n_peaks × K

  # Weighted mean |Z| per peak — high when even trusted configs spread widely
  Z_spread <- vapply(seq_len(n_peaks), function(p)
    sum(w * abs(Z[p, ]), na.rm = TRUE), 0.0)

  # Z_spread is bounded ≤ 1 (Cauchy-Schwarz: weighted MAD ≤ weighted RSD).
  # Reference: Gaussian expected value = sqrt(2/π) ≈ 0.798.
  # Values above this indicate heavier-than-Gaussian tails — bimodal config
  # distributions caused by genuine integration ambiguity (fused peaks,
  # indeterminate valley).  Penalty is 1 for Gaussian-like agreement,
  # proportionally > 1 for heavy-tailed spread.
  penalty <- pmax(1, Z_spread / sqrt(2 / pi))

  data.frame(peak_no   = seq_len(n_peaks),
             area_cons = area_cons,
             u_step1   = u1,
             Z_spread  = Z_spread,
             penalty   = penalty,
             u_cons    = u1 * penalty,
             U_cons    = 2 * u1 * penalty)
}
