# summarise.R
# Pipeline summary: detect peaks, integrate, estimate uncertainty and LOD/LOQ
# for every chromatogram column, then print a structured report.
#
# Requires (source before calling):
#   source("./atsa/peak_detect.R")
#   source("./integrate/peak_integrate.R")
#   source("./integrate/baseline_ensemble.R")   # only if run_ensemble = TRUE
#   source("./integrate/lod_loq.R")
#
# Main function:
#   result <- pipeline_summary(signal_mat, baseline_mat, RT, ...)
#
# All heavy work (peak detection, integration, ensemble, LOD/LOQ) is done
# inside this function.  Pass the objects you already have from read_chroma.R:
#   signal_mat   = df2        (aligned signal, before baseline subtraction)
#   baseline_mat = df3        (SASS baseline)
#   RT           = RT2
#   atsa_result, inh, RMSE_pred are optional but enrich the report.

source("./atsa/peak_detect.R")
source("./integrate/peak_integrate.R")
source("./integrate/baseline_ensemble.R")
source("./integrate/lod_loq.R")
source("./config/read_config.R")
source("./result/chrom_result.R")
library(parallel)


# =============================================================================
# Internal formatting helpers
# =============================================================================

.hline <- function(char = "-", width = 72) paste(rep(char, width), collapse = "")

.section <- function(title) {
  cat("\n", .hline("="), "\n", title, "\n", .hline("="), "\n", sep = "")
}

.subsection <- function(title) {
  cat("\n  ", .hline("-", 68), "\n  ", title, "\n  ", .hline("-", 68), "\n",
      sep = "")
}

.fmt <- function(x, digits = 3) {
  if (is.na(x)) return("NA")
  formatC(x, digits = digits, format = "g")
}

.print_df <- function(df, indent = "  ") {
  # Print a data frame as a fixed-width table with indentation.
  if (nrow(df) == 0L) { cat(indent, "<empty>\n"); return(invisible(NULL)) }
  lines <- capture.output(print(df, row.names = FALSE))
  cat(paste0(indent, lines, "\n"), sep = "")
}


# =============================================================================
# Main function
# =============================================================================

#' Run and summarise the full downstream pipeline for a set of chromatograms.
#'
#' @param signal_mat      n_samples × n_chroms matrix (aligned signal, e.g. df2).
#' @param baseline_mat    n_samples × n_chroms matrix (fitted baseline, e.g. df3).
#' @param RT              Retention time vector in minutes (length = n_samples).
#' @param atsa_result     List returned by atsa(). Optional; used for alignment stats.
#' @param inh             List returned by inhibit(). Optional; used to report regions.
#' @param RMSE_pred       n_samples × n_chroms RMSE matrix from read_chroma.R. Optional.
#' @param smooth_width    Smoothing window (samples) for peak detection.
#' @param amp_thresh      Minimum peak height for detection.
#' @param slope_thresh    Minimum derivative slope for peak detection.
#' @param integrate_methods  Character vector: any of "PD","TS","gauss","EGH".
#' @param run_ensemble    Logical. Run baseline ensemble to estimate area uncertainty.
#'                        Slow (~12 baselines × n_chroms). Set FALSE to skip.
#' @param ensemble_configs  Configuration list; default_ensemble_configs() if NULL.
#' @param ensemble_method   Integration method used inside the ensemble ("TS").
#' @param method_lod      Method LOD in signal units for comparison. NULL = skip.
#' @param method_loq      Method LOQ in signal units. NULL = skip.
#' @param cal_slope       Calibration slope (signal / conc unit). NULL = skip.
#' @param noise_method    Noise estimator: "mad", "sd", or "rms".
#' @param col_names       Optional character vector of column labels.
#' @param verbose         Print section headers while running.
#'
#' @return Named list:
#'   $data_info     – scalar summary of inputs
#'   $alignment     – per-column alignment stats (if atsa_result given)
#'   $baseline_rmse – per-column RMSE summary (if RMSE_pred given)
#'   $peaks         – list of per-column peak data frames
#'   $integration   – list of per-column integration data frames
#'   $uncertainty   – list of per-column ensemble summary data frames
#'   $lod_loq       – data frame, one row per column
#'   $report        – character vector of printed report lines (invisible)
pipeline_summary <- function(signal_mat,
                              baseline_mat,
                              RT,
                              signal_mat_comp    = NULL,
                              RT_comp            = NULL,
                              baseline_mat_comp  = NULL,
                              signal_presass_mat = NULL,
                              atsa_result       = NULL,
                              inh              = NULL,
                              RMSE_pred        = NULL,
                              smooth_width     = 11L,
                              amp_thresh       = 0,
                              slope_thresh     = 0,
                              integrate_methods = c("PD", "TS", "gauss", "EGH"),
                              run_ensemble     = TRUE,
                              ensemble_configs = NULL,
                              ensemble_method  = "TS",
                              method_lod       = NULL,
                              method_loq       = NULL,
                              cal_slope        = NULL,
                              noise_method     = "mad",
                              col_names        = NULL,
                              method_config    = NULL,
                              verbose          = TRUE) {

  n_samp  <- nrow(signal_mat)
  n_chrom <- ncol(signal_mat)
  dt      <- RT[2] - RT[1]

  if (is.null(col_names))
    col_names <- paste0("chrom_", seq_len(n_chrom))

  if (is.null(ensemble_configs))
    ensemble_configs <- default_ensemble_configs()

  # Capture all cat() output for $report
  report_lines <- character(0)
  .cat <- function(...) {
    txt <- paste0(...)
    cat(txt)
    report_lines <<- c(report_lines, txt)
  }

  # ---- 0. Compute baseline-corrected signal ---------------------------------
  # NA in baseline_mat marks inhibited regions; treat as zero correction there.
  baseline_safe <- baseline_mat
  baseline_safe[is.na(baseline_safe)] <- 0
  signal_residual <- signal_mat - baseline_safe   # raw residual (signed) for noise estimation
  signal_bc       <- signal_residual
  signal_bc[signal_bc < 0] <- 0                  # clipped signal for peak detection / integration
  # Zero out inhibited regions so peaks are only detected in SASS-corrected segments.
  if (!is.null(inh) && length(inh$block) == nrow(signal_bc))
    signal_bc[inh$block, ] <- 0

  # ---- 1. Peak detection per column ----------------------------------------
  if (verbose) message("[summary] Detecting peaks ...")
  peaks_list <- mclapply(seq_len(n_chrom), function(j) {
    pk <- find_peaks(signal_bc[, j], RT,
                     smooth_width = smooth_width,
                     amp_thresh   = amp_thresh,
                     slope_thresh = slope_thresh)
    # Drop peaks that fall in SASS baseline-overshoot regions.
    # Check is done on the compressed corrected signal (block-average faithfully
    # captures overshoot; full-res noise spikes above a smooth interpolated
    # baseline would otherwise pass the filter).
    if (nrow(pk) > 0L && !is.null(signal_mat_comp) && !is.null(baseline_mat_comp)) {
      bl_safe  <- baseline_mat_comp[, j]
      bl_safe[is.na(bl_safe)] <- signal_mat_comp[is.na(bl_safe), j]  # inhibit → 0 residual
      corr_comp <- signal_mat_comp[, j] - bl_safe
      pk <- pk[sapply(seq_len(nrow(pk)), function(i) {
        ci  <- which.min(abs(RT_comp - pk$RT[i]))
        win <- max(1L, ci - 1L):min(length(corr_comp), ci + 1L)
        max(corr_comp[win]) > 0
      }), ]
    }
    pk
  }, mc.cores = getOption("mc.cores", 4L))

  # ---- 2. Integration per column -------------------------------------------
  if (verbose) message("[summary] Integrating peaks ...")
  integ_list <- mclapply(seq_len(n_chrom), function(j)
    tryCatch(
      integrate_peaks(signal_bc[, j], RT, peaks_list[[j]],
                      baseline = NULL,
                      methods  = integrate_methods),
      error = function(e) data.frame()),
    mc.cores = getOption("mc.cores", 4L))

  # ---- 3. Baseline ensemble (uncertainty) ----------------------------------
  ens_list <- vector("list", n_chrom)
  if (run_ensemble) {
    if (verbose) message("[summary] Running baseline ensemble ...")
    # Ensemble baseline parameters (LMV span, airPLS/arPLS lambda) were
    # calibrated for compressed signal (~24 samples/min).  Running at full
    # resolution (~1200 samples/min) makes LMV spans absurdly tight — they
    # follow the peak shape and subtract peak area as baseline, producing a
    # ~9× spread in areas (CV ~50%).
    #
    # Strategy: fit baselines on COMPRESSED SASS-corrected signal; decompress
    # via approx() to full resolution; integrate on full-resolution signal.
    # This keeps span parameters meaningful while preserving integration accuracy.
    #
    # Running on raw compressed signal is avoided: the solvent-front spike
    # dominates and algorithms handle it very differently.  Instead we use
    # signal_bc_comp = compressed (raw − SASS baseline), clipped ≥ 0.

    use_comp_ens <- !is.null(signal_mat_comp) &&
                    !is.null(baseline_mat_comp) &&
                    !is.null(RT_comp)

    if (use_comp_ens) {
      # Build compressed SASS-corrected signal.
      # NA in baseline_mat_comp marks inhibit regions: treat as zero residual
      # (set bl = signal so signal − bl = 0, no need for a separate inh lookup).
      bl_comp_safe <- baseline_mat_comp
      bl_comp_safe[is.na(bl_comp_safe)] <-
        signal_mat_comp[is.na(bl_comp_safe)]
      signal_bc_comp <- signal_mat_comp - bl_comp_safe
      signal_bc_comp[signal_bc_comp < 0] <- 0
    }

    ens_list <- mclapply(seq_len(n_chrom), function(j) {
      if (verbose) message(sprintf("  column %d / %d", j, n_chrom))
      if (use_comp_ens) {
        tryCatch(
          baseline_ensemble(
            signal           = signal_bc_comp[, j],
            RT               = RT_comp,
            peaks            = peaks_list[[j]],
            integrate_method = ensemble_method,
            configs          = ensemble_configs,
            verbose          = FALSE,
            signal_full      = signal_bc[, j],
            RT_full          = RT),
          error = function(e) NULL)
      } else {
        # Fallback: no compressed inputs — run on full-res signal_bc.
        # CV may be inflated if LMV spans are too tight for full resolution.
        tryCatch(
          baseline_ensemble(
            signal           = signal_bc[, j],
            RT               = RT,
            peaks            = peaks_list[[j]],
            integrate_method = ensemble_method,
            configs          = ensemble_configs,
            verbose          = FALSE,
            signal_full      = NULL,
            RT_full          = NULL),
          error = function(e) NULL)
      }
    }, mc.cores = getOption("mc.cores", 4L))
  }

  # ---- 4. LOD / LOQ --------------------------------------------------------
  if (verbose) message("[summary] Computing LOD/LOQ ...")
  ll_df <- lod_loq_batch(signal_residual, RT, peaks_list,
                          method_lod   = method_lod,
                          method_loq   = method_loq,
                          cal_slope    = cal_slope,
                          noise_method = noise_method)
  ll_df$chrom_name <- col_names

  # ---- 4b. u_noise per peak ------------------------------------------------
  # Noise propagation into the integrated area: integrating N = width/dt
  # independent noise samples each with SD sigma_noise gives
  #   u_noise = sigma_noise * sqrt(width * dt)
  # Both width and dt are in minutes here; area_scale converts to counts·s.
  u_noise_list <- lapply(seq_len(n_chrom), function(j) {
    sigma_j <- ll_df$noise_sd[j]
    pk_j    <- peaks_list[[j]]
    if (is.na(sigma_j) || nrow(pk_j) == 0L)
      return(rep(NA_real_, nrow(pk_j)))
    sigma_j * sqrt(pk_j$width * dt)   # in signal × min; scaled below
  })

  # ---- 4c. u_sass per peak -------------------------------------------------
  # SASS denoising modifies the signal within peak windows.  The area of the
  # signed difference (pre-SASS minus post-SASS) integrated over each peak
  # is the direct area shift introduced by denoising.  Its absolute value is
  # treated as a one-sided uncertainty component; it scales with peak amplitude.
  u_sass_list <- lapply(seq_len(n_chrom), function(j) {
    pk_j <- peaks_list[[j]]
    if (is.null(signal_presass_mat) || nrow(pk_j) == 0L)
      return(rep(NA_real_, nrow(pk_j)))
    diff_j <- signal_presass_mat[, j] - signal_mat[, j]   # pre − post SASS
    vapply(seq_len(nrow(pk_j)), function(p) {
      lo <- if (!is.null(pk_j$lo_idx) && !is.na(pk_j$lo_idx[p])) pk_j$lo_idx[p]
            else max(1L, pk_j$idx[p] - round(pk_j$width[p] / dt / 2))
      hi <- if (!is.null(pk_j$hi_idx) && !is.na(pk_j$hi_idx[p])) pk_j$hi_idx[p]
            else min(n_samp, pk_j$idx[p] + round(pk_j$width[p] / dt / 2))
      abs(.trapz(RT, diff_j, lo, hi))   # signal × min; scaled below
    }, 0.0)
  })

  # ---- 5. Area unit conversion ---------------------------------------------
  # Integration is performed over RT in whatever unit RT is expressed (here:
  # minutes).  Empower and most CDS report peak areas in counts·s.
  # Apply area_scale (60 when RT_unit = "min") to all area columns so that
  # the printed report and the returned ChromResult are in counts·s.
  area_scale <- if (!is.null(method_config)) method_config$area_scale else 1

  if (area_scale != 1) {
    area_int_cols <- c("area_PD", "area_TS", "area_gauss", "area_EGH")
    integ_list <- lapply(integ_list, function(ir) {
      if (is.null(ir) || nrow(ir) == 0L) return(ir)
      for (col in intersect(area_int_cols, names(ir)))
        ir[[col]] <- ir[[col]] * area_scale
      ir
    })
    area_ens_cols <- c("area_cons", "u_step1", "u_cons", "U_cons",
                       "area_mean", "area_sd",
                       "area_min",  "area_max", "area_range")
    ens_list <- lapply(ens_list, function(en) {
      if (is.null(en) || is.null(en$summary) || nrow(en$summary) == 0L)
        return(en)
      for (col in intersect(area_ens_cols, names(en$summary)))
        en$summary[[col]] <- en$summary[[col]] * area_scale
      # area_cv_pct, Z_spread, penalty are dimensionless — not scaled
      en
    })
    u_noise_list <- lapply(u_noise_list, function(u) u * area_scale)
    u_sass_list  <- lapply(u_sass_list,  function(u) u * area_scale)
  }

  # ==========================================================================
  # Print report
  # ==========================================================================
  .section("CHROMATOGRAPHY PIPELINE SUMMARY")

  # ---- Data info -----------------------------------------------------------
  .subsection("1. Data")
  .cat(sprintf("  Chromatograms    : %d\n", n_chrom))
  .cat(sprintf("  Samples / chrom  : %d\n", n_samp))
  .cat(sprintf("  RT range         : %.3f – %.3f min\n", min(RT), max(RT)))
  .cat(sprintf("  Time step        : %.4f min  (%.1f Hz)\n", dt, 1 / dt))

  if (!is.null(inh)) {
    n_inh <- sum(inh$block)
    .cat(sprintf("  Inhibited points : %d  (%.1f%% of total)\n",
                 n_inh, 100 * n_inh / n_samp))
    if (length(inh$segments) > 0) {
      .cat(sprintf("  Good segments    : %d\n", length(inh$segments)))
    }
  }

  # ---- Alignment -----------------------------------------------------------
  if (!is.null(atsa_result)) {
    .subsection("2. Time-shift Alignment (ATSA)")
    ref  <- atsa_result$reference_idx
    shft <- atsa_result$shifts          # n_seg × n_chrom
    .cat(sprintf("  Reference column : %s  (%d peaks)\n",
                 col_names[ref], nrow(atsa_result$ref_peaks)))
    .cat(sprintf("  Segments         : %d\n", nrow(atsa_result$segments)))

    mean_shift_min <- apply(shft, 2, function(s) mean(abs(s), na.rm = TRUE)) * dt
    align_df <- data.frame(
      chrom          = col_names,
      mean_abs_shift_min = round(mean_shift_min, 4),
      max_abs_shift_min  = round(apply(shft, 2,
                                       function(s) max(abs(s), na.rm = TRUE)) * dt, 4)
    )
    .cat("\n")
    .print_df(align_df)
  }

  # ---- Baseline RMSE -------------------------------------------------------
  if (!is.null(RMSE_pred)) {
    .subsection("3. Baseline Fit Quality (SASS RMSE)")
    rmse_summary <- data.frame(
      chrom    = col_names,
      RMSE_mean = round(colMeans(RMSE_pred, na.rm = TRUE), 4),
      RMSE_max  = round(apply(RMSE_pred, 2, max, na.rm = TRUE), 4),
      RMSE_p95  = round(apply(RMSE_pred, 2,
                              function(x) quantile(x, 0.95, na.rm = TRUE)), 4)
    )
    .cat("\n")
    .print_df(rmse_summary)
  }

  # ---- Peak detection summary ---------------------------------------------
  .subsection("4. Peak Detection")
  n_peaks_vec <- vapply(peaks_list, nrow, 1L)
  pk_summary <- data.frame(
    chrom    = col_names,
    n_peaks  = n_peaks_vec,
    RT_first = vapply(peaks_list, function(p)
                 if (nrow(p) > 0) round(min(p$RT), 3) else NA_real_, 1.0),
    RT_last  = vapply(peaks_list, function(p)
                 if (nrow(p) > 0) round(max(p$RT), 3) else NA_real_, 1.0),
    height_max = vapply(peaks_list, function(p)
                 if (nrow(p) > 0) round(max(p$height), 3) else NA_real_, 1.0)
  )
  .cat("\n")
  .print_df(pk_summary)

  # Shared peak-RT table (reference column only, as anchor)
  if (!is.null(atsa_result)) {
    ref_pks <- atsa_result$ref_peaks
    if (nrow(ref_pks) > 0) {
      .cat(sprintf("\n  Reference peaks (%s):\n", col_names[atsa_result$reference_idx]))
      ref_tbl <- data.frame(
        peak_no = seq_len(nrow(ref_pks)),
        RT_min  = round(ref_pks$RT, 4),
        height  = round(ref_pks$height, 3),
        width_min = round(ref_pks$width, 4)
      )
      .print_df(ref_tbl)
    }
  }

  # ---- Integration ---------------------------------------------------------
  .subsection("5. Peak Integration")
  .cat(sprintf("  Methods: %s\n\n", paste(integrate_methods, collapse = ", ")))

  for (j in seq_len(n_chrom)) {
    ir <- integ_list[[j]]
    if (is.null(ir) || nrow(ir) == 0L) next
    .cat(sprintf("  %s  (%d peaks)\n", col_names[j], nrow(ir)))

    # Round numeric columns for display
    disp <- ir
    num_cols <- sapply(disp, is.numeric)
    disp[num_cols] <- lapply(disp[num_cols], round, digits = 4)
    .print_df(disp, indent = "    ")
    .cat("\n")
  }

  # ---- Baseline uncertainty ------------------------------------------------
  if (run_ensemble && !all(vapply(ens_list, is.null, TRUE))) {
    .subsection("6. Peak Area Uncertainty (Baseline Ensemble)")
    .cat(sprintf("  Configurations: %d  |  Integration method: %s\n\n",
                 length(ensemble_configs), ensemble_method))

    for (j in seq_len(n_chrom)) {
      en <- ens_list[[j]]
      if (is.null(en) || nrow(en$summary) == 0L) next
      .cat(sprintf("  %s\n", col_names[j]))
      disp <- en$summary[, c("peak_no", "RT", "area_cons", "U_cons",
                              "Z_spread", "penalty", "area_cv_pct")]
      u_noise_j <- u_noise_list[[j]]
      u_sass_j  <- u_sass_list[[j]]
      if (!is.null(u_noise_j) && length(u_noise_j) == nrow(disp)) {
        disp$u_noise <- u_noise_j
        disp$u_sass  <- if (!is.null(u_sass_j) && length(u_sass_j) == nrow(disp))
                          u_sass_j else NA_real_
        u_area_j     <- sqrt(en$summary$u_cons^2 +
                             u_noise_j^2 +
                             replace(u_sass_j^2, is.na(u_sass_j), 0))
        disp$U_area  <- round(2 * u_area_j, 3)
      }
      num_cols <- sapply(disp, is.numeric)
      disp[num_cols] <- lapply(disp[num_cols], round, digits = 3)
      .print_df(disp, indent = "    ")
      .cat("\n")
    }
  }

  # ---- LOD / LOQ -----------------------------------------------------------
  .subsection("7. LOD / LOQ")
  .cat(sprintf("  Noise estimator : %s\n", noise_method))
  if (!is.null(method_lod))
    .cat(sprintf("  Method LOD      : %s signal units\n", .fmt(method_lod)))
  if (!is.null(method_loq))
    .cat(sprintf("  Method LOQ      : %s signal units\n", .fmt(method_loq)))
  if (!is.null(cal_slope))
    .cat(sprintf("  Calibration slope: %s\n", .fmt(cal_slope)))
  .cat("\n")

  # LOD_conc, LOQ_conc, LOD_ok, LOQ_ok belong to the quality gate that runs
  # after calibration — omit them from the per-run report and result object.
  ll_disp <- ll_df[, intersect(
    c("chrom_name", "noise_sd", "LOD", "LOQ",
      "n_peaks_above_LOD", "n_peaks_above_LOQ"),
    names(ll_df))]
  num_cols <- sapply(ll_disp, is.numeric)
  ll_disp[num_cols] <- lapply(ll_disp[num_cols], round, digits = 4)
  .print_df(ll_disp)

  .cat("\n", .hline("="), "\n")

  # Strip quality-gate columns from ll_df before storing in ChromResult
  qg_cols <- c("LOD_conc", "LOQ_conc", "LOD_ok", "LOQ_ok",
                "method_LOD", "method_LOQ")
  ll_df <- ll_df[, setdiff(names(ll_df), qg_cols), drop = FALSE]

  # ==========================================================================
  # Build unified peaks table
  # ==========================================================================
  peaks_rows <- lapply(seq_len(n_chrom), function(j) {
    ir <- integ_list[[j]]
    if (is.null(ir) || nrow(ir) == 0L) return(NULL)

    df_out <- cbind(
      data.frame(chrom = col_names[j], stringsAsFactors = FALSE),
      ir)

    # u_noise: noise floor propagated through integration
    u_noise_j <- u_noise_list[[j]]
    df_out$u_noise <- if (!is.null(u_noise_j) && length(u_noise_j) == nrow(ir))
                        u_noise_j else rep(NA_real_, nrow(ir))

    # u_sass: area shift from SASS denoising within peak window
    u_sass_j <- u_sass_list[[j]]
    df_out$u_sass <- if (!is.null(u_sass_j) && length(u_sass_j) == nrow(ir))
                       u_sass_j else rep(NA_real_, nrow(ir))

    en <- ens_list[[j]]
    if (run_ensemble && !is.null(en) &&
        !is.null(en$summary) && nrow(en$summary) > 0L) {
      ens_cols <- intersect(
        c("area_cons", "u_step1", "Z_spread", "penalty", "u_cons", "U_cons",
          "area_mean", "area_sd", "area_cv_pct",
          "area_min",  "area_max", "area_range"),
        names(en$summary))
      df_out <- cbind(df_out, en$summary[, ens_cols, drop = FALSE])
    }

    # Combined uncertainty: baseline ambiguity + noise + SASS in quadrature
    if ("u_cons" %in% names(df_out)) {
      df_out$u_area <- sqrt(df_out$u_cons^2 +
                            df_out$u_noise^2 +
                            replace(df_out$u_sass^2, is.na(df_out$u_sass), 0))
      df_out$U_area <- 2 * df_out$u_area
    }

    df_out
  })
  peaks_combined <- do.call(rbind, peaks_rows)
  if (is.null(peaks_combined))
    peaks_combined <- data.frame()
  rownames(peaks_combined) <- NULL

  # Assign peak names from method config (if supplied)
  if (!is.null(method_config) && nrow(peaks_combined) > 0L)
    peaks_combined <- match_peak_names(peaks_combined, method_config)

  # ==========================================================================
  # Return ChromResult
  # ==========================================================================
  meta <- list(
    n_chrom           = n_chrom,
    n_samp            = n_samp,
    RT_range          = range(RT),
    dt_min            = dt,
    col_names         = col_names,
    method_name       = if (!is.null(method_config)) method_config$method$name else NULL,
    integrate_methods = integrate_methods,
    ensemble_method   = ensemble_method,
    n_ensemble_configs= length(if (!is.null(ensemble_configs))
                                 ensemble_configs else default_ensemble_configs())
  )

  invisible(.new_ChromResult(
    meta          = meta,
    peaks         = peaks_combined,
    lod_loq       = ll_df,
    alignment     = if (!is.null(atsa_result)) align_df    else NULL,
    baseline_rmse = if (!is.null(RMSE_pred))  rmse_summary else NULL,
    report        = report_lines
  ))
}
