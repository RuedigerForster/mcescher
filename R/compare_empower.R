# compare_empower.R
# Run mcescher on one sample directory and write results to CSV for comparison
# with Empower (CDS) peak table extracted by compare_empower.py.
#
# Usage (from R/ directory):
#   Rscript compare_empower.R <data_dir>  e.g.
#   Rscript compare_empower.R ../data/UNIS_FID_9000002_908779/

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "../data/UNIS_FID_9000002_908779/"
if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
data_dir <- paste0(normalizePath(data_dir, mustWork = TRUE), "/")  # absolute + trailing slash

setwd("/home/nvidiauser/Developments/mcescher/R")

options(mc.cores = parallel::detectCores(logical = FALSE))
suppressPackageStartupMessages({
  source("./config/read_config.R")
  source("./import/import_chromas.R")
  source("./compress/compress.R")
  source("./inhibit/inhibit.R")
  source("./atsa/atsa.R")
  source("./sass/sass.R")
  source("./arpls/doArPLS.R")
  source("./summarise/summarise.R")
})

cfg            <- read_method_config("./method.yaml")
x_factor       <- cfg$processing$x_factor
inhibit_regions <- cfg$inhibit
denoise        <- cfg$processing$denoise

cat(sprintf("Loading CDFs from: %s\n", data_dir))
df  <- import_chromas(data_dir, "[.]cdf$")
RT  <- attr(df, "retention_time")
src <- attr(df, "source_file")
cat(sprintf("  %d chromatograms, %d samples, %.2f-%.2f min\n",
            ncol(df), nrow(df), min(RT), max(RT)))

# ── Compress ──────────────────────────────────────────────────────────────────
lt  <- nrow(df) - nrow(df) %% x_factor
df2 <- compress_mat(df[1:lt, , drop = FALSE], c(x_factor, 1L))
RT2 <- compress_vec(RT[1:lt], x_factor)
rm(df); gc()

# ── ATSA ──────────────────────────────────────────────────────────────────────
cat("ATSA alignment...\n")
n_comp_per_min <- 1 / (RT2[2] - RT2[1])            # compressed sampling rate (samples/min)
lmv_span_min   <- 2.0                               # desired LMV baseline window (min)
atsa_result <- atsa(df2, RT2,
                    segment_size     = 3.0,
                    shift_max        = 0.5,
                    smooth_width     = max(11L, round(0.1 * n_comp_per_min)),
                    amp_thresh       = 0,
                    baseline_correct = FALSE,
                    lmv_span         = max(50L, round(lmv_span_min * n_comp_per_min)),
                    ref_peaks_RT     = vapply(cfg$peaks, `[[`, 0.0, "RT_ref"),
                    ref_peaks_width  = min(vapply(cfg$peaks, `[[`, 0.0, "RT_window")),
                    verbose          = FALSE)
df2 <- atsa_result$aligned; rm(atsa_result)

# ── Inhibit ───────────────────────────────────────────────────────────────────
inh <- inhibit(RT2, inhibit_regions = inhibit_regions)

# ── SASS denoising ────────────────────────────────────────────────────────────
df2_presass <- NULL
if (denoise) {
  cat("SASS denoising...\n")
  df2_presass <- df2
  df2_dn <- df2
  for (seg in inh$segments) {
    idx_seg <- seg["begin"]:seg["end"]
    for (col in seq_len(ncol(df2)))
      df2_dn[idx_seg, col] <- as.numeric(
        sass_L1(as.numeric(df2[idx_seg, col]),
                cfg$denoising$d, cfg$denoising$fc,
                cfg$denoising$K, cfg$denoising$lam)$x)
  }
  df2 <- df2_dn; rm(df2_dn)
}

# ── arPLS baseline ────────────────────────────────────────────────────────────
cat("arPLS baseline...\n")
df3 <- matrix(NA_real_, nrow = nrow(df2), ncol = ncol(df2))
for (seg in inh$segments)
  for (col in seq_len(ncol(df2)))
    df3[seg["begin"]:seg["end"], col] <-
      doArPLS(as.numeric(df2[seg["begin"]:seg["end"], col]),
              lambda = cfg$baseline$lambda, ratio = cfg$baseline$ratio)$z

# ── Decompress ────────────────────────────────────────────────────────────────
cat("Decompressing...\n")
df2_full <- decompress_mat(df2, lt)
df3_tmp  <- df3; df3_tmp[inh$block, ] <- df2[inh$block, ]
df3_full <- decompress_mat(df3_tmp, lt)
df2_presass_full <- if (!is.null(df2_presass)) decompress_mat(df2_presass, lt) else NULL
RT_full  <- RT[1:lt]
inh_full <- inhibit(RT_full, inhibit_regions = inhibit_regions)
# Keep df2 and df3 for ensemble (must run on compressed signal so LMV spans
# are meaningful; full-res LMV span=30 is ~1.5 s and follows peak shapes).
rm(df3_tmp, df2_presass); gc()

# ── Pipeline summary ──────────────────────────────────────────────────────────
cat("Integrating peaks...\n")
result <- pipeline_summary(
  signal_mat         = df2_full,
  baseline_mat       = df3_full,
  RT                 = RT_full,
  signal_mat_comp    = df2,
  RT_comp            = RT2,
  baseline_mat_comp  = df3,
  signal_presass_mat = df2_presass_full,
  atsa_result        = NULL,
  inh                = inh_full,
  RMSE_pred          = NULL,
  col_names          = src,
  amp_thresh         = cfg$detection$amp_thresh,
  smooth_width       = cfg$detection$smooth_width_factor * x_factor,
  integrate_methods  = cfg$integration$methods,
  run_ensemble       = cfg$integration$run_ensemble,
  ensemble_method    = cfg$integration$ensemble_method,
  method_config      = cfg,
  verbose            = FALSE
)

# ── Write mcescher peaks to CSV ───────────────────────────────────────────────
out_peaks <- paste0(data_dir, "mcescher_peaks.csv")
peaks_out <- result$peaks[, c("chrom", "RT", "name", "height",
                               "area_PD", "area_TS", "area_gauss", "area_EGH",
                               "area_AIC", "aic_winner",
                               "area_cons", "u_cons", "U_cons",
                               "u_noise", "u_sass", "u_area", "U_area",
                               "area_cv_pct", "Z_spread", "penalty")]
write.csv(peaks_out, out_peaks, row.names = FALSE)
cat(sprintf("mcescher peaks written to: %s\n", out_peaks))
cat(sprintf("Total peaks: %d across %d chromatograms\n",
            nrow(peaks_out), length(unique(peaks_out$chrom))))
