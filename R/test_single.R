# test_single.R
# Processing pipeline — operates on a pre-loaded list of chromatogram matrices.
#
# INPUT  chroma_list : a list of chromatogram matrices produced by an import
#                      pipeline (e.g. import_chromas.R).  Each element is an
#                      n_samples × n_channels numeric matrix with attributes:
#
#                        attr(., "retention_time") — numeric vector (minutes),
#                                                    length n_samples
#                        attr(., "source_file")    — character, file name(s)
#                        attr(., "sample_type")    — character (optional)
#
#        chroma_list must exist in the calling environment before this script
#        is sourced.  Example:
#
#          source("./import/import_chromas.R")
#          raw <- import_chromas("../data/", "^UNIS_FID1061008\\.cdf$")
#          chroma_list <- list(raw)
#          source("./test_single.R")
#
# Run from the R/ directory.

options(mc.cores = 1L)

source("./config/read_config.R")
source("./compress/compress.R")
source("./inhibit/inhibit.R")
source("./atsa/atsa.R")
source("./sass/sass.R")
source("./arpls/doArPLS.R")
source("./summarise/summarise.R")

# --- Load method config ------------------------------------------------------
cfg <- read_method_config("./method.yaml")
cat(sprintf("  Method : %s\n", cfg$method$name))

x_factor        <- cfg$processing$x_factor
inhibit_regions <- cfg$inhibit
denoise         <- cfg$processing$denoise

# --- Validate input ----------------------------------------------------------
if (!exists("chroma_list") || !is.list(chroma_list) || length(chroma_list) == 0L)
  stop("chroma_list not found or empty.\n",
       "  Provide a list of chromatogram matrices before sourcing this script.")

df          <- chroma_list[[1]]
RT          <- attr(df, "retention_time")
source_file <- attr(df, "source_file")

if (is.null(RT))
  stop("chroma_list[[1]] is missing the 'retention_time' attribute.")
if (is.null(source_file))
  source_file <- paste0("chrom_", seq_len(ncol(df)))

cat(sprintf("  Input  : %s\n", paste(source_file, collapse = ", ")))
cat(sprintf("  %d column(s), %d samples, %.2f \u2013 %.2f min\n",
            ncol(df), nrow(df), min(RT), max(RT)))

# --- Downsample --------------------------------------------------------------
lt  <- nrow(df) - nrow(df) %% x_factor
df2 <- compress_mat(df[1:lt, , drop = FALSE], c(x_factor, 1L))
RT2 <- compress_vec(RT[1:lt], x_factor)
rm(df); gc()
cat(sprintf("  Compressed to %d samples (x%d)\n", nrow(df2), x_factor))

# --- ATSA --------------------------------------------------------------------
message("=== ATSA ===")
atsa_result <- atsa(df2, RT2,
                    segment_size     = 3.0,
                    shift_max        = 0.5,
                    smooth_width     = 11L,
                    amp_thresh       = 0,
                    baseline_correct = FALSE,
                    ref_peaks_RT     = vapply(cfg$peaks, `[[`, 0.0, "RT_ref"),
                    ref_peaks_width  = min(vapply(cfg$peaks, `[[`, 0.0, "RT_window")),
                    verbose          = TRUE)
df2 <- atsa_result$aligned
rm(atsa_result)

# --- Inhibit -----------------------------------------------------------------
message("=== Inhibit ===")
inh <- inhibit(RT2, inhibit_regions = inhibit_regions)
cat("  Good segments:", length(inh$segments), "\n")

# --- Optional SASS denoising -------------------------------------------------
df2_presass <- NULL
if (denoise) {
  message("=== SASS denoising ===")
  d_dn <- cfg$denoising$d;  fc_dn <- cfg$denoising$fc
  K_dn <- cfg$denoising$K;  lam_dn <- cfg$denoising$lam
  df2_presass <- df2                     # save pre-SASS for u_sass computation
  df2_dn <- df2
  for (seg in inh$segments) {
    idx_seg <- seg["begin"]:seg["end"]
    for (col in seq_len(ncol(df2))) {
      res <- sass_L1(as.numeric(df2[idx_seg, col]), d_dn, fc_dn, K_dn, lam_dn)
      df2_dn[idx_seg, col] <- as.numeric(res$x)
    }
  }
  df2 <- df2_dn; rm(df2_dn)
  cat("  Denoising done.\n")
}

# --- arPLS baseline ----------------------------------------------------------
message("=== arPLS baseline ===")
df3 <- matrix(NA_real_, nrow = nrow(df2), ncol = ncol(df2))
for (seg in inh$segments) {
  idx_seg <- seg["begin"]:seg["end"]
  for (col in seq_len(ncol(df2))) {
    res <- doArPLS(as.numeric(df2[idx_seg, col]),
                   lambda = cfg$baseline$lambda,
                   ratio  = cfg$baseline$ratio)
    df3[idx_seg, col] <- res$z
  }
}

# Save compressed plot data (plot built after pipeline_summary — needs peak RTs)
plot_RT  <- RT2
plot_raw <- df2[, 1]
plot_bl  <- df3[, 1]

# --- RMSE --------------------------------------------------------------------
index     <- seq_len(nrow(df3))
RMSE      <- sapply(seq_len(ncol(df3)), function(i) sqrt((df3[, i] - df2[, i])^2))
RMSE_fit  <- lapply(seq_len(ncol(RMSE)), function(i)
  loess(RMSE[, i] ~ index, span = 0.05, na.action = na.omit))
RMSE_pred <- matrix(NA_real_, nrow = nrow(RMSE), ncol = ncol(RMSE))
RMSE_pred[!inh$block, ] <- sapply(seq_len(ncol(df3)), function(i)
  predict(RMSE_fit[[i]]))
inh_block_comp <- inh$block
rm(RMSE, RMSE_fit, index, inh); gc()

# --- Decompress to full resolution -------------------------------------------
message("=== Decompressing to full resolution ===")
df2_full <- decompress_mat(df2, lt)
df3_tmp  <- df3
df3_tmp[inh_block_comp, ] <- df2[inh_block_comp, ]
rm(inh_block_comp)
df3_full     <- decompress_mat(df3_tmp, lt)
df2_presass_full <- if (!is.null(df2_presass)) decompress_mat(df2_presass[1:lt, , drop = FALSE], lt) else NULL
rm(df2_presass)
RT_full  <- RT[1:lt]
rm(df3_tmp, RT); gc()
inh_full <- inhibit(RT_full, inhibit_regions = inhibit_regions)
cat("  Decompressed:", nrow(df2_full), "samples\n")

# --- Pipeline summary --------------------------------------------------------
message("=== Pipeline summary ===")
result <- pipeline_summary(
  signal_mat        = df2_full,
  baseline_mat      = df3_full,
  RT                = RT_full,
  signal_mat_comp   = df2,
  RT_comp           = RT2,
  baseline_mat_comp = df3,
  signal_presass_mat = df2_presass_full,
  atsa_result       = NULL,
  inh               = inh_full,
  RMSE_pred         = NULL,
  col_names         = source_file,
  amp_thresh        = cfg$detection$amp_thresh,
  smooth_width      = cfg$detection$smooth_width_factor * x_factor,
  integrate_methods = cfg$integration$methods,
  run_ensemble      = cfg$integration$run_ensemble,
  ensemble_method   = cfg$integration$ensemble_method,
  method_config     = cfg,
  verbose           = TRUE
)

print(result)

# --- Comparison plot ---------------------------------------------------------
message("=== Saving comparison plot ===")
peak_RTs  <- result$peaks$RT
corr_plot <- 100 * (plot_raw - replace(plot_bl, is.na(plot_bl),
                                        plot_raw[is.na(plot_bl)]))
png("test_single_comparison.png", width = 1200, height = 500)
matplot(plot_RT, cbind(plot_raw, plot_bl, corr_plot), type = "l", lty = 1,
        col = c("steelblue", "firebrick", "darkgreen"), lwd = c(1, 1.5, 1),
        ylim = c(-30, 300),
        xlab = "Retention time (min)", ylab = "Signal",
        main = sprintf("arPLS: raw (blue), baseline (red), 100\u00d7 corrected (green)  [%s]",
                       cfg$method$name))
abline(h = 0, col = "grey60", lty = 2)
if (length(peak_RTs) > 0) {
  abline(v = peak_RTs, col = "orange", lty = 2, lwd = 1)
  if ("name" %in% names(result$peaks)) {
    nms <- ifelse(is.na(result$peaks$name),
                  sprintf("RT%.2f", peak_RTs), result$peaks$name)
    axis(3, at = peak_RTs, labels = nms, cex.axis = 0.7, las = 2)
  }
}
legend("topright",
       legend = c("Raw (df2)", "arPLS baseline (df3)",
                  "100\u00d7 Corrected (df2-df3)", "Peak apex"),
       col = c("steelblue", "firebrick", "darkgreen", "orange"),
       lty = c(1, 1, 1, 2), lwd = c(1, 1.5, 1, 1))
dev.off()
rm(plot_RT, plot_raw, plot_bl, corr_plot)
cat("  Plot saved: test_single_comparison.png\n")

rm(df2, RT2, df3)
message("=== Done ===")
