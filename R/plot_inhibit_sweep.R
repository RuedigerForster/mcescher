# plot_inhibit_sweep.R
# Sweep the dead-time inhibit end point from the current value down to 0
# and record what arPLS does to the Methane / Ethane region.
#
# Usage (from R/ directory):
#   Rscript plot_inhibit_sweep.R [data_dir] [cdf_basename]

setwd("/home/nvidiauser/Developments/mcescher/R")

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else
              "../data/UNIS_FID_9000002_908779/"
cdf_file <- if (length(args) >= 2) args[2] else
              "UNIS_FID_9000002_909535.cdf"   # reference chromatogram
out_png  <- "../data/inhibit_sweep.png"

if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")

suppressPackageStartupMessages({
  source("./config/read_config.R")
  source("./import/import_chromas.R")
  source("./compress/compress.R")
  source("./inhibit/inhibit.R")
  source("./sass/sass.R")
  source("./arpls/doArPLS.R")
  source("./integrate/peak_integrate.R")
  source("./atsa/peak_detect.R")
})

cfg      <- read_method_config("./method.yaml")
x_factor <- cfg$processing$x_factor

# ── Load and compress ─────────────────────────────────────────────────────────
cat(sprintf("Loading: %s\n", cdf_file))
df  <- import_chromas(data_dir, paste0("^", cdf_file, "$"))
RT  <- attr(df, "retention_time")
lt  <- nrow(df) - nrow(df) %% x_factor
df2 <- compress_mat(df[1:lt, , drop = FALSE], c(x_factor, 1L))
RT2 <- compress_vec(RT[1:lt], x_factor)
RT_full <- RT[1:lt]
rm(df); gc()

# ── SASS denoising (fixed, using full inhibit) ────────────────────────────────
cat("SASS denoising...\n")
inh0   <- inhibit(RT2, inhibit_regions = cfg$inhibit)
sig_dn <- as.numeric(df2[, 1])
for (seg in inh0$segments) {
  idx <- seg["begin"]:seg["end"]
  sig_dn[idx] <- as.numeric(
    sass_L1(sig_dn[idx],
            cfg$denoising$d, cfg$denoising$fc,
            cfg$denoising$K, cfg$denoising$lam)$x)
}
# Also decompress SASS signal for full-resolution use
sig_dn_full <- as.numeric(decompress_mat(matrix(sig_dn, ncol = 1L), lt))

# ── Empower reference areas (for reference lines) ─────────────────────────────
EMP_METHANE <- 80742.8
EMP_ETHANE  <-  3550.4

# ── Inhibit sweep ─────────────────────────────────────────────────────────────
# Vary only the END of the first inhibit region (dead time / solvent front).
# Step from the current value down to 0 in 0.05 min increments.
inh_end_current <- cfg$inhibit[[1]][2]   # 1.10 min
inh_tail        <- cfg$inhibit[[2]]      # [20.38, 27.00] — fixed
inh_end_vals    <- seq(inh_end_current, 0, by = -0.05)

cat(sprintf("Sweeping inhibit end from %.2f to 0.00 min (%d steps)...\n",
            inh_end_current, length(inh_end_vals)))

results <- lapply(inh_end_vals, function(inh_end) {
  inh_regions <- list(c(0, inh_end), inh_tail)
  inh <- inhibit(RT2, inhibit_regions = inh_regions)

  # arPLS on compressed SASS-denoised signal
  bl <- rep(NA_real_, length(sig_dn))
  for (seg in inh$segments) {
    idx <- seg["begin"]:seg["end"]
    bl[idx] <- doArPLS(sig_dn[idx],
                       lambda = cfg$baseline$lambda,
                       ratio  = cfg$baseline$ratio)$z
  }
  bl[inh$block] <- sig_dn[inh$block]

  # Decompress baseline to full resolution
  bl_full <- as.numeric(decompress_mat(matrix(bl, ncol = 1L), lt))

  # Baseline-corrected full-res signal
  sig_bc  <- pmax(sig_dn_full - bl_full, 0)

  # Full-res inhibit for peak detection
  inh_full <- inhibit(RT_full, inhibit_regions = inh_regions)
  sig_bc[inh_full$block] <- 0

  # Detect peaks (using cfg params)
  peaks <- find_peaks(sig_bc, RT_full,
                      smooth_width = cfg$detection$smooth_width_factor * x_factor,
                      amp_thresh   = cfg$detection$amp_thresh)

  # Match Methane and Ethane by nearest RT
  get_peak <- function(rt_ref, win) {
    ok <- abs(peaks$RT - rt_ref) <= win
    if (!any(ok)) return(NULL)
    peaks[ok, ][which.min(abs(peaks$RT[ok] - rt_ref)), ]
  }
  me  <- get_peak(cfg$peaks[[1]]$RT_ref, cfg$peaks[[1]]$RT_window)
  eth <- get_peak(cfg$peaks[[2]]$RT_ref, cfg$peaks[[2]]$RT_window)

  # Integrate TS and PD
  area_ts <- function(pk) {
    if (is.null(pk)) return(NA_real_)
    res <- tryCatch(
      integrate_peaks(sig_bc, RT_full, pk, baseline = NULL, methods = "TS"),
      error = function(e) NULL)
    if (is.null(res)) NA_real_ else res$area_TS
  }

  # Baseline value at the (compressed) Methane RT
  me_RT  <- if (!is.null(me)) me$RT else cfg$peaks[[1]]$RT_ref
  eth_RT <- if (!is.null(eth)) eth$RT else cfg$peaks[[2]]$RT_ref
  bl_at_me  <- approx(RT2, bl, xout = me_RT,  rule = 2)$y
  bl_at_eth <- approx(RT2, bl, xout = eth_RT, rule = 2)$y

  list(
    inh_end    = inh_end,
    n_peaks    = nrow(peaks),
    me_RT      = if (!is.null(me)) me$RT else NA_real_,
    me_height  = if (!is.null(me)) me$height else NA_real_,
    me_area    = area_ts(me),
    eth_RT     = if (!is.null(eth)) eth$RT else NA_real_,
    eth_height = if (!is.null(eth)) eth$height else NA_real_,
    eth_area   = area_ts(eth),
    bl_at_me   = bl_at_me,
    bl_at_eth  = bl_at_eth
  )
})

res_df <- do.call(rbind, lapply(results, as.data.frame))
cat("\nResults:\n")
print(res_df, digits = 4)

# ── Plot ──────────────────────────────────────────────────────────────────────
cat(sprintf("\nSaving: %s\n", out_png))
png(out_png, width = 1600, height = 1000, res = 110)

layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(4, 4.5, 2.5, 1), mgp = c(2.4, 0.7, 0))

x <- res_df$inh_end

# 1. Methane area vs Empower
plot(x, res_df$me_area, type = "b", pch = 16, col = "steelblue",
     xlab = "Inhibit end (min)", ylab = "Methane area (TS)",
     main = "Methane area_TS")
abline(h = EMP_METHANE, col = "tomato", lty = 2, lwd = 1.5)
legend("topleft", c("mcescher", "Empower"), col = c("steelblue","tomato"),
       lty = c(1,2), pch = c(16,NA), bty = "n", cex = 0.8)

# 2. Ethane area vs Empower
plot(x, res_df$eth_area, type = "b", pch = 16, col = "steelblue",
     xlab = "Inhibit end (min)", ylab = "Ethane area (TS)",
     main = "Ethane area_TS")
abline(h = EMP_ETHANE, col = "tomato", lty = 2, lwd = 1.5)

# 3. arPLS baseline at Methane RT
plot(x, res_df$bl_at_me, type = "b", pch = 16, col = "darkorange",
     xlab = "Inhibit end (min)", ylab = "arPLS value at Methane RT (comp.)",
     main = "Baseline at Methane RT")
abline(h = 0, col = "grey50", lty = 2)

# 4. arPLS baseline at Ethane RT
plot(x, res_df$bl_at_eth, type = "b", pch = 16, col = "darkorange",
     xlab = "Inhibit end (min)", ylab = "arPLS value at Ethane RT (comp.)",
     main = "Baseline at Ethane RT")
abline(h = 0, col = "grey50", lty = 2)

# 5. Methane height (after baseline)
plot(x, res_df$me_height, type = "b", pch = 16, col = "forestgreen",
     xlab = "Inhibit end (min)", ylab = "Methane height",
     main = "Methane peak height")
abline(h = 0, col = "grey50", lty = 2)

# 6. Number of detected peaks
plot(x, res_df$n_peaks, type = "b", pch = 16, col = "purple",
     xlab = "Inhibit end (min)", ylab = "n peaks detected",
     main = "Total peaks detected")

dev.off()
cat("Done.\n")
