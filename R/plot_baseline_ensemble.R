# plot_baseline_ensemble.R
# Plot all ensemble baseline fits overlaid on one chromatogram.
#
# Usage (from R/ directory):
#   Rscript plot_baseline_ensemble.R [data_dir] [cdf_basename]
# Defaults to one sample from UNIS_FID_9000002_908779.

setwd("/home/nvidiauser/Developments/mcescher/R")

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "../data/UNIS_FID_9000002_908779/"
cdf_file <- if (length(args) >= 2) args[2] else "UNIS_FID_9000002_909527.cdf"
out_png  <- sub("\\.cdf$", "_baseline_ensemble.png", cdf_file, ignore.case = TRUE)

if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
cdf_path <- paste0(data_dir, cdf_file)
if (!file.exists(cdf_path)) stop("CDF not found: ", cdf_path)
options(mc.cores = parallel::detectCores(logical = FALSE))

suppressPackageStartupMessages({
  source("./config/read_config.R")
  source("./import/import_chromas.R")
  source("./compress/compress.R")
  source("./inhibit/inhibit.R")
  source("./sass/sass.R")
  source("./arpls/doArPLS.R")
  source("./integrate/baseline_ensemble.R")
  source("./atsa/peak_detect.R")
})

cfg             <- read_method_config("./method.yaml")
x_factor        <- cfg$processing$x_factor
inhibit_regions <- cfg$inhibit

# ── Load single CDF ──────────────────────────────────────────────────────────
cat(sprintf("Loading: %s\n", cdf_file))
df  <- import_chromas(data_dir, paste0("^", cdf_file, "$"))
RT  <- attr(df, "retention_time")
lt  <- nrow(df) - nrow(df) %% x_factor

# ── Compress ─────────────────────────────────────────────────────────────────
df2  <- compress_mat(df[1:lt, , drop = FALSE], c(x_factor, 1L))
RT2  <- compress_vec(RT[1:lt], x_factor)
rm(df); gc()

inh <- inhibit(RT2, inhibit_regions = inhibit_regions)

# ── SASS denoising ────────────────────────────────────────────────────────────
cat("SASS denoising...\n")
df2_dn <- df2
for (seg in inh$segments) {
  idx <- seg["begin"]:seg["end"]
  df2_dn[idx, 1] <- as.numeric(
    sass_L1(as.numeric(df2[idx, 1]),
            cfg$denoising$d, cfg$denoising$fc,
            cfg$denoising$K, cfg$denoising$lam)$x)
}
sig_sass <- as.numeric(df2_dn[, 1])   # compressed SASS-denoised signal

# ── arPLS primary baseline ────────────────────────────────────────────────────
cat("arPLS baseline...\n")
bl_arpls <- rep(NA_real_, length(sig_sass))
for (seg in inh$segments) {
  idx <- seg["begin"]:seg["end"]
  bl_arpls[idx] <- doArPLS(sig_sass[idx],
                            lambda = cfg$baseline$lambda,
                            ratio  = cfg$baseline$ratio)$z
}
bl_arpls[inh$block] <- sig_sass[inh$block]

# Baseline-corrected compressed signal for ensemble
sig_bc <- pmax(sig_sass - bl_arpls, 0)

# ── Decompress to full resolution for peak detection + integration ────────────
RT_full       <- RT[1:lt]
sig_sass_full <- as.numeric(decompress_mat(matrix(sig_sass,  ncol = 1L), lt))
bl_arpls_full <- as.numeric(decompress_mat(matrix(bl_arpls,  ncol = 1L), lt))
sig_bc_full   <- pmax(sig_sass_full - bl_arpls_full, 0)
inh_full      <- inhibit(RT_full, inhibit_regions = inhibit_regions)
sig_bc_full[inh_full$block] <- 0

# ── Detect peaks at full resolution (matches production pipeline) ─────────────
source("./atsa/peak_detect.R")
peaks <- find_peaks(sig_bc_full, RT_full,
                    smooth_width = cfg$detection$smooth_width_factor * x_factor,
                    amp_thresh   = cfg$detection$amp_thresh)
cat(sprintf("  %d peaks detected\n", nrow(peaks)))

# ── Baseline ensemble on compressed signal, integrate at full resolution ──────
cat("Running baseline ensemble (14 configs)...\n")
ens <- baseline_ensemble(sig_bc, RT2, peaks,
                         integrate_method = cfg$integration$ensemble_method,
                         configs          = default_ensemble_configs(),
                         verbose          = TRUE,
                         signal_full      = sig_bc_full,
                         RT_full          = RT_full)
cat("Done.\n")

# ── Config metadata ───────────────────────────────────────────────────────────
cfg_names <- ens$config_names
n_cfg     <- length(cfg_names)

fam <- sub("_liu", "", sub("_lam.*|_span.*|_p[0-9].*|_fc.*", "", cfg_names))
fam_pal <- c(airPLS  = "#e41a1c",
             arPLS   = "#ff7f00",
             LMV     = "#4daf4a",
             flatfit = "#984ea3",
             BEADS   = "#a65628")
cfg_col <- fam_pal[fam]

# Total baseline for each config (on compressed signal):
#   total_k = arPLS + ensemble_residual_k
bl_total <- sweep(ens$baselines, 1, bl_arpls, "+")   # n_comp × n_cfg
bl_total[inh$block, ] <- sig_sass[inh$block]          # inhibit → no correction

# Baseline-corrected signal per config
sig_corr <- sweep(-bl_total, 1, sig_sass, "+")        # sig_sass − total_k
sig_corr[sig_corr < 0] <- 0

# ── Named peak windows from method config ────────────────────────────────────
named_peaks <- cfg$peaks
peak_windows <- lapply(named_peaks, function(p)
  c(p$RT_ref - p$RT_window * 1.8,
    p$RT_ref + p$RT_window * 1.8))

# ── Plot ─────────────────────────────────────────────────────────────────────
cat(sprintf("Saving plot: %s\n", out_png))

# Layout: 1 wide top panel + 3 zoom panels below
png(out_png, width = 1600, height = 900, res = 110)

layout(matrix(c(1, 1, 1,
                2, 3, 4), nrow = 2, byrow = TRUE),
       heights = c(1.6, 1))
par(mar = c(3.5, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))

# ── Panel 1: full chromatogram ────────────────────────────────────────────────
rt_range <- c(0.8, 20)
mask     <- RT2 >= rt_range[1] & RT2 <= rt_range[2]

y_lim  <- c(-30, 100)

plot(RT2[mask], sig_sass[mask], type = "l", col = "steelblue", lwd = 1.5,
     xlim = rt_range, ylim = y_lim,
     xlab = "Retention time (min)", ylab = "Signal",
     main = sprintf("Ensemble baselines — %s", cdf_file))

# Draw ensemble baselines (thin, semi-transparent by family)
for (k in seq_len(n_cfg)) {
  lines(RT2[mask], bl_total[mask, k], col = adjustcolor(cfg_col[k], 0.7), lwd = 1)
}
# Primary arPLS on top
lines(RT2[mask], bl_arpls[mask], col = "black", lwd = 2, lty = 2)
abline(h = 0, col = "grey70", lty = 3)

# Peak RT markers
pk_RTs <- vapply(named_peaks, `[[`, 0.0, "RT_ref")
abline(v = pk_RTs, col = "grey50", lty = 3)
axis(3, at = pk_RTs,
     labels = vapply(named_peaks, `[[`, "", "name"),
     cex.axis = 0.7, las = 2, tick = FALSE, line = -0.5)

# Legend
uniq_fam <- unique(fam)
legend("topright",
       legend = c("Signal (SASS)", "arPLS (primary)", uniq_fam),
       col    = c("steelblue", "black", fam_pal[uniq_fam]),
       lty    = c(1, 2, rep(1, length(uniq_fam))),
       lwd    = c(1.5, 2, rep(1.5, length(uniq_fam))),
       cex    = 0.75, bg = "white")

# ── Panels 2-4: peak-region zoom ─────────────────────────────────────────────
zoom_regions <- list(
  "Methane / Ethane" = c(0.9, 2.5),
  "Propane"          = c(4.5, 7.0),
  "iso-Butane / n-Butane" = c(14.0, 17.5)
)

for (reg_name in names(zoom_regions)) {
  rng  <- zoom_regions[[reg_name]]
  mask <- RT2 >= rng[1] & RT2 <= rng[2]

  y_corr <- sig_corr[mask, ]
  y_lim  <- c(0, quantile(y_corr, 0.999, na.rm = TRUE))

  plot(NA, xlim = rng, ylim = y_lim,
       xlab = "RT (min)", ylab = "Corrected signal",
       main = reg_name)

  for (k in seq_len(n_cfg)) {
    lines(RT2[mask], sig_corr[mask, k],
          col = adjustcolor(cfg_col[k], 0.75), lwd = 1)
  }
  # Weighted ensemble mean in black
  if (!is.null(ens$weights)) {
    mean_corr <- as.numeric(sig_corr[mask, ] %*% ens$weights)
    lines(RT2[mask], mean_corr, col = "black", lwd = 2)
  }
  abline(h = 0, col = "grey70", lty = 3)

  # Mark named peak RTs that fall in this window
  in_window <- which(pk_RTs >= rng[1] & pk_RTs <= rng[2])
  if (length(in_window)) {
    abline(v = pk_RTs[in_window], col = "grey40", lty = 3)
    axis(3, at = pk_RTs[in_window],
         labels = vapply(named_peaks[in_window], `[[`, "", "name"),
         cex.axis = 0.65, las = 2, tick = FALSE, line = -0.5)
  }
}

dev.off()
cat("Plot saved:", out_png, "\n")
