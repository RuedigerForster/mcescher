source("./import/import_chromas.R")
source("./compress/compress.R")
source("./inhibit/inhibit.R")
source("./atsa/atsa.R")    # Zheng et al. 2017, Sci.Rep. 7:256
source("./sass/sass.R")
source("./summarise/summarise.R")

library(lubridate)

tm <- Sys.time()
sample_types <- factor(c("UNKNOWN", "CONTROL", "CALIB", "BLANK"), ordered = TRUE,
                       levels = c("UNKNOWN", "CONTROL", "CALIB", "BLANK"))

safe <- par(no.readonly = TRUE)

# --- load data ----------------------------------------------------------------
path    <- "../data/"
pattern <- "^UNIS_FID.*\\.cdf$"
df <- import_chromas(path, pattern)
# operator did not set sample_type correctly; treat all as BLANK baselines
attr(df, "sample_type") <- rep(as.character(which(sample_types == "BLANK")), ncol(df))

RT <- attr(df, "retention_time")   # minutes

# --- select blank / baseline runs in chronological order ---------------------
blanks     <- as.integer(attr(df, "sample_type")) == which(sample_types == "BLANK")
blanks_idx <- which(blanks)
idx        <- blanks_idx[order(lubridate::as_datetime(attr(df, "run_datetime"))[blanks_idx])]
df2    <- df[, idx, drop = FALSE]

# --- downsample ---------------------------------------------------------------
x_factor <- 20L; y_factor <- 1L
lt  <- nrow(df2) - nrow(df2) %% x_factor
df2 <- compress_mat(df2[1:lt, ], c(x_factor, y_factor))
RT2 <- compress_vec(RT[1:lt], x_factor)   # keep RT aligned with df2

# --- ATSA time-shift alignment -----------------------------------------------
# Aligns all chromatograms to a reference before baseline fitting.
# Set baseline_correct=FALSE here since we apply SASS separately below.
# Increase amp_thresh if too many noise peaks are detected.
atsa_result <- atsa(df2, RT2,
                    segment_size     = 3.0,
                    shift_max        = 0.5,
                    smooth_width     = 11L,
                    amp_thresh       = 0,
                    baseline_correct = FALSE,
                    verbose          = TRUE)
df2 <- atsa_result$aligned

# --- inhibit regions (minutes) -----------------------------------------------
# Adjust these windows to cover valve-switching artifacts visible in your data.
# List of c(t_start, t_end) pairs; leave empty list() if none needed.
inhibit_regions <- list(
  # c(0.0, 1.0)    # example: dead-time / solvent front
  # c(12.0, 12.5)  # example: valve switch artifact
)
inh <- inhibit(RT2, inhibit_regions = inhibit_regions)

# --- SASS smoothing per good segment -----------------------------------------
d <- 2; fc <- 0.022   # low-pass filter parameters
K <- 1; lam <- 1.2    # SASS sparsity parameters

df3 <- matrix(NA_real_, nrow = nrow(df2), ncol = ncol(df2))
for (seg in inh$segments) {
  idx_seg <- seg["begin"]:seg["end"]
  for (col in seq_len(ncol(df2))) {
    y <- as.numeric(df2[idx_seg, col])
    res <- sass_L1(y, d, fc, K, lam)
    df3[idx_seg, col] <- as.numeric(res$x)
  }
}

# --- RMSE of fitted baselines (LOESS smoothed) --------------------------------
index <- seq_len(nrow(df3))
RMSE  <- sapply(seq_len(ncol(df3)), function(i) sqrt((df3[, i] - df2[, i])^2))

RMSE_fit  <- lapply(seq_len(ncol(RMSE)), function(i)
  loess(RMSE[, i] ~ index, span = 0.05, na.action = na.omit))

RMSE_pred <- matrix(NA_real_, nrow = nrow(RMSE), ncol = ncol(RMSE))
good_idx  <- !inh$block
RMSE_pred[good_idx, ] <- sapply(seq_len(ncol(df3)), function(i)
  predict(RMSE_fit[[i]]))

# --- example plot (column 1) --------------------------------------------------
col_ex <- 1
areaplot::confplot(
  cbind(index,
        df3[, col_ex] - 5 * RMSE_pred[, col_ex],
        df3[, col_ex] + 5 * RMSE_pred[, col_ex]),
  ylim = range(df2[, col_ex], na.rm = TRUE))
lines(df3[, col_ex], type = "l")

# --- Decompress to full resolution for integration ---------------------------
# Integration on block-averaged data is not comparable to full-resolution methods.
df2_full <- decompress_mat(df2, lt)
df3_tmp  <- df3; df3_tmp[is.na(df3_tmp)] <- 0
df3_full <- decompress_mat(df3_tmp, lt)
RT_full  <- RT[1:lt]
inh_full <- inhibit(RT_full, inhibit_regions = inhibit_regions)

# --- pipeline summary ---------------------------------------------------------
# Set method_lod / method_loq to your validated method limits (signal units).
# Set run_ensemble = FALSE to skip the slow baseline uncertainty step.
result <- pipeline_summary(
  signal_mat   = df2_full,
  baseline_mat = df3_full,
  RT           = RT_full,
  atsa_result  = NULL,             # shifts are in compressed-sample units; skip display
  inh          = inh_full,
  RMSE_pred    = NULL,             # computed on compressed data; skip display
  col_names    = attr(df, "source_file")[idx],
  method_lod   = NULL,   # e.g. 5.0
  method_loq   = NULL,   # e.g. 15.0
  cal_slope    = NULL,
  smooth_width = 11L * x_factor,  # scale to preserve smoothing window in minutes
  run_ensemble = TRUE,
  verbose      = TRUE
)

print(paste("Duration:", difftime(Sys.time(), tm)))
