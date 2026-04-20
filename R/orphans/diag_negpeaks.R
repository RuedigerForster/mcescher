source("./import/import_chromas.R")
source("./compress/compress.R")
source("./inhibit/inhibit.R")
source("./atsa/atsa.R")
source("./sass/sass.R")

sample_types <- factor(c("UNKNOWN","CONTROL","CALIB","BLANK"), ordered=TRUE,
                       levels=c("UNKNOWN","CONTROL","CALIB","BLANK"))
df  <- import_chromas("../data/", "^UNIS_FID1061008\\.cdf$")
RT  <- attr(df, "retention_time")
x_factor <- 20L
lt  <- nrow(df) - nrow(df) %% x_factor
df2 <- compress_mat(df[1:lt, , drop = FALSE], c(x_factor, 1L))
RT2 <- compress_vec(RT[1:lt], x_factor)
rm(df)

atsa_r <- atsa(df2, RT2, segment_size = 3, shift_max = 0.5,
               smooth_width = 11L, amp_thresh = 0,
               baseline_correct = FALSE, verbose = FALSE)
df2 <- atsa_r$aligned; rm(atsa_r)

inh <- inhibit(RT2, list(c(0, 1.1), c(20.38, 27)))
df3 <- matrix(NA_real_, nrow = nrow(df2), ncol = 1)
for (seg in inh$segments) {
  idx <- seg["begin"]:seg["end"]
  res <- sass_L1(as.numeric(df2[idx, 1]), 2, 0.022, 1, 1.2)
  df3[idx, 1] <- as.numeric(res$x)
}

corr <- df2[, 1] - df3[, 1]
win  <- which(RT2 >= 1.0 & RT2 <= 3.5)

cat(sprintf("%-7s  %8s  %10s  %10s\n", "RT", "raw", "baseline", "corrected"))
for (i in win)
  cat(sprintf("%-7.3f  %8.4f  %10.4f  %10.4f\n",
              RT2[i], df2[i, 1], df3[i, 1], corr[i]))
