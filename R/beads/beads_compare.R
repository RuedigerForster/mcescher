source("./beads/beads.R")

y <- as.numeric(read.table(
  "C:/Users/WinUser/Developments/Integration/src/BEADS-baseline-R-code/Data_1D_Baseline_Deg2.txt",
  header = FALSE, sep = "\t")[, 1])

fc_vals  <- c(0.005, 0.01, 0.05)
lam_sets <- list(
  old = c(lam0 = 0.4, lam1 = 4.0, lam2 = 3.2),
  new = c(lam0 = 0.5, lam1 = 5.0, lam2 = 4.0)
)

png("beads_compare.png", width = 1400, height = 1400)
par(mfrow = c(6, 1), mar = c(2, 4, 2, 1))

for (lam_name in names(lam_sets)) {
  lam <- lam_sets[[lam_name]]
  for (fc in fc_vals) {
    res <- beads(y, d = 1, fc = fc, r = 6,
                 lam0 = lam["lam0"], lam1 = lam["lam1"], lam2 = lam["lam2"])
    title <- sprintf("fc=%.3f  lam0=%.1f lam1=%.1f lam2=%.1f  (%s)",
                     fc, lam["lam0"], lam["lam1"], lam["lam2"], lam_name)
    plot(y, type = "l", col = "grey70", lwd = 1, main = title, ylab = "signal", xlab = "")
    lines(res$f, col = "firebrick", lwd = 2)
  }
}

dev.off()
message("Saved: R/beads_compare.png")
