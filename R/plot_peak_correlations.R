# plot_peak_correlations.R
# Correlation matrix of named-peak areas across chromatograms:
# area_TS (Tangent Skim) vs area_cons (baseline ensemble consensus).
#
# Usage (from R/ directory):
#   Rscript plot_peak_correlations.R [peaks_csv]
# Default: ../data/mcescher_peaks.csv

setwd("/home/nvidiauser/Developments/mcescher/R")

args      <- commandArgs(trailingOnly = TRUE)
peaks_csv <- if (length(args) >= 1) args[1] else "../data/mcescher_peaks.csv"
out_png   <- sub("\\.csv$", "_correlations.png", peaks_csv)

if (!file.exists(peaks_csv)) stop("CSV not found: ", peaks_csv)

NAMED <- c("Methane", "Ethane", "Propane", "iso-Butane", "n-Butane")

pk <- read.csv(peaks_csv, stringsAsFactors = FALSE)
pk <- pk[pk$name %in% NAMED, ]
pk$file_base <- basename(pk$chrom)

# ── Wide matrices: rows = chromatograms, cols = peaks ─────────────────────────
to_wide <- function(df, area_col) {
  mat <- reshape(df[, c("file_base", "name", area_col)],
                 idvar    = "file_base",
                 timevar  = "name",
                 direction = "wide")
  rownames(mat) <- mat$file_base
  mat <- mat[, paste0(area_col, ".", NAMED), drop = FALSE]
  colnames(mat) <- NAMED
  as.matrix(mat)
}

mat_ts   <- to_wide(pk, "area_TS")
mat_cons <- to_wide(pk, "area_cons")

# ── Correlation matrices ──────────────────────────────────────────────────────
cor_ts   <- cor(mat_ts,   use = "pairwise.complete.obs")
cor_cons <- cor(mat_cons, use = "pairwise.complete.obs")

# ── Plot helpers ──────────────────────────────────────────────────────────────
heat_col <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)

plot_corrmat <- function(mat, title) {
  n   <- ncol(mat)
  nms <- colnames(mat)

  # colour limits -1 to +1
  image(1:n, 1:n, t(mat[n:1, ]),
        col  = heat_col, zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "", main = title)
  box()

  axis(1, at = 1:n, labels = nms, las = 2, cex.axis = 0.82, tick = FALSE)
  axis(2, at = 1:n, labels = rev(nms), las = 1, cex.axis = 0.82, tick = FALSE)

  for (i in seq_len(n))
    for (j in seq_len(n)) {
      val <- mat[n + 1 - j, i]
      col <- if (abs(val) > 0.6) "white" else "black"
      text(i, j, sprintf("%.2f", val), cex = 0.78, col = col)
    }
}

# ── Scatter panels: TS vs cons per peak ───────────────────────────────────────
plot_scatter <- function(mat_x, mat_y, peak, xlab, ylab) {
  x <- mat_x[, peak]; y <- mat_y[, peak]
  ok <- complete.cases(x, y)
  lims <- range(c(x[ok], y[ok]), na.rm = TRUE)
  plot(x[ok], y[ok], pch = 16, cex = 0.8,
       col = adjustcolor("steelblue", 0.8),
       xlim = lims, ylim = lims,
       xlab = xlab, ylab = ylab,
       main = peak)
  abline(0, 1, col = "grey50", lty = 2)
  r <- cor(x[ok], y[ok])
  legend("topleft", legend = sprintf("r = %.3f", r),
         bty = "n", cex = 0.8)
}

# ── Assemble plot ─────────────────────────────────────────────────────────────
# Layout: row 1 = two correlation matrices (4 cols each) + colourbar (2 cols)
#         row 2 = five scatter plots (2 cols each) — 10 cols total
# 8 panels: 1=TS corr, 2=cons corr, 3=colourbar, 4-8=scatter per peak
cat(sprintf("Saving: %s\n", out_png))
png(out_png, width = 1800, height = 1100, res = 100)

layout(rbind(c(1,1,1,1, 2,2,2,2, 3,3),
             c(4,4, 5,5, 6,6, 7,7, 8,8)),
       heights = c(1.4, 1))
par(mar = c(4.5, 4.5, 2.5, 1))

plot_corrmat(cor_ts,   "area_TS  (Tangent Skim)")
plot_corrmat(cor_cons, "area_cons  (Baseline Ensemble)")

# Colour bar
par(mar = c(5, 0.5, 3, 2.5))
bar_y <- seq(-1, 1, length.out = 101)
image(1, bar_y, matrix(bar_y, nrow = 1),
      col = heat_col, axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(-1, -0.5, 0, 0.5, 1), las = 1, cex.axis = 0.8)
box()

# Scatter row
par(mar = c(4, 4, 2.5, 0.5))
for (pk_nm in NAMED)
  plot_scatter(mat_ts, mat_cons, pk_nm, "area_TS", "area_cons")

dev.off()
cat("Done.\n")
