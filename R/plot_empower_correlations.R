# plot_empower_correlations.R
# Correlation matrix of named-peak areas: Empower CDS vs mcescher ensemble.
#
# Usage (from R/ directory):
#   Rscript plot_empower_correlations.R [comparison_csv]
# Default: ../data/comparison_table.csv

setwd("/home/nvidiauser/Developments/mcescher/R")

args   <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[1] else "../data/comparison_table.csv"
out_png <- sub("\\.csv$", "_correlations.png", in_csv)

if (!file.exists(in_csv)) stop("CSV not found: ", in_csv)

NAMED <- c("Methane", "Ethane", "Propane", "iso-Butane", "n-Butane")

ct <- read.csv(in_csv, stringsAsFactors = FALSE)
ct <- ct[ct$name %in% NAMED, ]

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

mat_emp  <- to_wide(ct, "empower_area")
mat_cons <- to_wide(ct, "area_cons")

# ── Correlation matrices ──────────────────────────────────────────────────────
cor_emp  <- cor(mat_emp,  use = "pairwise.complete.obs")
cor_cons <- cor(mat_cons, use = "pairwise.complete.obs")

# ── Plot helpers ──────────────────────────────────────────────────────────────
heat_col <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)

plot_corrmat <- function(mat, title) {
  n   <- ncol(mat)
  nms <- colnames(mat)
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

# ── Scatter panels: Empower vs cons per peak ──────────────────────────────────
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
  r   <- cor(x[ok], y[ok])
  lm1 <- lm(y[ok] ~ x[ok])
  abline(lm1, col = "tomato", lty = 1, lwd = 1.2)
  legend("topleft",
         legend = c(sprintf("r = %.3f", r),
                    sprintf("slope = %.3f", coef(lm1)[2])),
         bty = "n", cex = 0.75)
}

# ── Assemble plot ─────────────────────────────────────────────────────────────
cat(sprintf("Saving: %s\n", out_png))
png(out_png, width = 1800, height = 1100, res = 100)

layout(rbind(c(1,1,1,1, 2,2,2,2, 3,3),
             c(4,4, 5,5, 6,6, 7,7, 8,8)),
       heights = c(1.4, 1))
par(mar = c(4.5, 4.5, 2.5, 1))

plot_corrmat(cor_emp,  "Empower  (CDS peak table)")
plot_corrmat(cor_cons, "area_cons  (Baseline Ensemble)")

# Colour bar
par(mar = c(5, 0.5, 3, 2.5))
bar_y <- seq(-1, 1, length.out = 101)
image(1, bar_y, matrix(bar_y, nrow = 1),
      col = heat_col, axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(-1, -0.5, 0, 0.5, 1), las = 1, cex.axis = 0.8)
box()

# Scatter row: Empower (x) vs area_cons (y)
par(mar = c(4, 4, 2.5, 0.5))
for (pk_nm in NAMED)
  plot_scatter(mat_emp, mat_cons, pk_nm, "Empower area", "area_cons")

dev.off()
cat("Done.\n")
