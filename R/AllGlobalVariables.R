# Gaussian peak area coefficient: area = GAUSS_AREA_FACTOR * height * FWHM.
# Derived as sqrt(2*pi) / 2.3548 (2.3548 = 2*sqrt(2*log(2))).
GAUSS_AREA_FACTOR <- sqrt(2 * pi) / (2 * sqrt(2 * log(2)))

# Default colour palette for plot.ChromResult.
.chrom_plot_palette <- list(
  signal   = "steelblue",
  baseline = "firebrick",
  band     = "#CC000033",
  peak_ref = "orange"
)
