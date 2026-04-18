# flatfit.R
# Pure-R port of the FlatFit baseline algorithm from the MOCCA package (Bayer).
# Reference: github.com/bayer-science-for-a-better-life/mocca
#
# Finds baseline z that minimises:
#   (y - z)' W (y - z)  +  smoothness * L^4 * z' D'D z
#
# where W = diag(w) and weights  w = 1 / (slope^2 + curvature^2 + eps).
# High-slope / high-curvature regions (peak flanks) receive low weight, so
# the baseline is pulled away from peaks without iterative sign logic.
#
# The L^4 scaling of the smoothness penalty makes the algorithm invariant to
# signal length / sampling frequency (matches the original Python).
#
# Parameters:
#   y          : numeric vector (signal segment, baseline-free regions only)
#   smoothness : smoothness penalty scale (default 1e-3; try 1e-4 to 1e-2)
#   p          : Savitzky-Golay window as fraction of signal length (default 0.05)
#
# Returns:
#   Numeric vector of baseline estimates, same length as y.

library(Matrix)
library(pracma)   # for savgol()

flatfit_baseline <- function(y, smoothness = 1e-3, p = 0.05) {
  y <- as.numeric(y)
  L <- length(y)
  if (L < 5L) return(rep(mean(y), L))   # degenerate segment

  # ---- Savitzky-Golay window: odd, >= 5, <= L-2 ----------------------------
  fl <- max(5L, as.integer(L * p))
  if (fl %% 2L == 0L) fl <- fl + 1L       # pracma::savgol requires odd
  fl <- min(fl, L - 2L)
  if (fl %% 2L == 0L) fl <- fl - 1L       # re-enforce odd after min()
  fl <- max(fl, 5L)                        # safety floor

  # ---- First derivative (slope) via Savitzky-Golay -------------------------
  slope <- pracma::savgol(y, fl = fl, forder = 3L, dorder = 1L)

  # ---- Second derivative (curvature) via central differences of slope ------
  curvature        <- numeric(L)
  curvature[1L]    <- slope[2L]   - slope[1L]
  curvature[L]     <- slope[L]    - slope[L - 1L]
  if (L > 2L)
    curvature[2L:(L-1L)] <- (slope[3L:L] - slope[1L:(L-2L)]) / 2

  # ---- Normalised inverse-derivative weights --------------------------------
  s2 <- slope^2;     s2 <- s2 / (sum(s2) + .Machine$double.eps)
  c2 <- curvature^2; c2 <- c2 / (sum(c2) + .Machine$double.eps)
  w  <- 1 / (s2 + c2 + 1e-10)

  # ---- Second-difference smoothness penalty: H = smoothness * L^4 * D'D ----
  D  <- diff(Diagonal(L), differences = 2L)   # (L-2) x L sparse matrix
  H  <- (smoothness * L^4) * crossprod(D)     # L x L

  # ---- Single weighted least-squares solve ----------------------------------
  W  <- Diagonal(x = w)
  z  <- as.numeric(solve(W + H, w * y))

  z
}
