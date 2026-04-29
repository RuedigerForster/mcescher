# airPLS_liu.R
# Pure-R port of airPLS by Shixuan Liu (2014).
# Reference Python: https://github.com/zmzhang/airPLS
#
# Key differences from the Zhang (2011) port in airPLS.R:
#   - First-order difference penalty (porder = 1)
#   - No endpoint anchoring (wep / p removed)
#   - Weights for above-baseline points zeroed; endpoint fix prevents
#     the first/last sample from being permanently zeroed too
#   - Weight update: exp(j * |d| / dssn) for below-baseline points only
#
# Interface identical to airPLS() for drop-in replacement.
# Input:  X  m  n row matrix (one row per chromatogram)
# Output: list(Xc = corrected matrix, Z = baseline matrix)

.airPLS_liu <- function(X, lambda = 1e7, porder = 1L, itermax = 15L) {
  m <- nrow(X)
  n <- ncol(X)

  D  <- diff(Matrix::Diagonal(n), differences = porder)
  DD <- lambda * crossprod(D)

  Z <- matrix(0.0, m, n)

  for (i in seq_len(m)) {
    x <- X[i, ]
    w <- rep(1.0, n)

    z <- x
    for (j in seq_len(itermax)) {
      # Small floor keeps W + DD positive definite when weights hit zero
      W    <- Matrix::Diagonal(x = pmax(w, 1e-10))
      z    <- as.numeric(solve(W + DD, w * x))
      d    <- x - z
      dssn <- sum(abs(d[d < 0]))

      if (dssn < 0.001 * sum(abs(x))) break

      neg      <- d < 0
      w[]      <- 0.0
      w[neg]   <- exp(j * abs(d[neg]) / dssn)
      # Endpoint fix: anchors first/last point so the baseline cannot
      # diverge at the boundaries (Liu 2014, lines 19-20)
      w_end  <- exp(j * max(d[neg]) / dssn)
      w[1L]  <- w_end
      w[n]   <- w_end
    }
    Z[i, ] <- z
  }

  list(Xc = X - Z, Z = Z)
}
