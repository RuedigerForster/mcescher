## Description
# DOARPLS Estimate baseline with arPLS in MATLAB
# Function published in :
# Baek, S.-J., Park, A., Ahn, Y.-J., Choo, J. (2015)
# "Baseline correction using asymmetrically reweighted penalized least
# squares smoothing" Analyst, 140, 250-257. doi: 10.1039/c4an01061b
#
##

# input: y y values,
#        z baseline model
# output:
#       bslPts false if peak pts, true in baseline pts


#' Baseline correction using arPLS (asymmetrically reweighted penalized least squares).
#'
#' @param y      Numeric signal vector.
#' @param lambda Smoothness penalty (larger = smoother). Typical range 1e51e9.
#' @param ratio  Convergence criterion on weight change. Default 1e-6 works well.
#' @return List: \code{$z} baseline vector, \code{$bslPts} logical baseline mask.
#' @export
doArPLS <- function(y, lambda, ratio) {
  N <- length(y)
  # Sparse difference matrix  avoids materialising a dense NN identity matrix
  D <- diff(Matrix::Diagonal(N), differences = 2)
  H <- lambda * crossprod(D)
  w <- rep(1, N)
  iterMax <- 20

  for (ii in 1:iterMax) {
    W <- Matrix::Diagonal(x = w)          # sparse diagonal  not diag(w) (dense NN)
    z <- as.numeric(solve(W + H, w * y))
    d <- y - z
    dn <- d[d < 0]
    m <- mean(dn)
    s <- sd(dn)
    overflow <- log(.Machine$double.xmax) / 2
    arg <- 2 * (abs(d) - (2 * s - m)) / s
    wt <- ifelse(arg >= overflow, 0, 1 / (1 + exp(arg)))

    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < ratio) break
    w <- wt
  }

  bslPts <- w > 0.05
  list(z = z, bslPts = bslPts)
}
