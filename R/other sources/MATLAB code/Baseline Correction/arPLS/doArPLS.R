## Description
# DOARPLS Estimate baseline with arPLS in MATLAB
# Function published in :
# Baek, S.-J., Park, A., Ahn, Y.-J., Choo, J. (2015) 
# "Baseline correction using asymmetrically reweighted penalized least 
# squares smoothing" Analyst, 140, 250-257. doi: 10.1039/c4an01061b
#
##

# input: y y values, 
#        z basline model
# output:
#       bslPts false if peak pts, true in baseline pts


doArPLS <- function(y, lambda, ratio) {
  N <- length(y)
  D <- diff(diag(N), differences = 2)
  H <- lambda * t(D) %*% D
  w <- rep(1, N)
  iterMax <- 20
  
  for (ii in 1:iterMax) {
    W <- diag(w)
    # Cholesky decomposition
    C <- chol(W + H)
    z <- backsolve(C, forwardsolve(t(C), w * y))
    d <- y - z
    dn <- d[d < 0]
    m <- mean(dn)
    s <- sd(dn)
    wt <- 1 / (1 + exp(2 * (abs(d) - (2 * s - m)) / s))
    
    if (sqrt(sum((w - wt)^2)) / sqrt(sum(w^2)) < ratio) break
    w <- wt
  }
  
  bslPts <- w > 0.05
  list(z = z, bslPts = bslPts)
}
