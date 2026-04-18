#  Baseline correction using adaptive iteratively reweighted Penalized Least Squares;		
#  Input 
#         X:row matrix of spectra or chromatogram (size m*n, m is sample and n is variable)
#         lambda: lambda is an adjustable parameter, it can be adjusted by user. The larger lambda is, the smoother z will be 
#         order: an integer indicating the order of the difference of penalties
#         wep: weight exception proportion at both the start and end
#         p: asymmetry parameter for the start and end
#         itermax: maximum iteration times
#  Output
#         Xc: the corrected spectra or chromatogram vector (size m*n)
#         Z: the fitted vector (size m*n)
#  Examples:
#         Xc=airPLS(X);
#         [Xc,Z]=airPLS(X,10e5,2,0.1,0.5,20);
#  Reference:
#         (1) Eilers, P. H. C., A perfect smoother. Analytical Chemistry 75 (14), 3631 (2003).
#         (2) Eilers, P. H. C., Baseline Correction with Asymmetric Least
#         Squares Smoothing, http://www.science.uva.nl/~hboelens/publications/draftpub/Eilers_2005.pdf
#         (3) Gan, Feng, Ruan, Guihua, and Mo, Jinyuan, Baseline correction by improved iterative polynomial fitting with automatic threshold. Chemometrics and Intelligent Laboratory Systems 82 (1-2), 59 (2006).
# 
#  zhimin zhang @ central south university on Mar 30,2011


airPLS <- function(X, lambda=1e8, order=2, wep=0.1, p=0.05, itermax=20) {
  library(Matrix)
  
  if (missing(X)) stop("airPLS:NotEnoughInputs - Not enough input arguments. See airPLS.")
  
  m <- nrow(X)
  n <- ncol(X)
  
  wi <- c(seq_len(ceiling(n * wep)), floor(n - n * wep + 1):n)
  
  # Construct difference matrix D of given order
  D <- diff(diag(n), differences = order)
  DD <- lambda * crossprod(D)
  
  Z <- matrix(0, m, n)
  
  for (i in seq_len(m)) {
    w <- rep(1, n)
    x <- X[i, ]
    
    for (j in seq_len(itermax)) {
      W <- Diagonal(x = w)
      # Solve (W + DD) z = W * x
      C <- chol(W + DD)
      z <- solve(C, solve(t(C), w * x))
      d <- x - z
      dssn <- abs(sum(d[d < 0]))
      
      if (dssn < 0.001 * sum(abs(x))) {
        break
      }
      
      w[d >= 0] <- 0
      w[wi] <- p
      w[d < 0] <- exp(j * abs(d[d < 0]) / dssn)
    }
    Z[i, ] <- z
  }
  
  Xc <- X - Z
  list(Xc = Xc, Z = Z)
}
