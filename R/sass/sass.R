library(Matrix)

# INPUT
#   d  : degree of filter is 2d
#   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
#   N  : length of signal
#   K  : order of difference matrix D (need 1 <= K <= 2*d) (default K = 1)

ABfilt <- function(deg, fc, N, K = 1) {
  if (K > 2 * deg) {
    stop("ABfilt: K > 2*d")
  }
  
  omc <- 2 * pi * fc
  t <- ((1 - cos(omc)) / (1 + cos(omc)))^deg
  
  # Define p such that P(z)P(1/z) = B(z), i.e., P'*P = B
  p <- 1
  for (k in 1:deg) {
    p <- convolve(p, c(-1, 1), type = "open")
  }
  
  # Construct banded matrix P
  rows_P <- N - deg
  cols_P <- N
  P <- bandSparse(rows_P, cols_P, k = 0:deg, diagonals = lapply(0:deg, function(i) rep(p[i + 1], rows_P)))
  
  B <- t(P) %*% P
  
  q <- sqrt(t)
  for (i in 1:deg) {
    q <- convolve(q, c(1, 1), type = "open")
  }
  
  Q <- bandSparse(rows_P, cols_P, k = 0:deg, diagonals = lapply(0:deg, function(i) rep(q[i + 1], rows_P)))
  
  A <- t(P) %*% P + t(Q) %*% Q
  
  if (K <= deg) {
    d <- 1
    for (i in 1:K) {
      d <- convolve(d, c(-1, 1), type = "open")
    }
    rows_D <- N - K
    cols_D <- N
    D <- bandSparse(rows_D, cols_D, k = 0:K, diagonals = lapply(0:K, function(i) rep(d[i + 1], rows_D)))
    
    # Polynomial deconvolution p1 = p / d
    p1 <- pracma::deconv(p, d)$q
    
    rows_P1 <- N - deg
    cols_P1 <- N - K
    P1 <- bandSparse(rows_P1, cols_P1, k = 0:(deg - K), diagonals = lapply(0:(deg - K), function(i) rep(p1[i + 1], rows_P1)))
    
    B1 <- t(P) %*% P1
    
    b1 <- convolve(p1, rev(p), type = "open")
  } else {
    K2 <- 2 * deg - K
    d <- 1
    for (i in 1:K2) {
      d <- convolve(d, c(-1, 1), type = "open")
    }
    rows_B1 <- N - K2
    cols_B1 <- N
    B1 <- t(bandSparse(rows_B1, cols_B1, k = 0:K2, diagonals = lapply(0:K2, function(i) rep(d[i + 1], rows_B1))))
    
    p1 <- pracma::deconv(p, d)$q
    
    rows_D1 <- N - deg
    cols_D1 <- N - K2
    D1 <- bandSparse(rows_D1, cols_D1, k = 0:(deg - K2), diagonals = lapply(0:(deg - K2), function(i) rep(p1[i + 1], rows_D1)))
    
    D <- t(D1) %*% P
    
    b1 <- d
  }
  
  a <- convolve(p, rev(p), type = "open") + convolve(q, rev(q), type = "open")
  b <- convolve(p, rev(p), type = "open")
  
  # verify that B = B1*D
  err <- B - B1 %*% D
  mae <- max(abs(err@x))
  if (mae > 1e-10) {
    message("Error in ABfilt (B1*D not equal to B)")
  }
  
  # Calculate filter norms
  imp <- numeric(dim(B1)[2])
  imp[round(N / 2)] <- 1
  
  h1 <- solve(A, B1 %*% imp)
  H1norm <- sqrt(sum(abs(h1)^2))
  
  hh <- t(B) %*% solve(A %*% t(A), B1 %*% imp)
  HTH1norm <- sqrt(sum(abs(hh)^2))
  
  list(A = A, B = B, B1 = B1, D = D, a = a, b = b, b1 = b1, H1norm = H1norm, HTH1norm = HTH1norm)
}


sass_L1 <- function(y, d, fc, K, lam, u_init = NULL) {
  # Sparsity-assisted signal smoothing (SASS)
  # using L1 norm as sparsity penalty.
  #
  # Signal model: y = f + g + noise
  #   f : low-pass signal
  #   g : sparse order-K derivative
  #
  # INPUT
  #   y - noisy data
  #   d - filter degree parameter (d = 1, 2, 3)
  #   fc - cut-off frequency (normalized, 0 < fc < 0.5)
  #   K - order of sparse derivative (1 <= K <= 2d)
  #   lam - regularization parameter
  #
  # OUTPUT
  #   x - output of SASS algorithm
  #   u - sparse signal
  #   cost - cost function history
  #   v - filtered sparse signal, A\(B1*u)
  
  MAX_ITER <- 10000
  TOL_STOP <- 1e-4
  y <- as.matrix(y)
  if (ncol(y) > 1) y <- t(y) # ensure column vector
  y <- as.numeric(y)
  N <- length(y)
  cost <- numeric(MAX_ITER)
  
  # ABfilt function needs to be defined or imported
  # It should return list(A, B, B1, D)
  filt <- ABfilt(d, fc, N, K)
  A <- filt$A
  B <- filt$B
  B1 <- filt$B1
  D <- filt$D
  
  # Define H function: H(x) = solve(A, B %*% x)
  H <- function(x) {
    solve(A, B %*% x)
  }
  
  AAT <- A %*% t(A)
  Hy <- H(y)
  b <- t(B1) %*% solve(AAT, B %*% y)
  
  if (!is.null(u_init)) {
    u <- u_init
  } else {
    u <- D %*% y
  }
  L <- ncol(B1)
  
  iter <- 0
  old_u <- u
  last_iter <- FALSE
  while (!last_iter) {
    iter <- iter + 1
    
    # browser()
    
    Lam <- Matrix::Diagonal(x = as.matrix(abs(u)) / lam)
    Q <- AAT + B1 %*% Lam %*% t(B1)
    
    # Solve Q system efficiently
    Qinv_B1_Lam_b <- solve(Q, B1 %*% (Lam %*% b))
    u <- Lam %*% (b - t(B1) %*% Qinv_B1_Lam_b)
    
    cost[iter] <- 0.5 * sum(abs(Hy - solve(A, B1 %*% u))^2) + lam * sum(abs(u))
    
    delta_u <- max(abs(u - old_u)) / (max(abs(old_u)) + .Machine$double.eps)
    old_u <- u
    
    if ((delta_u <= TOL_STOP) || (iter >= MAX_ITER)) {
      last_iter <- TRUE
    }
  }
  
  v <- solve(A, B1 %*% u)
  x <- y - Hy + v
  
  r <- t(B1) %*% solve(AAT, B %*% y - B1 %*% u) / lam
  k <- abs(r) > 1.02
  NZL <- sum(k)
  cat(NZL, "incorrectly locked zeros detected.\n")
  
  cost <- cost[1:iter]
  
  list(x = x, u = u, cost = cost, v = v)
}
