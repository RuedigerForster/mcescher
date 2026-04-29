# beads.R
# Baseline Estimation and Denoising with Sparsity (BEADS)
#
# Pure-R port of the MATLAB implementation by Xiaoran Ning (student of I. Selesnick).
# Original reference:
#   Ning, Selesnick, Duval (2014). Chromatogram baseline estimation and denoising
#   using sparsity (BEADS). Chemometrics and Intelligent Laboratory Systems.
#   doi: 10.1016/j.chemolab.2014.09.014
#
# Dependencies: Matrix

.phi_v1  <- function(x, EPS1) sqrt(x^2 + EPS1)
.wfun_v1 <- function(x, EPS1) 1 / sqrt(x^2 + EPS1)
.phi_v2  <- function(x, EPS1) abs(x) - EPS1 * log(abs(x) + EPS1)
.wfun_v2 <- function(x, EPS1) 1 / (abs(x) + EPS1)

.theta <- function(x, EPS0, r) {
  pos  <- x >  EPS0
  neg  <- x < -EPS0
  mid  <- !pos & !neg
  sum(x[pos]) -
    r * sum(x[neg]) +
    sum((1 + r) / (4 * EPS0) * x[mid]^2 + (1 - r) / 2 * x[mid] + EPS0 * (1 + r) / 4)
}

.H <- function(x, A, B) B %*% solve(A, x)

.BAfilt <- function(d, fc, N) {
  b1 <- c(1, -1)
  if (d > 1)
    for (i in seq_len(d - 1))
      b1 <- convolve(b1, rev(c(-1, 2, -1)), type = "open")

  b   <- convolve(b1, rev(c(-1, 1)), type = "open")
  omc <- 2 * pi * fc
  t   <- ((1 - cos(omc)) / (1 + cos(omc)))^d

  a <- 1
  for (i in seq_len(d))
    a <- convolve(a, rev(c(1, 2, 1)), type = "open")
  a <- b + t * a

  NbDiag <- d + 1L
  DiagA  <- vector("list", NbDiag)
  DiagB  <- vector("list", NbDiag)
  for (i in 0:d) {
    DiagA[[i + 1]] <- rep(a[d - i + 1], N - i)
    DiagB[[i + 1]] <- rep(b[d - i + 1], N - i)
  }

  A <- Matrix::bandSparse(N, k = 0:d, diagonals = DiagA, symm = TRUE)
  B <- Matrix::bandSparse(N, k = 0:d, diagonals = DiagB, symm = TRUE)
  list(A = A, B = B)
}

#' BEADS baseline estimation.
#'
#' @param y    Numeric vector  raw chromatogram signal.
#' @param d    Filter order (1 or 2).
#' @param fc   Filter cut-off frequency (cycles/sample, 0 < fc < 0.5).
#' @param r    Asymmetry ratio.
#' @param lam0 Regularisation parameter (sparsity of x).
#' @param lam1 Regularisation parameter (sparsity of first derivative).
#' @param lam2 Regularisation parameter (sparsity of second derivative).
#' @param pen  Penalty: "L1_v1" (smooth L1) or "L1_v2" (log-barrier L1).
#' @param Nit  Number of iterations (default 30).
#'
#' @return List: $x (denoised signal), $f (baseline), $cost (cost history).
#' @export
beads <- function(y, d = 1, fc = 0.01, r = 6,
                  lam0 = 0.5, lam1 = 5, lam2 = 4,
                  pen = "L1_v2", Nit = 30L) {

  EPS0 <- 1e-6
  EPS1 <- 1e-6

  if (pen == "L1_v1") {
    phi  <- function(x) .phi_v1(x,  EPS1)
    wfun <- function(x) .wfun_v1(x, EPS1)
  } else if (pen == "L1_v2") {
    phi  <- function(x) .phi_v2(x,  EPS1)
    wfun <- function(x) .wfun_v2(x, EPS1)
  } else {
    stop("pen must be 'L1_v1' or 'L1_v2'")
  }

  y <- as.numeric(y)
  x <- y
  N <- length(y)

  filt  <- .BAfilt(d, fc, N)
  A     <- filt$A
  B     <- filt$B
  BTB   <- crossprod(B)

  e   <- rep(1, N - 1)
  D1  <- Matrix::bandSparse(N - 1, N, k = c(0, 1), diagonals = list(-e,  e))
  D2  <- Matrix::bandSparse(N - 2, N, k = c(0, 1, 2),
                    diagonals = list(e, -2 * e, e))
  D   <- rbind(D1, D2)

  w <- c(rep(lam1, N - 1), rep(lam2, N - 2))
  b <- rep((1 - r) / 2, N)
  d_vec <- as.numeric(BTB %*% solve(A, y)) - lam0 * as.numeric(t(A) %*% b)

  gamma <- rep(1, N)
  cost  <- numeric(Nit)

  for (i in seq_len(Nit)) {
    Dxw    <- as.numeric(D %*% x)
    lam_d  <- w * wfun(Dxw)
    Lambda <- Matrix::Diagonal(x = lam_d)

    k <- abs(x) > EPS0
    gamma[!k] <- (1 + r) / (4 * EPS0)
    gamma[ k] <- (1 + r) / (4 * abs(x[k]))
    Gamma <- Matrix::Diagonal(x = gamma)

    M <- 2 * lam0 * Gamma + t(D) %*% Lambda %*% D
    x <- as.numeric(A %*% solve(BTB + t(A) %*% M %*% A, d_vec))

    resid <- y - x
    cost[i] <- 0.5 * sum(.H(resid, A, B)^2) +
               lam0  * .theta(x, EPS0, r) +
               lam1  * sum(phi(diff(x))) +
               lam2  * sum(phi(diff(x, differences = 2)))
  }

  f <- y - x - as.numeric(.H(y - x, A, B))
  list(x = x, f = f, cost = cost)
}
