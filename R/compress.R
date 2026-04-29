.rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  return(newrange[1]+(x-xrange[1])*mfac)
}

.interp_bilinear <- function(obj, loc) {
  x <- obj$x; y <- obj$y; z <- obj$z
  px <- loc[, 1]; py <- loc[, 2]
  ix <- pmax(1L, pmin(findInterval(px, x, rightmost.closed = TRUE), length(x) - 1L))
  iy <- pmax(1L, pmin(findInterval(py, y, rightmost.closed = TRUE), length(y) - 1L))
  tx <- (px - x[ix]) / (x[ix + 1L] - x[ix])
  ty <- (py - y[iy]) / (y[iy + 1L] - y[iy])
  z[cbind(ix,      iy     )] * (1-tx) * (1-ty) +
  z[cbind(ix + 1L, iy     )] *    tx  * (1-ty) +
  z[cbind(ix,      iy + 1L)] * (1-tx) *    ty  +
  z[cbind(ix + 1L, iy + 1L)] *    tx  *    ty
}

.ResizeMat <- function(mat, ndim=dim(mat)){
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = .rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = .rescale(ncord[,2], c(1,odim[2]))
  ans[ncord] <- .interp_bilinear(obj, loc)
  ans
}


.ResizeVec <- function(vec, ndim=length(vec)){
  seq(from = vec[1], to = vec[length(vec)], length.out = ndim)
}

.compress_mat <- function(df, compression_factor) {
  xf <- compression_factor[1]
  n  <- nrow(df) %/% xf
  nc <- ncol(df)
  # Block-average each column: reshape each column into (xf  n) and take column means.
  # This correctly averages all samples within each block, preserving narrow peak areas.
  # The old multi-column path used ResizeMat (bilinear sampling), which aliases narrow
  # peaks  a peak narrower than xf samples would land between sample grid points and
  # effectively vanish.
  m <- vapply(seq_len(nc), function(j)
    colMeans(matrix(as.numeric(df[1:(n * xf), j]), nrow = xf)),
    numeric(n))
  matrix(m, nrow = n, ncol = nc)
}

.compress_vec <- function(vec, compression_factor) {
  ndim <- trunc(length(vec)/compression_factor)
  seq(from = vec[1], to = vec[length(vec)], length.out = ndim)
}

.decompress_mat <- function(mat, n_full) {
  xout <- seq(1, nrow(mat), length.out = n_full)
  x_in <- seq_len(nrow(mat))
  m <- vapply(seq_len(ncol(mat)), function(j)
    approx(x_in, mat[, j], xout = xout, rule = 2)$y,
    numeric(n_full))
  matrix(m, nrow = n_full, ncol = ncol(mat))
}
