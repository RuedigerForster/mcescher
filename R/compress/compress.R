library(fields)

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  return(newrange[1]+(x-xrange[1])*mfac)
}

ResizeMat <- function(mat, ndim=dim(mat)){
  if(!require(fields)) stop("`fields` required.")
  
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  
  # interpolation
  ans[ncord] <- interp.surface(obj, loc)
  
  ans
}


ResizeVec <- function(vec, ndim=length(vec)){
  seq(from = vec[1], to = vec[length(vec)], length.out = ndim)
}

compress_mat <- function(df, compression_factor) {
  if (ncol(df) == 1L) {
    xf <- compression_factor[1]
    n  <- nrow(df) %/% xf
    m  <- matrix(colMeans(matrix(as.numeric(df[1:(n * xf), 1]), nrow = xf)), ncol = 1L)
    return(m)
  }
  return(
    ResizeMat(df, ndim = c(nrow(df)/compression_factor[1], ncol(df)/compression_factor[2]))
  )
}

compress_vec <- function(vec, compression_factor) {
  ndim <- trunc(length(vec)/compression_factor)
  seq(from = vec[1], to = vec[length(vec)], length.out = ndim)
}

decompress_mat <- function(mat, n_full) {
  if (ncol(mat) == 1L) {
    xout <- seq(1, nrow(mat), length.out = n_full)
    return(matrix(approx(seq_len(nrow(mat)), mat[, 1], xout = xout)$y, ncol = 1L))
  }
  ResizeMat(mat, ndim = c(n_full, ncol(mat)))
}
