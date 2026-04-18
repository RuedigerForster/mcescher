MAD <- function(y, k = 1.4826) {
  # median absolute deviation
  d <- sapply(seq_len(ncol(y)), function(i) c(0, diff(y[, i])))
  d_med <- sapply(seq_len(ncol(y)), function(i) median(d[, i]))
  return(k * d_med * abs(d - d_med))
}


swarbled <- function(df, divs = 5, thresh = 3.29, blind = 200, exclude = NULL) {
  # SWitching ARtifacts BLock Enable / Disable
  # 1st derivatives of baselines
  # X <- scale(abs(rowMeans(sapply(seq_len(ncol(df)), function(x) c(0, diff(df[, x]))))))
  X <- scale(sapply(seq_len(ncol(df)), function(x) c(0, diff(df[, x]))))
  # select baselines which correlate best with 1st baseline
  idx <- sapply(seq_len(ncol(X)), function(x) cor(X[, x], X[, 1], method = "pearson", use = "complete.obs")) > 0.5
  MAD_X <- rowSums(MAD(X[, idx]))

  block <- vector(mode = "logical", length = nrow(df))
  block[MAD_X >= (quantile(MAD_X)[4] * thresh)] <- TRUE
  # if (!is.null(exclude) && (exclude[2] > exclude[1])) {
  #   block[exclude[1]:exclude[2]] <- FALSE
  # }

  runs <- c(0, rle(diff(block)))
  falling <- which(runs$values == -1)
  rising <- which(runs$values == 1)
  
  if (!identical(length(falling), length(rising))) browser()
  
  # eliminate orphans
  for (i in 1:(length(falling)-1)) {
    lo <- (sum(runs$lengths[1:falling[i]]) ) : (sum(runs$lengths[1:rising[i+1]]))
    if (length(lo) < blind) {
      block[lo] <- TRUE
    }
  }
  runs <- c(0, rle(diff(block)))
  falling <- which(runs$values == -1)
  rising <- which(runs$values == 1)
  good <- list()
  # good segments start with falling and end with rising
  for (i in 1:(length(rising)-1)) {
    begin <- sum(runs$lengths[1:falling[i]])
    end = sum(runs$lengths[1:rising[i+1]])
    # print(paste0(i, ": ", end - begin))
    if ((end - begin) > blind) {
      good[[i]] <- c(begin = begin, end = end)
    }
  }
  # last segment
  good[[i+1]] <- c(begin = sum(runs$lengths[1:falling[i+1]]), end = nrow(df))
  # first segment
  end <- sum(runs$lengths[1:rising[1]])
  if (end > blind) {
    good <- c(list(c(begin = 1, end = end)), good)
    block[1:end] <- TRUE
  }
  # reconfigure block
  block <- rep(TRUE, nrow(df))
  for (i in seq_along(good)) {
    block[good[[i]]['begin']:good[[i]]['end']] <- FALSE
  }
  
  return(list(block = block, segments = good))
}


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
  return(
    ResizeMat(df, ndim = c(nrow(df)/compression_factor[1], ncol(df)/compression_factor[2]))
  )
}

compress_vec <- function(vec, compression_factor) {
  ndim <- trunc(length(vec)/compression_factor)
  seq(from = vec[1], to = vec[length(vec)], length.out = ndim)
}
  

