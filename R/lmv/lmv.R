library(MASS)
library(nlme)
library(stats)
library(segmented)


# alocate outliers in the chromatogram based on SNR
# calculate SNR:
# SNR[i] = max(abs(x[i] − median(x[w])) / sigma), abs(x[i] − x[i−1]) / sigma
# outlier = SNR[i] > 2.5
#
is_outlier <- function(arr, len, crit = 2.5) {
  if (length(arr) != len) stop(sprintf("is_outlier: length mismatch (%d vs %d)", length(arr), len))
  sigma <- stats::mad(arr, na.rm = TRUE)
  med <- median(arr)
  SNR1 <- (abs(arr - med) / sigma)[-1]
  SNR2 <- diff(arr) / sigma
  idx <- pmax(SNR1, SNR2) > crit
  idx <- c(FALSE, idx)
  return(idx)
}

# allocate outliers in the 1st dericative of the chromatogram based on MAD
is_outlier_diff <- function(arr, crit = 2.5) {
  arr <- diff(arr)
  sigma <- stats::mad(arr, na.rm = TRUE)
  med <- median(arr)
  SNR1 <- (abs(arr - med) / sigma)[-1]
  SNR2 <- diff(arr) / sigma
  idx <- pmax(SNR1, SNR2) > crit
  idx <- c(FALSE, FALSE, idx)
  return(idx)
}

# interpolate outliers by segmented fitting
interpolate_outlier <- function(ind, arr) {
  m <- min(arr)
  arr <- arr - m
  st <- 1:length(arr)
  fit.glm <- glm(arr ~ st, family=gaussian)
  fit.seg<-segmented(fit.glm)
  arr[ind] <- fit.seg$fitted.values[ind]
  arr <- arr + m
  return(arr)
}

# outlier is replaced by the median of the window
substitute_outlier <- function(ind, arr, len) {
  if((ind - len/2) > 0) {
    begin <- ind - len/2
    shift <- 0
  } else {
    begin <- 0
    shift <- ind - len/2
  }
  if((ind + len/2 + shift) < length(arr)) {
    end <- ind - shift + len/2
  } else {
    end <- length(arr)
    shift <- length(arr) - (ind + len/2 + shift)
    begin <- begin + shift
  }
  # if ((end - begin) != len) {
  #   print(paste(ind, "window size:", end - begin))
  # }
  ex <- median(arr[begin:end])
  ## R function for this is:
  ## replace(x, list, values)
  arr[ind] <- ex
  return(arr)
}


lmv_rsa <- function(signal, baseline = NULL, crit = 2.5, span = 50, max_iter = 20) {
  # stage 1: extract local minimum values and record their positions
  # establish the minimum vector, where:
  # 𝑥[𝑖−1] > 𝑥[𝑖], and
  # 𝑥[𝑖] < 𝑥[𝑖+1]
  #
  # browser()
  l <- length(signal)
  x <- signal
  index <- seq(1:length(x))
  falling <- c(TRUE, x[-length(x)] > x[-1])
  rising  <- c(x[-length(x)] < x[-1], TRUE)
  min_vec <- which(c(rising & falling))
  x <- signal[min_vec]
  if (length(min_vec) == 0) stop("Cannot detect local minima.")
  # stage 2a: Detect outliers points by using moving window
  # and replace the outlier points with median value in the window
  #
  o <- lapply(1:(length(min_vec)-span), function(i) seq(from = i, to = i + span - 1, by = 1)) 
  ot <- list()
  ###
  ## this snippet does not belong to the LMV algorithm
  library(runner)
  x <- runner::runner(
    x = signal,
    k = span,
    lag = span, 
    f = function(x) min(x)
  )[min_vec]
  idx <- which(is.na(x))
  x[idx] <- signal[idx]
  if (sum(is.na(x)) > 0) stop("lmv_rsa: NA values remain after runner imputation")
  if (length(x) != length(min_vec)) stop("lmv_rsa: length mismatch between x and min_vec")
  ###
  new_x1 <- x
  iter <- 0
  repeat {
    for (i in seq_along(o)) {
      ot[[i]] <- which(is_outlier(new_x1[o[[i]]], span))
      for (j in seq_along(ot[[i]])) {
        new_x1 <- substitute_outlier(o[[i]][ot[[i]][j]], new_x1, span)
      }
    }
    frobenius_norm <-
      sum((new_x1-x)^2, na.rm = TRUE) / sum(x^2, na.rm = TRUE)
    cat(".")
    if((frobenius_norm < 1e-4) | (iter > max_iter)) break
    x <- new_x1
  }
  cat("\n")
  
  # stage 2b: Detect outlier points in the ﬁrst-order derivative of x
  # and replace outlier points with new values by linear interpolation
  #
  o <- lapply((0:trunc((length(min_vec)-span) / span)), function(i) 
    seq(from = i * span + 1, to = (i+1) * span, by = 1))
  o[[length(o) + 1]] <- max(unlist(o)):length(min_vec)
  x <- signal[min_vec]
  new_x2 <- x
  pt <- list()
  iter <- 1
  repeat {
    for (i in seq_along(o)) {
      # print(paste0("[[", i, "]]"))
      pt[[i]] <- which(is_outlier_diff(x[o[[i]]]))
      new_x2[o[[i]]] <- interpolate_outlier(pt[[i]], x[o[[i]]]) 
    }
    frobenius_norm <-
      sum((new_x2-x)^2, na.rm = TRUE) / sum(x^2, na.rm = TRUE)
    cat(".")
    if((frobenius_norm < 1e-4) | (iter > max_iter)) break
    x <- new_x2
    iter <- iter + 1
  }
  cat("\n")
  
  # obtain a minimum vector from results of stage 2a and 2b
  #
  idx <- new_x1 < new_x2
  new_x2[idx] <- new_x1[idx]

  # fill in the gaps by linear interpolation
  # (x[i] - x[i-1]) / (index[i] - index[i-1]) * (j - p[i-1]) + x[i-1]
  #
  if (any(rle(min_vec)$lengths > 1)) stop() # this should never happen
  x <- approx(min_vec, new_x2, xout = seq_len(l))$y
  idx <- is.na(signal)
  x[idx] <- NA
  return(x)
}


# # input parameters
# crit <- 2.5
# span <- 50
# dat <- df2[, 35] - df3[, 35]
# l <- length(x)
# plot(dat, type = "l", ylim = c(-0.2, 0.2))
# x <- lmv_rsa(dat, crit = 2.5, span = 50)
# col <- "red"
# lines(x, col = col, lwd = 1)
# 
# plot(x - dat, type = "l", ylim = c(-0.5, 0.1))
# 
# 

