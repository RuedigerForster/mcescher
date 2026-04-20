dpa <- function(x, max.shift = 30) {
  
  accu <- function(w) {
    # browser()
    A <- rep(0, l)
    belly <- function(x, sigma, mu) {1 / sigma / sqrt(2 * 3.14159269) * exp(-(x - mu)^2  / 2 / sigma^2)}
    weights <- belly(1:(w+1), w/8, w/2)
    
    for(i in seq(w, l - w, 1)) {
      # s <- ifelse((i-w) > 0, w, 0)
      # gaussian weighted sum
      A[i] <- sum(P[(i-w):i] * weights) + sum(-N[i:(i+w)] * rev(weights))
      
      # A[i] <- (sum(P[(i-s):i], na.rm = TRUE) / SP + sum(-N[i:(i+s)], na.rm = TRUE) / SN) * (SP + SN)
    }
    return(A)
  }
  
  l <- length(x)
  index <- 1:l
  N <- P <- rep(0, l)
  
  dx <- c(0, diff(x))
  dx[is.na(dx)] <- 0
  pos <- dx > 0
  neg <- dx < 0
  P[pos] <- dx[pos]
  N[neg] <- dx[neg]
  SP <- max(cumsum(P))
  SN <- max(cumsum(-N))
  
  # print(paste(SP, SN, max(cumsum(P)), max(cumsum(-N))))
  
  O <- numeric()
  for (i in seq(max.shift)) {
    O[i] <- max(accu(i)) # accu(i)[29067] # 
  }
  # plot(O)
  first_max <- which(diff(log(O)) <= 0)[1] + 1
  print(first_max)
  # f <- sum(P, na.rm = T) / sum(-N, na.rm = T)
  return(accu(first_max) / 2)
}


dat <- df[, 1] # - df3[, 15]
plot(dat, type = "l", ylim = c(-10, 250), xlim = c(26000, 29500))
lines(dpa(dat, 30), col = "red")
plot(dpa(dat, 30), col = "red", type = "l", ylim = c(-0.01, 0.2))
# sigma <- 4
# mu <- 8
# X <- 1:16 
# f <- function(x) {1 / sigma / sqrt(2 * 3.14159269) * exp(-(x - mu)^2  / 2 / sigma^2)}
# plot(f(X))
