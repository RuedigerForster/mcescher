# parameter w: max shiftng width
w <- 10

x <- dat <- df2[, 10] - df3[, 15]
l <- length(x)
plot(x, type = "l", ylim = c(0, 2.2), xlim = c(2600, 2800))

dx <- c(0, diff(x))
dx[is.na(dx)] <- 0
index <- 1:l
# plot(index, dx, type = "l")

A <- N <- P <- rep(0, l)
pos <- dx > 0
neg <- dx < 0
P[pos] <- dx[pos]
N[neg] <- dx[neg]
 

lines(index, P, col = "red", type = "l", ylim = c(-0.2, 0.2), xlim = c(1000, 2000))
lines(index, -N, col = "darkgreen")

accu <- function(w) {
  A <- rep(0, l)
  for(i in index) {
    s <- ifelse((i-w) > 0, w, 0)
    A[i] <- sum(P[(i-s):i], na.rm = T) + sum(-N[i:(i+s)], na.rm = T)
  }
  lines(A)
  return(A)
}

# plot(1:trunc(l/8), rep(0, l/2), type = "l", ylim = c(-1, 5))
O <- numeric()
for (i in 1:30) {
  O[i] <- max(accu(i))
  print(paste0(i, ": ", O[i]))
  lines(O[i])
}

which(diff(O) <= 0)
 plot(O)

plot(accu(7) / 2, col = "red", type = "l", ylim = c(-2, 3))
lines(x)






plot(index[pos], dx[pos], col = "red", type = "l", ylim = c(-0.2, 0.2))
lines(index[neg], dx[neg], col = "darkgreen")

