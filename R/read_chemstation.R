library(chromConverter)

source("./compress/compress.R")
source("./swarbled/new_swarbled.R")
source("./sass/sass.R")
source("./lmv/lmv.R")


bas <- list()
sass_res <- list()
path <- "/home/nvidiauser/Developments/OGE_UDE/green stick/April/OGE_AB 2022-04-28 14-01-55/0428_0_A.D/FID3A.ch"
if (file.exists(path)) {
  bas[[1]] <- read_chemstation_ch(path)
}
path <- "/home/nvidiauser/Developments/OGE_UDE/green stick/April/OGE_AB 2022-04-28 14-01-55/0428_0_B.D/FID3A.ch"
if (file.exists(path)) {
  bas[[2]] <- read_chemstation_ch(path)
}
bas <- sapply(bas, cbind)
block <- swarbled(rowMeans(bas), block.dead.time = TRUE)

# Define low-pass filter parameters
d <- 2; fc <- 0.022
# SASS - run algorithm
K <- 1; lam <- 1.2
sass_res <- sass_L1(rowMeans(bas), d, fc, K, lam)


path <- "/home/nvidiauser/Developments/OGE_UDE/green stick/April/OGE_AB 2022-04-28 14-01-55/0428_8_B.D/FID3A.ch"
if (file.exists(path)) {
  data <- read_chemstation_ch(path)
}

data <- (data - sass_res$x@x)[!block]
# input parameters
crit <- 2.5
window_size <- 100
# plot(data, type = "l", ylim = c(-1e-5, 1e-4))
x <- lmv_rsa(data, span = window_size, crit = 2.5)
plot(bas[!block, 1] - sass_res$x@x[!block], type = "l", ylim = c(-2e-5, 2e-5))
col <- "red"
plot(x - bas[!block, 2] + sass_res$x@x[!block], type = "l", col = col)
abline(h = 0, lty = "dashed")
col <- "green"
lines(x, col = col, lwd = 1)

