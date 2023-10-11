library("Rssa")
library("dplyr")
library("pracma")

Lcov <- function(f1, f2, K) {
  f1 <- as.vector(f1)
  f2 <- as.vector(f2)
  N <- length(f1)
  c1 <- sapply(0:(K-1), function(i) sum(f1 * lag(f2, i), na.rm = TRUE) / (N - i))
  c2 <- sapply(0:(K-1), function(i) sum(lag(f1, i) * f2, na.rm = TRUE) / (N - i))
  Toeplitz(c1, c2)
}

toeplitz.mssa <- function(ts, L, D, method = c("sum", "block")) {
  #browser()
  N <- dim(ts)[1]
  #D <- dim(ts)[2]
  K <- N - L + 1
  this <- list("F" = ts, N = N, L = L, K = K, D = D)
  traj.mat <- new.hbhmat(ts, L = c(L, 1))
  if (method == "block") {
    C <- rbind()
    for (i in 1:D) {
      c <- cbind()
      for (j in 1:D) {
        c <- cbind(c, Lcov(ts[, i], ts[, j], K))
      }
      C <- rbind(C, c)
    }
    traj.mat <- t(traj.mat)
  } else if (method == "sum") {
    C <- matrix(0, nrow = L, ncol = L)
    for (i in 1:D) {
      C <- C + Lcov(ts[, i], ts[, i], L)
    }
  }
  else
    stop('method should be one of "block", "sum"')
  S <- eigen(C, symmetric = TRUE)
  U <- S$vectors
  Z <- crossprod(traj.mat, U)
  sigma <- apply(Z, 2, function(x) sqrt(sum(x^2)))
  V <- sweep(Z, 2, sigma, FUN = "/")
  o <- order(sigma[seq_len(min(50, L, K))], decreasing = TRUE)
  sigma <- sigma[o]
  U <- U[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]
  if (method == "block") {
    this$U <- V
    this$V <- U
  }
  else {
    this$U <- U
    this$V <- V
  }
  this$sigma <- sigma
  this
}

diag.avg <- function(x, group) {
  X <- matrix(0, nrow = x$L, ncol = x$D * x$K)
  for (i in group) {
    X <- X + x$sigma[i] * x$U[, i] %*% t(x$V[, i])
  }
  series <- matrix(0, nrow = x$N, ncol = x$D)
  for (i in 1:x$D) {
    M <- X[, ((i-1)*x$K + 1):(i*x$K)]
    if (x$K < x$L)
      M <- t(M)
    series[, i] <- sapply(0:(x$N-1), function(k) mean(M[row(M) + col(M) == k + 2]))
  }
  ts(series)
}

toeplitz.reconstruct <- function(x, groups) {
  residuals <- x$F
  out <- list()
  for (i in seq_along(groups)) {
    out[[names(groups)[i]]] <- diag.avg(x, groups[[i]])
    residuals <- residuals - out[[names(groups)[i]]]
  }
  out[['Residuals']] <- residuals
  out
}