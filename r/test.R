library(foreach)
library(doParallel)
library(doRNG)
source("toeplitz_mssa.R")
source("mcmssa_utils.R")

N <- 100
varphi <- 0.7
delta <- 1
omega <- 0.075
N <- 100
Ls <- c(10, 20, 50, 80, 90)
L_idx <- 1:length(Ls)
D <- 2 
G <- 1000
M <- 1000
model <- list(list(varphi = varphi,
                   delta = delta,
                   N = N),
              list(varphi = varphi,
                   delta = delta,
                   N = N))
signal <- sin(2*pi*(1:N) * omega)
plot(signal, type = "l")



signal <- replicate(D, 0, simplify=F)

set.seed(5)
p.values_noise.me1block.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", toeplitz.kind = "block", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1block.ev[[idx]] <- pvals
}

p.values_noise.me1block.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", toeplitz.kind = "block", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1block.fa[[idx]] <- pvals
}

p.values_noise.me1sum.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", toeplitz.kind = "sum", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1sum.ev[[idx]] <- pvals
}

p.values_noise.me1sum.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", toeplitz.kind = "sum", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1sum.fa[[idx]] <- pvals
}

signal <- list(signal.one.channel(model[[1]]$N, omega), 
               signal.one.channel(model[[1]]$N, omega))
p.values_signal.me1block.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", toeplitz.kind = "block", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1block.ev[[idx]] <- pvals
}

p.values_signal.me1block.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", toeplitz.kind = "block", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1block.fa[[idx]] <- pvals
}

p.values_signal.me1sum.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", toeplitz.kind = "sum", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1sum.ev[[idx]] <- pvals
}

p.values_signal.me1sum.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", toeplitz.kind = "sum", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1sum.fa[[idx]] <- pvals
}

alphas <- 0:1000/1000
alphas_idx <- 1:length(alphas)
clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

roc.me1sum.fa <- list()
alpha_1.me1sum.fa <- list()
beta.me1sum.fa <- list()
for (l in L_idx)
{
  roc.me1sum.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa[[l]] < alpha)/M
    beta.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

plot(roc.me1sum.fa[[4]]$fpr, roc.me1sum.fa[[4]]$tpr, type = 'l', lwd = lwds[4], xlim=c(0, 1), col = clrs[4], ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
lines(roc.me1sum.fa[[3]]$fpr[-1], roc.me1sum.fa[[3]]$tpr[-1], col = clrs[3], type = 'l', lwd = lwds[3])
for (l in L_idx[c(-1,-2,-3)])
  lines(roc.me1sum.fa[[l]]$fpr, roc.me1sum.fa[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)

plot(alpha_1.me1sum.fa[[1]]$alpha, alpha_1.me1sum.fa[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.fa[[l]]$alpha, alpha_1.me1sum.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)


roc.me1sum.ev <- list()
alpha_1.me1sum.ev <- list()
beta.me1sum.ev <- list()
for (l in L_idx)
{
  roc.me1sum.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev[[l]] < alpha)/M
    beta.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

plot(roc.me1sum.ev[[1]]$fpr, roc.me1sum.ev[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[-1])
  lines(roc.me1sum.ev[[l]]$fpr, roc.me1sum.ev[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)

plot(alpha_1.me1sum.ev[[1]]$alpha, alpha_1.me1sum.ev[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.ev[[l]]$alpha, alpha_1.me1sum.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)


roc.me1block.ev <- list()
alpha_1.me1block.ev <- list()
beta.me1block.ev <- list()
for (l in L_idx)
{
  roc.me1block.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    alpha_1.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev[[l]] < alpha)/M
    beta.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

plot(roc.me1block.ev[[1]]$fpr, roc.me1block.ev[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[c(-1,-4,-5)])
  lines(roc.me1block.ev[[l]]$fpr, roc.me1block.ev[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(roc.me1block.ev[[4]]$fpr[-1], roc.me1block.ev[[4]]$tpr[-1], col = clrs[4], type = 'l', lwd = lwds[4])
lines(roc.me1block.ev[[5]]$fpr[-1], roc.me1block.ev[[5]]$tpr[-1], col = clrs[5], type = 'l', lwd = lwds[5])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)

plot(alpha_1.me1block.ev[[1]]$alpha, alpha_1.me1block.ev[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.ev[[l]]$alpha, alpha_1.me1block.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)


roc.me1block.fa <- list()
alpha_1.me1block.fa <- list()
beta.me1block.fa <- list()
for (l in L_idx)
{
  roc.me1block.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    alpha_1.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa[[l]] < alpha)/M
    beta.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

plot(roc.me1block.fa[[1]]$fpr[-1], roc.me1block.fa[[1]]$tpr[-1], type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
lines(roc.me1block.fa[[2]]$fpr[-1], roc.me1block.fa[[2]]$tpr[-1], col = clrs[2], type = 'l', lwd = lwds[2])
lines(roc.me1block.fa[[3]]$fpr[-1], roc.me1block.fa[[3]]$tpr[-1], col = clrs[3], type = 'l', lwd = lwds[3])
lines(roc.me1block.fa[[4]]$fpr, roc.me1block.fa[[4]]$tpr, col = clrs[4], type = 'l', lwd = lwds[4])
lines(roc.me1block.fa[[5]]$fpr, roc.me1block.fa[[5]]$tpr, col = clrs[5], type = 'l', lwd = lwds[5])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)

plot(alpha_1.me1block.fa[[1]]$alpha, alpha_1.me1block.fa[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.fa[[l]]$alpha, alpha_1.me1block.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
