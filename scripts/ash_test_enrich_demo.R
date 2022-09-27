library(ashr)
source("../code/ash.R")
set.seed(1)
b <- c(rep(1,10),rep(0,1990))
se <- abs(rnorm(2000))
fit0 <- ash(b,se,mixcompdist = "normal",pointmass = FALSE,
            outputlevel = 1)
g0 <- fit0$fitted_g
logLR <- NULL
coefs <- NULL
for (i in seq(0,1990,10)) {
  j <- i + (1:10)
  out <- ash_test_enrich(b[j],se[j],g0,prior = 1.01 + g0$pi)
  logLR <- c(logLR,out$logLR)
  coefs <- c(coefs,out$coef)
}
layout(matrix(1:2,2,1))
plot(1:200,logLR,pch = 20)
plot(1:200,coefs,pch = 20)
