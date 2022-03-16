library(fastTopics)
fit <- readRDS("fit-Buenrostro2018-binarized-scd-ex-k=10.rds")$fit
fit <- poisson2multinom(fit)
heatmap(cor(fit$F))
k1 <- 5
k2 <- 8
i <- which(with(fit,F[,k1] > 1e-5 | F[,k2] > 1e-5))
plot(fit$F[i,k1],fit$F[i,k2],pch = 20)
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

D <- matrix(0,10,10)
for (i in 1:10)
  for (j in 1:10)
    D[i,j] <- quantile(abs(log2(fit$F[,i] + 1e-6) - log2(fit$F[,j] + 1e-6)),
                       0.99)
heatmap(cor(fit$F),symm = TRUE)
