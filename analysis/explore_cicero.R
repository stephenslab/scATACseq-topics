# Load the K = 10 topic model fit.
fit <- readRDS(
  file.path("../output/Cusanovich_2018/tissues",
            "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit
fit <- poisson2multinom(fit)

