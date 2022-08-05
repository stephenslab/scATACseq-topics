library(Matrix)
library(fastTopics)
load("data/Cusanovich_2018/processed_data/Cusanovich_2018_Kidney.RData")
fit <- readRDS(file.path("output/Cusanovich_2018/tissues",
                         "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit
