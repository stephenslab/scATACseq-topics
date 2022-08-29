# Run ashr on these results.
pdat <- transform(pdat,se = postmean/z)
fit <- ash(pdat$postmean,pdat$se,mixcompdist = "normal",method = "shrink")
pdat$ash <- fit$result$PosteriorMean
ggplot(pdat,aes(x = postmean,y = ash,color = se)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dashed") +
  theme_cowplot()
