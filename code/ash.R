# TOO DO: Explain here what this function is for, and how to use it.
#
ash_test_enrich <- function (b, se, g, prior = "uniform") {
  
  # Set up the z-scores output.
  z <- b
  z[is.na(b) | is.na(se)] <- as.numeric(NA)
  
  # Run adaptive shrinkage.
  i   <- which(!(is.na(b) | is.na(se)))
  out <- ash(b[i],se[i],mixcompdist = "normal",pointmass = FALSE,
             g = g,prior = prior,outputlevel = 2)
  b1  <- out$result$PosteriorMean
  se1 <- out$result$PosteriorSD
  
  # Extract the posterior estimates and their standard errors. 
  b[i]  <- b1
  se[i] <- se1
  z[i]  <- b[i]/se[i]
  z[i[b1 == 0]]  <- 0
  z[i[se1 == 0]] <- as.numeric(NA)

  # Extract the lfsr estimates.
  lfsr <- out$result$lfsr

  # Output the revised estimates (b), the standard errors (se), the
  # z-scores (z), the local false sign rates (lfsr), the likelihood
  # ratio (logLR), the weighted posterior coefficient (coef), and the
  # updated prior weights (pi).
  i <- which(lfsr < 0.05)
  if (length(i) == 0)
    coef <- 0
  else
    coef <- mean(b[i])
  return(list(b = b,se = se,z = z,lfsr = lfsr,logLR = out$logLR,
              coef = coef,pi = out$fitted_g$pi))
}
