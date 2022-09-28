# The idea with this function ash_test_enrich is to use ashr in a way
# so that it is more suitable as a test for enrichment. The key
# outputs are the likelihood ratio (logLR) which gives the strength of
# evidence for enrichment, and the mean posterior coefficient (coef),
# which summarizes whether the effects are consistent in direction (if
# they aren't, this will be close to zero).
ash_test_enrich <- function (b, se, g,
                             prior = rep(1,length(g$pi)),
                             lfsr.threshold = 0.05) {
  
  # Set up the z-scores output.
  n <- length(b)
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
  lfsr <- rep(0.5,n)
  lfsr[i] <- out$result$lfsr

  # Get the mean posterior enrichment coefficient among the
  # co-ordinates satisfying the lfsr threshold.
  i <- which(lfsr < lfsr.threshold)
  if (length(i) == 0)
    coef <- 0
  else
    coef <- mean(b[i])
  
  # Output the revised estimates (b), the standard errors (se), the
  # z-scores (z), the local false sign rates (lfsr), the likelihood
  # ratio (logLR), the mean posterior coefficient (coef), and the
  # updated prior weights (pi).
  return(list(b = b,se = se,z = z,lfsr = lfsr,logLR = out$logLR,
              coef = coef,pi = out$fitted_g$pi))
}
