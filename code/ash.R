# Perform adaptive shrinkage on the matrix of effect estimates b and
# their standard errors, then output the revised effect estimates (b),
# standard errors (se), z-scores (z) and local false sign rates
# (lfsr). All effects i in which either b[i] or se[i] is missing (NA)
# are not revised.
shrink_estimates <- function (b, se, mixsd, prior = rep(4,length(mixsd))) {
  
  # Set up the z-scores output.
  z <- b
  z[is.na(b) | is.na(se)] <- as.numeric(NA)
  
  # Run adaptive shrinkage.
  i   <- which(!(is.na(b) | is.na(se)))
  out <- ash(b[i],se[i],mixcompdist = "normal",mixsd = mixsd,
             prior = prior,...)
  b1  <- out$result$PosteriorMean
  se1 <- out$result$PosteriorSD
  
  # Extract the posterior estimates and their standard errors. 
  b[i]  <- b1
  se[i] <- se1
  z[i]  <- b[i]/se[i]
  z[i[b1 == 0]]  <- 0
  z[i[se1 == 0]] <- as.numeric(NA)

  # Extract the lfsr estimates.
  m       <- nrow(b)
  k       <- ncol(b)
  lfsr    <- matrix(as.numeric(NA),m,k)
  lfsr[i] <- out$result$lfsr

  # Output the revised estimates (b), the standard errors (se), the
  # z-scores (z), and the local false sign rates (lfsr).
  return(list(b = b,se = se,z = z,lfsr = lfsr))
}
