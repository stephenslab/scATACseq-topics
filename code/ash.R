# Perform adaptive shrinkage on the matrix of effect estimates b and
# their standard errors, then output the revised effect estimates (b),
# standard errors (se), z-scores (z) and local false sign rates
# (lfsr). All effects i in which either b[i] or se[i] is missing (NA)
# are not revised.
#
# The specific ash options used here are intended to work well for the
# shrinkage of the de_analysis results for peaks nearby a gene.
#
shrink_estimates <- function (b, se, g, fixg = TRUE) {
  
  # Set up the z-scores output.
  z <- b
  z[is.na(b) | is.na(se)] <- as.numeric(NA)
  
  # Run adaptive shrinkage.
  i   <- which(!(is.na(b) | is.na(se)))
  out <- ash(b[i],se[i],mixcompdist = "normal",method = "shrink", 
             g = g,fixg = fixg)
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
  rownames(lfsr) <- rownames(b)
  colnames(lfsr) <- colnames(b)

  # Output the revised estimates (b), the standard errors (se), the
  # z-scores (z), the local false sign rates (lfsr) and the raw ash
  # output (ash).
  return(list(b = b,se = se,z = z,lfsr = lfsr,ash = out))
}
