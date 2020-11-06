
#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#

# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P.
# Microsimulation modeling for health decision sciences using R: A tutorial.
# Med Decis Making. 2018;38(3):400–22. https://www.ncbi.nlm.nih.gov/pubmed/29587047

################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS.
# All rights reserved in Canada, the United States and worldwide.
# Copyright, trademarks, trade names and any and all associated intellectual property
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating
# institutions and may not be used, reproduced, modified, distributed or adapted
# in any way without appropriate citation.

################################################################################
# Developed by Petros Pechlivanoglou

samplev <- function(m.Probs) {
  # Arguments
  # m.Probs: matrix with probabilities (n.i * n.s)
  # Return
  # ran: n.i x m matrix filled with sampled health state(s) per individual
  d <- dim(m.Probs) # dimensions of the matrix filled with the multinomical probabilities for the health states
  n <- d[1] # first dimension - number of rows (number of individuals to sample for)
  k <- d[2] # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]] # extract the names of the health states considered for sampling
  if (!length(lev)) # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k # create a sequence from 1:k (number of health states considered)
  # create a matrix
  ran <- rep(lev[1], n) # create the matrix ran, filled with the first health state of the levels
  U <- t(m.Probs) # transposed m.Probs matrix n.i x n.s --> n.s x n.i
  for(i in 2:k) { # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual
  }
  if (any((U[k, ] - 1) > 1e-05)) # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  un <- rep(runif(n), rep(k, n)) # sample from a uniform distribution of length n*k
  ran <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  ran # return the new health state per individual n.i x m
} # close the function
