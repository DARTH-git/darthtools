#' Fit multiple survival models survival data
#'
#' \code{fit.fun} fits multiple survival models to survival data using flexsurv.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param data dataframe containing the time and status variables.
#' @param extrapolate extrapolate beyond model time horizon.
#' Default = FALSE.
#' @param times time horizon the extrapolation is done over.
#' @return
#' a list containing all survival model fits and a plot of superimposed survival curves.
#' @export
fit.fun <- function(time, status, data = data, extrapolate = FALSE, times) {
  require(survHE)
  # Extract the right data columns
  data$time   <- data[,   time]
  data$status <- data[, status]

  # Specify time horizon based on whether to extrapolate
  if (extrapolate == TRUE)  {
    plot.times <- times
  } else if  (extrapolate == FALSE) {
    plot.times <- seq(min(times),max(data$time), by = diff(times)[1])
  }

  # Fit parametric survival models
  # Define the vector of models to be used
  mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz")
  # Run the models using MLE via flexsurv
  fit.survHE <- fit.models(formula = Surv(time, status) ~ 1, data = data, distr = mods)

  # Run spline model via flexsurvspline
  fit.survspline <- flexsurvspline(formula = Surv(time, status) ~ 1, data = data, k = 2)

  # Extrapolate all models beyond the KM curve and plot
  # print(plot(fit.survHE, add.km = T, t = plot.times))
  KM.fit <- survfit(Surv(time, status) ~ 1, data = data) # fit Kaplan-Meier curve
  plot(KM.fit, ylab = "Survival", xlab = "Time", ylim = c(0,1), xlim = c(0, max(plot.times)), conf.int = F, mark.time = T)
  lines(fit.survHE$models$Exponential,     t = plot.times, col = 2, ci = F)
  lines(fit.survHE$models$`Weibull (AFT)`, t = plot.times, col = 3, ci = F)
  lines(fit.survHE$models$Gamma,           t = plot.times, col = 4, ci = F)
  lines(fit.survHE$models$`log-Normal`,    t = plot.times, col = 5, ci = F)
  lines(fit.survHE$models$`log-Logistic`,  t = plot.times, col = 6, ci = F)
  lines(fit.survHE$models$Gompertz,        t = plot.times, col = 7, ci = F)
  lines(fit.survspline,                    t = plot.times, col = 8, ci = F)
  # add a legend
  legend("topright", cex = 0.7, c("Kaplan-Meier", names(fit.survHE$models), "Spline"), col = 1:(length(mods)+2), lty = rep(1, (length(mods)+2)), bty="n")

  # Compare AIC values
  AIC <- c(fit.survHE$model.fitting$aic,  # parametric models
           AIC(fit.survspline))           # spline model
  AIC <- round(AIC,3)

  # Compare BIC values
  BIC <- c(fit.survHE$model.fitting$bic,  # parametric models
           BIC(fit.survspline))           # spline model
  BIC <- round(BIC,3)

  names(AIC) <- names(BIC) <- c(names(fit.survHE$models),
                                "Spline")

  # Store and return results
  res <- list(fit.survHE = fit.survHE,
              models     = append(fit.survHE$models, c("Spline" = list(fit.survspline))),
              AIC        = AIC,
              BIC        = BIC)
  return(res)
}

#' Fit multiple mixture cure models survival data
#'
#' \code{fit.fun.cure} fits multiple mixure cure models to survival data using flexsurv.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param data dataframe containing the time and status variables.
#' @param extrapolate extrapolate beyond model time horizon.
#' @param times time horizon the extrapolation is done over.
#' @return
#' a list containing all survival model fits and a plot of superimposed survival curves.
#' @export
fit.fun.cure <- function(time, status, data = data, extrapolate = FALSE, times) {
  require(survHE)
  # Extract the right data columns
  data$time   <- data[,   time]
  data$status <- data[, status]

  # Specify time horizon based on whether to extrapolate
  if (extrapolate == TRUE)  {
    plot.times <- times
  } else if  (extrapolate == FALSE) {
    plot.times <- seq(min(times),max(data$time), by = diff(times)[1])
  }

  # Fit parametric survival models
  # Define the vector of models to be used
  mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz")

  # Run mixture cure models via flexsurvcure
  fit.survcure <- lapply(mods, function(x) flexsurvcure(formula = Surv(time, status) ~ 1, data = data, dist = x, mixture = T))
  names(fit.survcure) <- paste(c("Exponential", "Weibull (AFT)", "Gamma", "log-Normal", "log-Logistic", "Gompertz"), "Cure", sep = " ")

  # Extrapolate all models beyond the KM curve and plot - cure models and spline models
  KM.fit <- survfit(Surv(time, status) ~ 1, data = data) # fit Kaplan-Meier curve
  plot(KM.fit, ylab = "Survival", xlab = "Time", ylim = c(0,1), xlim = c(0, max(plot.times)), conf.int = F, mark.time = T)
  lines(fit.survcure$`Exponential Cure`,   t = plot.times, col = 2, ci = F)
  lines(fit.survcure$`Weibull (AFT) Cure`, t = plot.times, col = 3, ci = F)
  lines(fit.survcure$`Gamma Cure`,         t = plot.times, col = 4, ci = F)
  lines(fit.survcure$`log-Normal Cure`,    t = plot.times, col = 5, ci = F)
  lines(fit.survcure$`log-Logistic Cure`,  t = plot.times, col = 6, ci = F)
  lines(fit.survcure$`Gompertz Cure`,      t = plot.times, col = 7, ci = F)
  # add a legend
  legend("topright", cex = 0.7, c("Kaplan-Meier", names(fit.survHE$models), "Spline"), col = 1:(length(mods)+1), lty = rep(1, (length(mods)+1)), bty="n")

  # Compare AIC values
  AIC <- unlist(lapply(fit.survcure, function(x) AIC(x)))  # cure models
  AIC <- round(AIC,3)

  # Compare BIC values
  BIC <- unlist(lapply(fit.survcure, function(x) BIC(x)))  # cure models
  BIC <- round(BIC,3)

  names(AIC) <- names(BIC) <- names(fit.survcure)

  # Store and return results
  res <- list(fit.survcure = fit.survcure,
              AIC          = AIC,
              BIC          = BIC)
  return(res)
}

#' Fit partitioned survival model to survival data
#'
#' \code{partsurv} fits partitioned survival model to survival data.
#'
#' @param pfs_survHE survHE obj fitting PFS.
#' @param os_survHE survHE obj fitting OS.
#' @param choose_PFS preferred PFS distribution.
#' @param choose_OS preferred OS distribution.
#' @param time numeric vector of time to estimate probabilities.
#' @param PA run probabilistic analysis.
#' Default = FALSE.
#' @param n_sim number of PA simulations.
#' Default = 100.
#' @param seed seed for random number generation.
#' Default = 421.
#' @return
#' a list containing Markov trace, expected survival, survival probabilities, transition probabilities.
#' @export
partsurv <- function(pfs_survHE, os_survHE, choose_PFS, choose_OS, time = times, PA = FALSE, n_sim = 100, seed = 421){
  set.seed(seed)
  deter <- ifelse(PA == 1, 0, 1) # determine if analysis is deterministic or probabilistic
  chosen_models <- paste0("PFS: ", choose_PFS, ", ", "OS: ", choose_OS) # chosen model names

  # fit parametric survival models for PFS and OS
  mod.pfs <- names(pfs_survHE$fit.survHE$models)
  mod.os  <- names(os_survHE$fit.survHE$models)

  # Calculate survival probabilities
  if (deter == 0) { # probabilistic
    fit_PFS <- make.surv(pfs_survHE$fit.survHE,
                         mod = which(mod.pfs == choose_PFS),
                         nsim = n_sim,
                         t = times)
    fit_OS  <- make.surv(os_survHE$fit.survHE,
                         mod = which(mod.os == choose_OS),
                         nsim = n_sim,
                         t = times)
    pfs.surv <- surv_prob(fit_PFS, PA = TRUE)
    os.surv  <- surv_prob( fit_OS, PA = TRUE)
  } else { # deterministic
    pfs.surv <- surv_prob(pfs_survHE$models[[which(mod.pfs == choose_PFS)]])
    os.surv  <- surv_prob( os_survHE$models[[which(mod.os  == choose_OS)]])
  }

  # Calculate area under survival curve (expected survival)
  if (deter == 0) { # probabilistic
    pfs.expected.surv <- mean(apply(pfs.surv, 2, function(x){expected_surv(time = times, surv = x)}))
    os.expected.surv  <- mean(apply(os.surv,  2, function(x){expected_surv(time = times, surv = x)}))
  } else { # deterministic
    pfs.expected.surv <- expected_surv(time = times, surv = pfs.surv)
    os.expected.surv  <- expected_surv(time = times, surv = os.surv)
  }

  # Calculate transition probabilities
  if (deter == 0) { # probabilistic
    pfs.trans.prob <- rowMeans(apply(pfs.surv, 2, trans_prob))
    os.trans.prob  <- rowMeans(apply( os.surv, 2, trans_prob))
  } else { # deterministic
    pfs.trans.prob <- trans_prob(pfs.surv)
    os.trans.prob  <- trans_prob(os.surv)
  }

  # Calculate state occupation proportions
  sick                 <- os.surv - pfs.surv    # estimate the probability of remaining in the progressed state
  sick[sick < 0]       <- 0                     # in cases where the probability is negative replace with zero
  healthy              <- pfs.surv              # probability of remaining stable
  dead                 <- 1 - os.surv           # probability of being dead
  trace <- abind(healthy,
                 sick,
                 dead, rev.along = 0)

  # Calculate Markov trace
  if (deter == 0){ # probabilistic
    trace      <- aperm(trace, perm = c(1,3,2)) # Markov trace of all simulations
    mean.trace <- apply(trace, 1:2, mean)       # average Markov trace across simulations
    CI         <- apply(trace, 1:2, quantile, probs = c(0.025, 0.975)) # trace confidence intervals
    CI         <- aperm(CI, perm = c(2,3,1))
    dimnames(mean.trace)[[2]] <- c("healthy", "sick", "dead")
    dimnames(CI)[[3]] <- c("low", "high")
    dimnames(CI)[[2]] <- c("healthy", "sick", "dead")
    dimnames(trace)[[2]] <- c("healthy", "sick", "dead")

    # List of items to return
    res <- list(trace = trace, CI = CI, Mean = mean.trace, # Markov trace
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  } else { # deterministic
    dimnames(trace)[[2]] <- c("healthy", "sick",  "dead")

    # List of items to return
    res <- list(trace = trace, # Markov trace
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  }

  return(res)
}

#' Calculate survival probabilities.
#'
#' \code{surv_prob} calculates survival probabilities from survival models.
#'
#' @param model survival model.
#' @param PA run probabilistic analysis.
#' Default = FALSE.
#' @return
#' vector of survival probabilities.
#' @export
surv_prob <- function(model, PA = FALSE) {
  if (PA) {
    surv <- model$mat[[1]][,-1]
  } else {
    surv <- summary(model, t = times, ci = F)[[1]]$est
  }
  return(surv)
}

#' Calculate transition probabilities.
#'
#' \code{trans_prob} calculates transition probabilities using survival probabilities.
#'
#' @param surv vector of survival probabilities.
#' @return
#' matrix of transition probabilities
#' @export
trans_prob <- function(surv){
  t.p <- 1- surv[-1]/(surv[-length(surv)])
  return(t.p = t.p)
}

#' Calculate expected survival.
#'
#' \code{expected_surv} calculates expected survival (area under survival curve).
#'
#' @param time vector of time to estimate probabilities.
#' @param surv vector of survival probabilities.
#' @return
#' expected survival.
#' @export
expected_surv <- function(time, surv) {
  require(zoo)
  sum(diff(time[order(time)])*rollmean(surv[order(time)],2))
}

#' Plot Markov trace from a partitioned survival model.
#'
#' \code{plot_trace_PSM} plots Markov trace from a partitioned survival model.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param partsurv.model partitioned survival model.
#' @param PA run probabilistic analysis.
#' Default = FALSE.
#' @return
#' a plot of the cohort trace.
#' @export
plot_trace_PSM <- function(time, partsurv.model, PA=F) {
  if (PA) {
    matplot(time, partsurv.model$Mean, type = 'l', lty = 1, ylab = "Markov trace")
    title(main = partsurv.model$chosen_models)
    matlines(time, partsurv.model$CI[,,1], lty = 2)
    matlines(time, partsurv.model$CI[,,2], lty = 2)
  } else {
    matplot(time, partsurv.model$trace, type = "l", lty = 1, ylab = "Markov trace")
    title(main = partsurv.model$chosen_models)
  }
}

#' Fit partitioned survival model on all combinations of chosen PFS and OS parametric survival functions.
#'
#' \code{all_partsurv} fits partitioned survival model on all combinations of chosen PFS and OS parametric survival functions.
#'
#' @param pfs_survHE survHE obj fitting PFS.
#' @param os_survHE survHE obj fitting OS.
#' @param choose_PFS preferred PFS distribution.
#' @param choose_OS preferred OS distribution.
#' @param time numeric vector of time to estimate probabilities.
#' @param PA run probabilistic analysis.
#' Default = FALSE.
#' @param n_sim number of PA simulations.
#' Default = 100.
#' @param seed seed for random number generation.
#' Default = 421.
#' @return
#' a list containing Markov trace, expected survival, survival probabilities, transition probabilities.
#' @export
all_partsurv <- function(pfs_survHE, os_survHE, choose_PFS, choose_OS, time = times, PA = FALSE, n_sim = 100, seed = 421) {
  all_m_M_PSM <- list() # empty list to store Markov traces
  model_names <- c(choose_PFS, choose_OS) # all parametric model names

  for (i in 1:length(model_names)) { # loop through model for PFS
    for (j in 1:length(model_names)) { # loop through model for OS
      # Construct a partitioned survival model out of the fitted models
      # fit a partitioned survival model for DRUG
      # fit partitioned survival model
      m_M_PSM <- partsurv(pfs_survHE = pfs_survHE,
                          os_survHE  = os_survHE,
                          choose_PFS = model_names[i],
                          choose_OS  = model_names[j],
                          time = time, PA = PA, n_sim = n_sim, seed = seed)
      # store model outputs
      all_m_M_PSM <- list.append(all_m_M_PSM, m_M_PSM)
    }
  }
  # Name the outputs
  model_names_combos <- expand.grid(model_names, model_names)
  names(all_m_M_PSM) <- paste0("PFS: ", model_names_combos$Var2, ", ",
                               "OS: ", model_names_combos$Var1)
  return(all_m_M_PSM)
}

#' Generate data for partitioned survival models (PFS and OS data).
#'
#' \code{gen_data} generates survival data for overall survival (OS) and progression-free survival (PFS).
#'
#' @param n_pat number of patients.
#' @param n_years follow-up period.
#' @return
#' generated survival data.
#' @export
gen_data <- function(n_pat, n_years){
  # specification of hazard functions to generate data from
  hazardf <- gems::generateHazardMatrix(n_s)
  colnames(hazardf@list.matrix) <-
    rownames(hazardf@list.matrix) <- v_n

  # specifying the transition hazard from healthy -> sick
  hazardf[["healthy","sick"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }

  # specifying the transition hazard from healthy -> dead
  hazardf[["healthy","dead"]] <- function (t, r1, r2){
    flexsurv::hgompertz(t,r1, r2)
  }

  # specifying the transition hazard from sick -> dead
  hazardf[["sick","dead"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }

  # list of parameters for the hazard functions defined above
  mu        <- gems::generateParameterMatrix(hazardf)
  rownames(mu@list.matrix) <-
    colnames(mu@list.matrix) <- v_n

  mu[["healthy", "sick"]] <- list(1.5, 6)      #  the Weibull parameters for H -> S
  mu[["healthy", "dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["sick",    "dead"]] <- list(0.5,4)       #  the Weibull parameters for S -> D

  # simulate the cohort
  cohort <- gems::simulateCohort(
    transitionFunctions = hazardf,
    parameters = mu,
    cohortSize = n_pat,
    to = n_years)

  # extract the simulated true data
  true_data <- cohort@time.to.state
  colnames(true_data) <- v_n

  true_data$dead[is.na(true_data$dead)] <- n_years
  true_data$sick[is.na(true_data$sick)] <- true_data$dead[is.na(true_data$sick)]


  # create a status variable that will capture the transition events
  true_status         <- matrix(NA, nrow = n_pat, ncol = n_s, dimnames = list(1:n_pat,v_n))
  true_status         <- as.data.frame(true_status)
  true_status$healthy <- ifelse(is.na(true_data$healthy),0,1)
  true_status$dead    <- ifelse(true_data$dead == n_years, 0, 1)
  true_status$sick    <- ifelse(true_data$dead == true_data$sick, 0, 1)


  censtime <- runif(n = n_pat, 0, n_years)

  censored_sick <- ifelse(censtime      <= true_data$sick |
                            true_data$sick >  5, 1, 0)
  censored_dead <- ifelse(censtime <= true_data$dead|
                            true_data$dead >5, 1, 0)

  sim_data <- true_data

  sim_data$sick[censored_sick == 1] <-  censtime[censored_sick == 1]
  sim_data$sick[sim_data$sick >5 ]  <-  5

  sim_data$dead[censored_dead == 1] <-  censtime[censored_dead == 1]
  sim_data$dead[sim_data$dead >5] <-  5

  status <- true_status

  status$sick[censored_sick == 1] = 0
  status$dead[censored_dead == 1] = 0

  # Usually trials report OS/PFS outcomes so we will recreate those

  OS_PFS_data <- data.frame(row.names = 1:n_pat)

  OS_PFS_data$PFS_time        <- apply(sim_data[, c("sick","dead")], 1, min)
  OS_PFS_data$PFS_status      <- ifelse(status$dead == 1 | status$sick == 1, 1, 0 )

  OS_PFS_data$OS_time         <- sim_data$dead
  OS_PFS_data$OS_status       <- status$dead
  list(cohort = cohort, true_data = true_data, true_status = true_status,
       sim_data =  sim_data, status = status, OS_PFS_data = OS_PFS_data)
}

#' Fit multi-state model.
#'
#' \code{fit.mstate} fits multi-state model.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param trans matrix of transition history.
#' @param data dataframe containing the time and status variables.
#' @param add add superimposed curves to the KM plot.
#' Default = FALSE.
#' @param extrapolate extrapolate beyond model time horizon.
#' Default = FALSE.
#' @param times time horizon the extrapolation is done over.
#' @return
#' Multi-state model fit.
#' @export
fit.mstate <- function(time, status, trans,  data = data , add = FALSE, extrapolate = FALSE, times) {
  data$time  <- data[, time  ]
  data$tatus <- data[, status]

  if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }

  # Progression free survival
  KM.fit     <-     survfit(Surv(time, status) ~ trans , data = data)                                 # fit Kaplan-Meier curve
  fit.llogis <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, status) ~ trans + sdlog(trans), data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "gamma"  ) # fit model with gamma distribution

  # extrapolate all models beyond the KM curve
  if(add){lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  if(!add){plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal"), col = 1:5, lty = rep(1, 5), bty="n")

  # compare AIC values
  AIC <- c(    Loglogistic = AIC(fit.llogis),
               Weibull     = AIC(fit.weib),
               Lognormal   = AIC(fit.lnorm),
               Gamma       = AIC(fit.gamma))
  AIC= round(AIC,3)

  # compare BIC values
  BIC <- c(    Loglogistic = BIC(fit.llogis),
               Weibull     = BIC(fit.weib),
               Lognormal   = BIC(fit.lnorm),
               Gamma       = BIC(fit.gamma))

  BIC <- round(BIC,3)

  res <- list( Loglogistic = fit.llogis,
               Weibull     = fit.weib,
               Lognormal   = fit.lnorm,
               Gamma       = fit.gamma,
               AIC         = AIC,
               BIC         = BIC)
  res
}

#' Compute Markov trace out of a multi-state model using DES.
#'
#' \code{trace.DES} computes Markov trace out of a multi-state model using DES.
#'
#' @param msm_sim multi-state model.
#' @param tmat matrix of transition history.
#' @param n_i number of individuals.
#' @param times time horizon the extrapolation is done over.
#' @return
#' Matrix of Markov trace
#' @export
trace.DES <- function(msm_sim = des_sim, tmat, n_i, times) {
  # Restructure the data to extract markov trace
  data.mstate.sim <- data.frame(cbind(matrix(t(msm_sim$st), ncol=1),
                                      matrix(t(msm_sim$t) , ncol=1)))
  colnames(data.mstate.sim) <- c("state","time")
  data.mstate.sim$subject <- rep(1:n_i, each = ncol(msm_sim$st))

  data.mstate.sim = na.omit(data.mstate.sim)
  data.mstate.sim = data.mstate.sim[!duplicated(data.mstate.sim), ] # remove duplicate entries in the dataset

  # create transition intensitiy matrix with initial values based on the structure of tmat
  twoway7.q               <- tmat
  twoway7.q[!is.na(tmat)] <- 0.5
  twoway7.q[is.na(tmat)]  <- 0
  # fit msm model only so that we can extract the prevalence (i.e. trace) thrrough the prevalence.msm function

  fit.msm.sim <- msm(state ~ time,subject = subject, data = data.mstate.sim, qmatrix = twoway7.q,
                     exacttimes = T, use.deriv = TRUE, analyticp = FALSE, fixedpars = TRUE, hessian = F)

  M.tr.des <- prevalence.msm(fit.msm.sim, times = times) # Markov trace when DES model is used

  return(M.tr.des[[3]]/100)
}

#' Plot cohort trace from microsimulation model
#'
#' \code{plot_trace_microsim} computes Markov trace out of a multi-state model using DES.
#'
#' @param m_M cohort trace array from microsimulation model.
#' @return
#' Plot of the cohort trace.
#' @export
plot_trace_microsim <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE))))
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)

}
