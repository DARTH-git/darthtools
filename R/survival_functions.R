#' Fit multiple survival models survival data
#'
#' \code{fit.fun} fits multiple survival models to survival data using survHE.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param data dataframe containing the time and status variables.
#' @param extrapolate extrapolate beyond model time horizon.
#' Default = FALSE.
#' @param times time horizon the extrapolation is done over.
#' @return
#' a list containing all survival model objects.
#' @export
fit.fun <- function(time, status, data = data, extrapolate = FALSE, times, k = 2) {
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
  mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "rps")
  # Run the models using MLE via flexsurv
  fit.survHE <- fit.models(formula = Surv(time, status) ~ 1, data = data, distr = mods, k = k)

  # Extrapolate all models beyond the KM curve and plot
  # KM.fit <- survfit(Surv(time, status) ~ 1, data = data) # fit Kaplan-Meier curve
  # plot(KM.fit, ylab = "Survival", xlab = "Time", ylim = c(0,1), xlim = c(0, max(plot.times)), conf.int = F, mark.time = T)
  # lines(fit.survHE$models$Exponential,      t = plot.times, col = 2, ci = F)
  # lines(fit.survHE$models$`Weibull (AFT)`,  t = plot.times, col = 3, ci = F)
  # lines(fit.survHE$models$Gamma,            t = plot.times, col = 4, ci = F)
  # lines(fit.survHE$models$`log-Normal`,     t = plot.times, col = 5, ci = F)
  # lines(fit.survHE$models$`log-Logistic`,   t = plot.times, col = 6, ci = F)
  # lines(fit.survHE$models$Gompertz,         t = plot.times, col = 7, ci = F)
  # lines(fit.survHE$models$`Royston-Parmar`, t = plot.times, col = 8, ci = F)
  # # add a legend
  # legend("topright", cex = 0.7, c("Kaplan-Meier", names(fit.survHE$models)), col = 1:(length(mods)+1), lty = rep(1, (length(mods)+1)), bty="n")
  print(plot(fit.survHE,
       t               = plot.times,
       add.km          = T,
       legend.position = c(0.75, 0.6),
       legend.text     = element_text(size = 9),
       legend.title    = element_text(size = 11)) +
       theme(plot.margin = unit(c(1,3.5,0.5,0.5), "cm")) + # top, right, bottom, left
       theme(legend.position = c(1.14,0.5)) +
       labs(title      = "Fitted survival curves"))

  # Compare AIC values
  AIC <- fit.survHE$model.fitting$aic
  AIC <- round(AIC,3)

  # Compare BIC values
  BIC <- fit.survHE$model.fitting$bic
  BIC <- round(BIC,3)

  names(AIC) <- names(BIC) <- names(fit.survHE$models)

  # Store and return results
  res <- list(model.objects = fit.survHE,
              AIC           = AIC,
              BIC           = BIC)
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
#' a list containing all survival model objects.
#' @export
fit.fun.cure <- function(time, status, data = data, extrapolate = FALSE, times) {
  require(flexsurvcure)
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

  # Extrapolate all models beyond the KM curve and plot - cure models
  KM.fit <- survfit(Surv(time, status) ~ 1, data = data) # fit Kaplan-Meier curve
  plot(KM.fit, ylab = "Survival", xlab = "Time", ylim = c(0,1), xlim = c(0, max(plot.times)), conf.int = F, mark.time = T)
  lines(fit.survcure$`Exponential Cure`,   t = plot.times, col = 2, ci = F)
  lines(fit.survcure$`Weibull (AFT) Cure`, t = plot.times, col = 3, ci = F)
  lines(fit.survcure$`Gamma Cure`,         t = plot.times, col = 4, ci = F)
  lines(fit.survcure$`log-Normal Cure`,    t = plot.times, col = 5, ci = F)
  lines(fit.survcure$`log-Logistic Cure`,  t = plot.times, col = 6, ci = F)
  lines(fit.survcure$`Gompertz Cure`,      t = plot.times, col = 7, ci = F)
  # add a legend
  legend("topright", cex = 0.7, c("Kaplan-Meier", names(fit.survcure)), col = 1:(length(mods)+1), lty = rep(1, (length(mods)+1)), bty="n")

  # Compare AIC values
  AIC <- unlist(lapply(fit.survcure, function(x) AIC(x)))  # cure models
  AIC <- round(AIC,3)

  # Compare BIC values
  BIC <- unlist(lapply(fit.survcure, function(x) BIC(x)))  # cure models
  BIC <- round(BIC,3)

  names(AIC) <- names(BIC) <- names(fit.survcure)

  # Store and return results
  res <- list(models = fit.survcure,
              AIC    = AIC,
              BIC    = BIC)
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
#' @param v_names_states vector of state names.
#' @param PA run probabilistic analysis.
#' Default = FALSE.
#' @param n_sim number of PA simulations.
#' Default = 100.
#' @param seed seed for random number generation.
#' Default = 421.
#' @return
#' a list containing Markov trace, expected survival, survival probabilities, transition probabilities.
#' @export
partsurv <- function(pfs_survHE, os_survHE, choose_PFS, choose_OS, time = times, v_names_states, PA = FALSE, n_sim = 100, seed = 421){
  set.seed(seed)
  deter <- ifelse(PA == 1, 0, 1) # determine if analysis is deterministic or probabilistic
  chosen_models <- paste0("PFS: ", choose_PFS, ", ", "OS: ", choose_OS) # chosen model names

  # Model-setup
  # model objects
  pfs_survHE <- pfs_survHE$model.objects
   os_survHE <-  os_survHE$model.objects
  # model names
  mod.pfs <- names(pfs_survHE$models)
   mod.os <- names(os_survHE$models)
  # chosen model index based on name
  mod.pfs.chosen <- which(mod.pfs == choose_PFS)
   mod.os.chosen <- which(mod.pfs == choose_OS)

  # Calculate survival probabilities
  if (deter == 0) { # probabilistic
    fit_PFS <- make.surv(pfs_survHE,
                         mod = mod.pfs.chosen,
                         nsim = n_sim,
                         t = times)
    fit_OS  <- make.surv(os_survHE,
                         mod = mod.os.chosen,
                         nsim = n_sim,
                         t = times)
    pfs.surv <- surv_prob(fit_PFS, PA = TRUE)
    os.surv  <- surv_prob( fit_OS, PA = TRUE)
  } else { # deterministic
    pfs.surv <- surv_prob(pfs_survHE$models[[which(mod.pfs == choose_PFS)]], time = times)
    os.surv  <- surv_prob( os_survHE$models[[which(mod.os  ==  choose_OS)]], time = times)
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
  Sick                 <- os.surv - pfs.surv    # estimate the probability of remaining in the progressed state
  check_PFS_OS(Sick)                            # print warning message if PFS > OS
  Sick[Sick < 0]       <- 0                     # in cases where the probability is negative replace with zero
  Healthy              <- pfs.surv              # probability of remaining stable
  Dead                 <- 1 - os.surv           # probability of being Dead
  trace <- abind(Healthy,
                 Sick,
                 Dead, rev.along = 0)

  # Calculate Markov trace
  if (deter == 0){ # probabilistic
    trace      <- aperm(trace, perm = c(1,3,2)) # Markov trace of all simulations
    mean.trace <- apply(trace, 1:2, mean)       # average Markov trace across simulations
    CI         <- apply(trace, 1:2, quantile, probs = c(0.025, 0.975)) # trace confidence intervals
    CI         <- aperm(CI, perm = c(2,3,1))
    dimnames(mean.trace)[[2]] <- v_names_states
    dimnames(CI)[[3]] <- c("low", "high")
    dimnames(CI)[[2]] <- v_names_states
    dimnames(trace)[[2]] <- v_names_states

    # List of items to return
    res <- list(trace = trace, CI = CI, Mean = mean.trace, # Markov trace
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  } else { # deterministic
    dimnames(trace)[[2]] <- v_names_states

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
surv_prob <- function(model, times = NULL, PA = FALSE) {
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
                          time = time, PA = PA, v_names_states = v_names_states, n_sim = n_sim, seed = seed)
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
#' @param n_years follow-up period in years.
#' @return
#' generated survival data.
#' @export
gen_data <- function(n_pat, n_years){
  # specification of hazard functions to generate data from
  hazardf <- gems::generateHazardMatrix(n_states)
  colnames(hazardf@list.matrix) <-
    rownames(hazardf@list.matrix) <- v_names_states

  # specifying the transition hazard from Healthy -> Sick
  hazardf[["Healthy","Sick"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }

  # specifying the transition hazard from Healthy -> Dead
  hazardf[["Healthy","Dead"]] <- function (t, r1, r2){
    flexsurv::hgompertz(t,r1, r2)
  }

  # specifying the transition hazard from Sick -> Dead
  hazardf[["Sick","Dead"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }

  # list of parameters for the hazard functions defined above
  mu        <- gems::generateParameterMatrix(hazardf)
  rownames(mu@list.matrix) <-
    colnames(mu@list.matrix) <- v_names_states

  mu[["Healthy", "Sick"]] <- list(1.5, 6)      # the Weibull parameters for H -> S
  mu[["Healthy", "Dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["Sick",    "Dead"]] <- list(0.5,4)       # the Weibull parameters for S -> D

  # simulate the cohort
  cohort <- gems::simulateCohort(
    transitionFunctions = hazardf,
    parameters = mu,
    cohortSize = n_pat,
    to = n_years)

  # extract the simulated true data
  true_data <- cohort@time.to.state
  colnames(true_data) <- v_names_states

  true_data$Dead[is.na(true_data$Dead)] <- n_years
  true_data$Sick[is.na(true_data$Sick)] <- true_data$Dead[is.na(true_data$Sick)]

  # create a status variable that will capture the transition events
  true_status         <- matrix(NA, nrow = n_pat, ncol = n_states, dimnames = list(1:n_pat,v_names_states))
  true_status         <- as.data.frame(true_status)
  true_status$Healthy <- ifelse(is.na(true_data$Healthy),0,1)
  true_status$Dead    <- ifelse(true_data$Dead == n_years, 0, 1)
  true_status$Sick    <- ifelse(true_data$Dead == true_data$Sick, 0, 1)

  censtime <- runif(n = n_pat, 0, n_years)

  censored_Sick <- ifelse(censtime <= true_data$Sick |
                            true_data$Sick >  5, 1, 0)
  censored_Dead <- ifelse(censtime <= true_data$Dead|
                            true_data$Dead >5, 1, 0)

  sim_data <- true_data

  sim_data$Sick[censored_Sick == 1] <- censtime[censored_Sick == 1]
  sim_data$Sick[sim_data$Sick > 5 ] <- 5

  sim_data$Dead[censored_Dead == 1] <- censtime[censored_Dead == 1]
  sim_data$Dead[sim_data$Dead > 5] <- 5

  status <- true_status

  status$Sick[censored_Sick == 1] <- 0
  status$Dead[censored_Dead == 1] <- 0

  # Usually trials report OS/PFS outcomes so we will recreate those

  OS_PFS_data <- data.frame(row.names = 1:n_pat)

  OS_PFS_data$PFS_time   <- apply(sim_data[, c("Sick","Dead")], 1, min)
  OS_PFS_data$PFS_status <- ifelse(status$Dead == 1 | status$Sick == 1, 1, 0 )

  OS_PFS_data$OS_time    <- sim_data$Dead
  OS_PFS_data$OS_status  <- status$Dead
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

#' Converts cumulative hazards to hazard rates
#'
#' \code{cumhaz_to_haz} converts cumulative hazards to hazard rates across time points.
#'
#' @param cumhaz vector of cumulative hazards
#' @return
#' vector of hazard rates
#' @export
cumhaz_to_haz <- function(cumhaz) {
  c(cumhaz[1], diff(cumhaz))
}

#' Calculates hazards
#' @export
hazard.fn = function( x, t, start, ...) {
  ret <- x$dfns$h(t, ...) * (1 - x$dfns$p(start, ...))
  ret[t < start] <- 0
  ret
}

#' Bootstrap hazards
#' @export
normboot.haz <- function (x, t, start, newdata = NULL, X = NULL, fn, B) {
  sim <- normboot.flexsurvreg(x, B, newdata = newdata, X = X)
  X <- attr(sim, "X")
  if (!is.list(sim))
    sim <- list(sim)
  ret <- array(NA_real_, dim = c(nrow(X), B, length(t)))
  for (k in 1:nrow(X)) {
    for (i in seq(length = B)) {
      fncall <- list(t, start)
      fncall[["x"]] = x

      for (j in seq(along = x$dlist$pars)) fncall[[x$dlist$pars[j]]] <- sim[[k]][i,
                                                                                 j]
      for (j in seq_along(x$aux)) fncall[[names(x$aux)[j]]] <- x$aux[[j]]
      ret[k, i, ] <- do.call(fn, fncall)
    }
  }
  if (nrow(X) == 1)
    ret[1, , , drop = FALSE]
  else ret
}

#' Bootstrap hazards
#' @export
boot.haz <- function (x, t, start = 0 ,X = NULL, newdata =NULL, B = 1000) {
  dat <- x$data
  Xraw <- model.frame(x)[, unique(attr(model.frame(x), "covnames.orig")),
                         drop = FALSE]
  isfac <- sapply(Xraw, function(x) {
    is.factor(x) || is.character(x)
  })

  if (is.null(newdata)) {
    if (is.vector(X))
      X <- matrix(X, nrow = 1)
    if (x$ncovs > 0 && is.null(X)) {
      if (!all(isfac)) {
        nd <- colMeans(model.matrix(x))
        X <- matrix(nd, nrow = 1, dimnames = list(NULL,
                                                  names(nd)))
        attr(X, "newdata") <- as.data.frame(X)
      }
      else {
        X <- unique(model.matrix(x))
        nam <- as.matrix(unique(Xraw))
        for (i in 1:ncol(nam)) nam[, i] <- paste(colnames(nam)[i],
                                                 nam[, i], sep = "=")
        rownames(X) <- apply(nam, 1, paste, collapse = ",")
        attr(X, "newdata") <- unique(Xraw)
      }
    }
    else if (is.null(X))
      X <- as.matrix(0, nrow = 1, ncol = max(x$ncoveffs,
                                             1))
    else if (!is.matrix(X) || (is.matrix(X) && ncol(X) !=
                               x$ncoveffs)) {
      plural <- if (x$ncoveffs > 1)
        "s"
      else ""
      stop("expected X to be a matrix with ", x$ncoveffs,
           " column", plural, " or a vector with ", x$ncoveffs,
           " element", plural)
    }else {
      attr(X, "newdata") <- X
      colnames(attr(X, "newdata")) <- colnames(model.matrix(x))
    }
  }

  hazard.fn <- flexsurv:::expand.summfn.args(hazard.fn)
  fncall <- list( t, start)
  beta <- if (x$ncovs == 0) {
    0
  }else x$res[x$covpars, "est"]
  dlist <- x$dlist
  ret <- vector(nrow(X), mode = "list")
  if (!is.null(newdata)) {
    nd <- attr(X, "newdata")
    covnames <- apply(as.data.frame(nd), 1, function(x) paste0(names(nd),
                                                               "=", x, collapse = ", "))
  }else covnames <- rownames(X)
  names(ret) <- covnames
  for (i in 1:nrow(X)) {
    basepars.mat <- flexsurv:::add.covs(x, x$res.t[dlist$pars, "est"],
                                        beta, X[i, , drop = FALSE], transform = FALSE)
    basepars <- as.list(as.data.frame(basepars.mat))
    fncall[dlist$pars] <- basepars
    for (j in seq_along(x$aux)) {
      fncall[[names(x$aux)[j]]] <- x$aux[[j]]
    }

    fncall[["x"]] <- x
    # browser()

    y <- do.call(hazard.fn, fncall)
    ret <- normboot.haz(x = x, t = t, start = start,
                        X = X, fn = hazard.fn, B = B)

  }
  ret
}

#' Bootstrap hazards ratios
#' @export
boot_hr <- function(surv_model1, surv_model2, times, B = 100){
  # bootstrap hazards of survival model 1
  boot.haz1 = boot.haz(x = surv_model1, B = 100, t = times)
  boot.haz1[boot.haz1==0] = 0.01
  boot.haz1 <- boot.haz1[1,,]
  # take log hazards
  boot.log.haz1 <- log(boot.haz1)

  # bootstrap hazards of survival model 2
  boot.haz2 = boot.haz(x = surv_model2, B = 100, t = times)
  boot.haz2[boot.haz2==0] = 0.01
  boot.haz2 <- boot.haz2[1,,]
  # take log hazards
  boot.log.haz2 <- log(boot.haz2)

  # calculate log hazard ratios of model 1 vs. 2
  boot.log.hr <- boot.log.haz1 - boot.log.haz2
  # return median hazard ratios with 95% confidence intervals
  HR <- apply(boot.log.hr, 2, function(x)exp(quantile(x, probs = c(0.025, 0.5, 0.975))))
  HR <- as.data.frame(t(HR))
  HR$times <- times
  colnames(HR) <- c("lcl", "med", "ucl", "time")
  return(HR)
}

#' Print a warning message if PFS > OS
#'
#' \code{check_PFS_OS} prints a warning message if PFS > OS.
#'
#' @param Sick vector (or matrix) of PFS - OS probabilities
#' @return
#' a warning message if PFS > OS
#' @export
check_PFS_OS <- function(Sick){
  Sick1 <- as.matrix(Sick)
  for (j in 1:ncol(Sick1)) {
    if (length(which(Sick1[,j] < 0) > 0)) {
      print(paste0("Warning: PFS > OS starting cycle ", min(which(Sick1[,j] < 0)), " at simulation ", j))
    }
  }
}

