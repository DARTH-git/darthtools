#' Fit multiple survival models on survival data
#'
#' \code{fit.fun} fits multiple survival models to survival data using survHE.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param covariate logical value indicating whether treatment is being used as a covariate in parameteric survival models.
#' Default = FALSE.
#' @param rx numerical value indicating the treatment variable used as as covariate in parameteric survival models.
#' @param data dataframe containing the time and status variables.
#' @param extrapolate extrapolate beyond model time horizon.
#' Default = FALSE.
#' @param times time horizon the extrapolation is done over.
#' @param k number of knots in Royston-Parmar spline model.
#' Default = 2.
#' @param legend_position position of the legend.
#' Default = "bottom".
#' @param xlow time horizon the extrapolation is done over.
#' Default = min(time).
#' @param xhigh time horizon the extrapolation is done over.
#' Default = max(time).
#' @param ylow time horizon the extrapolation is done over.
#' Default = 0.
#' @param yhigh time horizon the extrapolation is done over.
#' Default = 1.
#' @param risktable time horizon the extrapolation is done over.
#' Default = F.
#' @param mods a vector of models to fit.
#' Choose from = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "rps", "gengamma").
#' @return
#' a list containing all survival model objects.
#' @export
fit.fun <- function(time, status, covariate = F, rx = NULL, data, extrapolate = FALSE, times, k = 2,
                    legend_position = "bottom", xlow = min(times), xhigh = max(times), ylow = 0, yhigh = 1, risktable = F,
                    mods = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "rps", "gengamma")) {
  require(survHE)
  require(survminer)
  require(tm)

  ## Set up
  # Extract the right data columns
  data$time   <- data[,   time]
  data$status <- data[, status]

  # Specify time horizon based on whether to extrapolate
  if (extrapolate == TRUE)  {
    plot.times <- times
  } else if  (extrapolate == FALSE) {
    plot.times <- seq(min(times),max(data$time), by = diff(times)[1])
  }

  ## Fit parametric survival models
  # Define the vector of models to be used
  mods <- mods
  # Run the models using MLE via flexsurv
  if (covariate) {
    fit.survHE <- fit.models(formula = Surv(time, status) ~ rx, data = data, distr = mods, k = k)
  } else {
    fit.survHE <- fit.models(formula = Surv(time, status) ~ 1, data = data, distr = mods, k = k)
  }

  ## Plots
  if (covariate) { # fit Kaplan-Meier curve
    KM.fit <- survfit(Surv(time, status) ~ rx, data = data)
  } else {
    KM.fit <- survfit(Surv(time, status) ~ 1, data = data)
  }

  S_superimpose <- ggsurvplot(
    KM.fit,
    data       = data,           # specify dataset
    size       = 1,              # change line size
    palette = "Dark2",
    conf.int   = TRUE,           # add confidence interval
    pval       = F,              # add p-value
    risk.table = risktable,      # add risk table
    risk.table.height = 0.25,    # useful to change when you have multiple groups
    ggtheme    = theme_bw(),     # change ggplot2 theme
    xlab       = 'Time',         # change x-axis label
    ylab       = 'Survival',     # change y-axis label
    xlim       = c(xlow, xhigh),
    ylim       = c(ylow, yhigh),
    title      = "Fitted survival curves", # add title
    legend     = legend_position, # change legend position
    legend.title = "KM",
  )
  surv_probs <- list()

  mod.aesthetics <- data.frame(mod = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "rps", "gengamma"),
                               color = c("#DB8681", "#B69B34", "#68A232", "#35A58C", "#30A8C3", "#249de5", "#9F95D4", "#E271C9"),
                               label = c("b", "c", "d", "e", "f", "g", "h", "i"))

  fit.fun.colors <- mod.aesthetics %>% filter(mod %in% mods) %>% select(color) %>% pull()
  fit.fun.labels <- mod.aesthetics %>% filter(mod %in% mods) %>% select(label) %>% pull()

  for (i in 1:length(fit.survHE$models)) {
    model_output <- fit.survHE$models[[i]]
    # extract the estimated survival probabilities and the confidence intervals
    surv_probs[[i]] <- as.data.frame(summary(model_output, t = plot.times))
    surv_probs[[i]]$Model <- paste0(fit.fun.labels[i], names(fit.survHE$models)[i])
    # superimpose the survival probabilities
  }
  surv_probs       <- rbindlist(surv_probs)
  surv_probs$Model <- as.factor(surv_probs$Model)

  if (covariate) { # if using treatment as a covariate
    surv_probs1 <- as.data.frame(surv_probs)
    colnames(surv_probs1)[1] <- "time"
    surv_probs1 <- surv_probs1[, colnames(surv_probs1) %in%
                                 c("time",
                                   colnames(surv_probs1)[grep("est", colnames(surv_probs1))],
                                   "Model")]
    surv_probs1 <- surv_probs1 %>% pivot_longer(cols = !c(Model, time), names_to = "rx",
                                                values_to = "est") %>% as.data.frame()
    # surv_probs1$rx <- removeWords(surv_probs1$rx, c("rx.", ".est"))
    surv_probs1$rx <- as.factor(surv_probs1$rx)
    surv_probs <- surv_probs1

    filler_colour <- rep("white", length(unique(surv_probs1$rx)))
    filler_label  <- rep("", length(unique(surv_probs1$rx)))

    S_superimpose$plot  <- S_superimpose$plot +
      geom_line(aes(x=time, y=est, color = Model,  linetype = rx),
                size = 0.75, alpha = 1, key_glyph = "path",
                data=surv_probs) +
      scale_color_manual(name = "", values = c(filler_colour, fit.fun.colors),
                         labels = c(filler_label, names(fit.survHE$models))) +
      scale_linetype_manual(name = "Strategy", values = 1:length(unique(surv_probs1$rx))) +
      scale_x_continuous(limits = c(0, max(plot.times))) +
      coord_cartesian(xlim = c(xlow, xhigh), ylim=c(ylow, yhigh)) +
      guides(color=guide_legend(override.aes = list(size=1.2)))
    print(S_superimpose)

  } else { # if not using treatment as a covariate
    S_superimpose$plot  <- S_superimpose$plot +
      geom_line(aes(x=time, y=est, color = Model),
                size = 0.75, alpha = 1, key_glyph = "path",
                data=surv_probs) +
      scale_color_manual(name = "Models", values = c("#988791", fit.fun.colors),
                         labels = c("Kaplan-Meier", names(fit.survHE$models))) +
      scale_x_continuous(limits = c(0, max(plot.times))) +
      coord_cartesian(xlim = c(xlow, xhigh), ylim=c(ylow, yhigh)) +
      guides(color=guide_legend(override.aes = list(size=1.2)))
    print(S_superimpose)
  }

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

#' Fit multiple mixture cure models on survival data
#'
#' \code{fit.fun.cure} fits multiple mixture cure models to survival data using flexsurv.
#'
#' @param time numeric vector of time to estimate probabilities.
#' @param status numeric vector of event status.
#' @param data dataframe containing the time and status variables.
#' @param covariate logical value indicating whether treatment is being used as a covariate in parameteric survival models.
#' Default = FALSE.
#' @param rx numerical value indicating the treatment variable used as as covariate in parameteric survival models.
#' @param extrapolate extrapolate beyond model time horizon.
#' Default = FALSE.
#' @param times time horizon the extrapolation is done over.
#' @param legend_position position of the legend.
#' Default = "top".
#' @param xlow time horizon the extrapolation is done over.
#' Default = min(time).
#' @param xhigh time horizon the extrapolation is done over.
#' Default = max(time).
#' @param ylow time horizon the extrapolation is done over.
#' Default = 0.
#' @param yhigh time horizon the extrapolation is done over.
#' Default = 1.
#' @param risktable time horizon the extrapolation is done over.
#' Default = F.
#' @param mods a vector of models to fit.
#' Default = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "gengamma").
#' @return
#' a list containing all survival model objects.
#' @export
fit.fun.cure <- function(time, status, covariate = F, rx = "rx", data = data, extrapolate = FALSE, times,
                         legend_position = "bottom", xlow = min(times), xhigh = max(times), ylow = 0, yhigh = 1, risktable = F,
                         mods = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "gengamma")) {
  require(flexsurvcure)
  require(survHE)
  require(survminer)
  # require(tm)

  ## Set up
  # Extract the right data columns
  data$time   <- data[,   time]
  data$status <- data[, status]

  # Specify time horizon based on whether to extrapolate
  if (extrapolate == TRUE)  {
    plot.times <- times
  } else if  (extrapolate == FALSE) {
    plot.times <- seq(min(times),max(data$time), by = diff(times)[1])
  }

  ## Fit parametric survival models
  # Define the vector of models to be used
  mods <- mods
  # Run the models using MLE via flexsurvcure
  if (covariate) {
    fit.survHE <- fit.models.cure(formula = Surv(time, status) ~ rx, data = data, distr = mods)
  } else {
    fit.survHE <- fit.models.cure(formula = Surv(time, status) ~ 1, data = data, distr = mods)
  }

  names(fit.survHE$models) <- paste0(names(fit.survHE$models), " Cure")

  ## Plots
  if (covariate) { # fit Kaplan-Meier curve
    KM.fit <- survfit(Surv(time, status) ~ rx, data = data)
  } else {
    KM.fit <- survfit(Surv(time, status) ~ 1, data = data)
  }

  S_superimpose <- ggsurvplot(
    KM.fit,
    data       = data,           # specify dataset
    size       = 1,              # change line size
    palette = "Dark2",
    conf.int   = TRUE,           # add confidence interval
    pval       = F,              # add p-value
    risk.table = risktable,      # add risk table
    risk.table.height = 0.25,    # useful to change when you have multiple groups
    ggtheme    = theme_bw(),     # change ggplot2 theme
    xlab       = 'Time',         # change x-axis label
    ylab       = 'Survival',     # change y-axis label
    xlim       = c(xlow, xhigh),
    ylim       = c(ylow, yhigh),
    title      = "Fitted survival curves", # add title
    legend     = legend_position, # change legend position
    legend.title = "KM",
  )
  surv_probs <- list()

  mod.aesthetics <- data.frame(mod = c("exp", "weibull", "gamma", "lnorm", "llogis", "gompertz", "rps", "gengamma"),
                               color = c("#DB8681", "#B69B34", "#68A232", "#35A58C", "#30A8C3", "#249de5", "#9F95D4", "#E271C9"),
                               label = c("b", "c", "d", "e", "f", "g", "h", "i"))

  fit.fun.colors <- mod.aesthetics %>% filter(mod %in% mods) %>% select(color) %>% pull()
  fit.fun.labels <- mod.aesthetics %>% filter(mod %in% mods) %>% select(label) %>% pull()

  for (i in 1:length(fit.survHE$models)) {
    model_output <- fit.survHE$models[[i]]
    # extract the estimated survival probabilities and the confidence intervals
    surv_probs[[i]] <- as.data.frame(summary(model_output, t = plot.times))
    surv_probs[[i]]$Model <- paste0(fit.fun.labels[i], names(fit.survHE$models)[i])
    # superimpose the survival probabilities
  }
  surv_probs       <- rbindlist(surv_probs)
  surv_probs$Model <- as.factor(surv_probs$Model)

  if (covariate) { # if using treatment as a covariate
    surv_probs1 <- as.data.frame(surv_probs)
    colnames(surv_probs1)[1] <- "time"
    surv_probs1 <- surv_probs1[, colnames(surv_probs1) %in%
                                 c("time",
                                   colnames(surv_probs1)[grep("est", colnames(surv_probs1))],
                                   "Model")]
    surv_probs1 <- surv_probs1 %>% pivot_longer(cols = !c(Model, time), names_to = "rx",
                                                values_to = "est") %>% as.data.frame()
    # surv_probs1$rx <- removeWords(surv_probs1$rx, c("rx.", ".est"))
    surv_probs1$rx <- as.factor(surv_probs1$rx)
    surv_probs <- surv_probs1

    filler_colour <- rep("white", length(unique(surv_probs1$rx)))
    filler_label  <- rep("", length(unique(surv_probs1$rx)))

    S_superimpose$plot  <- S_superimpose$plot +
      geom_line(aes(x=time, y=est, color = Model,  linetype = rx),
                size = 0.75, alpha = 1, key_glyph = "path",
                data=surv_probs) +
      scale_color_manual(name = "", values = c(filler_colour, fit.fun.colors),
                         labels = c(filler_label, names(fit.survHE$models))) +
      scale_linetype_manual(name = "Strategy", values = 1:length(unique(surv_probs1$rx))) +
      scale_x_continuous(limits = c(0, max(plot.times))) +
      coord_cartesian(xlim = c(xlow, xhigh), ylim=c(ylow, yhigh)) +
      guides(color=guide_legend(override.aes = list(size=1.2)))
    print(S_superimpose)

  } else { # if not using treatment as a covariate
    S_superimpose$plot  <- S_superimpose$plot +
      geom_line(aes(x=time, y=est, color = Model),
                size = 0.75, alpha = 1, key_glyph = "path",
                data=surv_probs) +
      scale_color_manual(name = "Models", values = c("#988791", fit.fun.colors),
                         labels = c("Kaplan-Meier", names(fit.survHE$models))) +
      scale_x_continuous(limits = c(0, max(plot.times))) +
      coord_cartesian(xlim = c(xlow, xhigh), ylim=c(ylow, yhigh)) +
      guides(color=guide_legend(override.aes = list(size=1.2)))
    print(S_superimpose)
  }

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

#' Fit partitioned survival model to survival data
#'
#' \code{partsurv} fits partitioned survival model to survival data.
#'
#' @param pfs_survHE survHE obj fitting PFS.
#' @param os_survHE survHE obj fitting OS.
#' @param l_d.data list of mean parameter estimates (for PFS and OS)
#' @param l_vc.data list of variance-covariance matrices (or their Cholesky decomposition) of parameter estimates (for PFS and OS)
#' @param par use parameter mean estimates and variance-covariance matrix instead of IPD
#' Default = FALSE
#' @param chol if the Cholesky decomposition of the variance-covariance matrix is used instead of the variance-covariance matrix
#' @param choose_PFS preferred PFS distribution.
#' @param choose_OS preferred OS distribution.
#' @param dist_name_PFS specify PFS distribution when par = T.
#' Choose from: Exponential, Weibull (AFT), Gamma, log-Normal, log-Logistic, Gompertz, Exponential Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure, log-Logistic Cure, Gompertz Cure.
#' @param dist_name_OS specify OS distribution when par = T.
#' Choose from: Exponential, Weibull (AFT), Gamma, log-Normal, log-Logistic, Gompertz, Exponential Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure, log-Logistic Cure, Gompertz Cure.
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
partsurv <- function(pfs_survHE = NULL, os_survHE = NULL, l_d.data = NULL, l_vc.data = NULL, par = FALSE, chol = FALSE,
                     dist_name_PFS = NULL, dist_name_OS = NULL, choose_PFS = NULL, choose_OS = NULL,
                     time = times, v_names_states, PA = FALSE, n_sim = 100, seed = 421){
  set.seed(seed)
  deter <- ifelse(PA == 1, 0, 1) # determine if analysis is deterministic or probabilistic

  if (par == TRUE) {
    dist_PFS <- dist_name_PFS
    dist_OS  <- dist_name_OS
    # deal with Cholesky decomposition

    if (chol == TRUE) {
      for (i in 1:length(l_vc.data)) {
        l_vc.data[[i]] <- crossprod(l_vc.data[[i]])
      }
    }
  } else {
    dist_PFS <- choose_PFS
    dist_OS  <- choose_OS
  }

  chosen_models <- paste0("PFS: ", dist_PFS, ", ", "OS: ", dist_OS) # chosen model names

  # Calculate survival probabilities
  if (deter == 0) { # probabilistic
    if (par == TRUE) { # if choose to use parameter mean estimates and variance-covariance matrix instead of IPD
      # randomly draw parameter values from multivariate normal distribution
      param_draws_PFS <- model.rmvnorm(dist.v  = dist_PFS,
                                       d.data  = l_d.data$PFS,
                                       vc.data = l_vc.data$PFS,
                                       n_sim   = n_sim)
      param_draws_OS  <- model.rmvnorm(dist.v  = dist_OS,
                                       d.data  = l_d.data$OS,
                                       vc.data = l_vc.data$OS,
                                       n_sim   = n_sim)
      # obtain survival probabilities
      pfs.surv <- os.surv <- matrix(NA, nrow = length(time), ncol = n_sim)
      for (j in 1:n_sim) {
        pfs.surv[, j] <- model.dist(dist.v = dist_PFS, d.data = param_draws_PFS[j, ], t = time)
        os.surv [, j] <- model.dist(dist.v = dist_OS,  d.data = param_draws_OS[j, ],  t = time)
      }
    } else {
      # Model-setup
      # model objects
      pfs_survHE <- pfs_survHE$model.objects
      os_survHE <-  os_survHE$model.objects
      # model names
      mod.pfs <- names(pfs_survHE$models)
      mod.os <- names(os_survHE$models)
      # chosen model index based on name
      mod.pfs.chosen <- which(mod.pfs == dist_PFS)
      mod.os.chosen <- which(mod.os == dist_OS)
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
    }
  } else { # deterministic
    if (par == TRUE) { # if choose to use parameter mean estimates and variance-covariance matrix instead of IPD
      # randomly draw parameter values from multivariate normal distribution
      param_draws_PFS <- model.rmvnorm(dist.v  = dist_PFS,
                                       d.data  = l_d.data[[1]],
                                       vc.data = l_vc.data[[1]],
                                       n_sim   = 1)
      param_draws_OS  <- model.rmvnorm(dist.v  = dist_OS,
                                       d.data  = l_d.data[[2]],
                                       vc.data = l_vc.data[[2]],
                                       n_sim   = 1)
      # obtain survival probabilities
      pfs.surv <- model.dist(dist.v = dist_PFS, d.data = param_draws_PFS[1, ], t = time)
      os.surv  <- model.dist(dist.v = dist_OS,  d.data =  param_draws_OS[1, ], t = time)
    } else {
      pfs.surv <- surv_prob(pfs_survHE$models[[which(mod.pfs == dist_PFS)]], time = times)
      os.surv  <- surv_prob( os_survHE$models[[which(mod.os  ==  dist_OS)]], time = times)
    }
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
    pfs.trans.prob <- apply(pfs.surv, 2, trans_prob)
    os.trans.prob  <- apply( os.surv, 2, trans_prob)
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
#' @param times time horizon to calculate survival probabilities.
#' @param PA run probabilistic analysis.
#' @param rx determines which treatment arm (same order as factor levels of treatment variable).
#' Default = FALSE.
#' @return
#' vector of survival probabilities.
#' @export
surv_prob <- function(model, times = NULL, PA = FALSE, rx = 1) {
  if (PA) {
    surv <- model$mat[[rx]][,-1]
  } else {
    surv <- summary(model, t = times, ci = F)[[rx]]$est
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

#' Bootstrap hazards ratios of two survival models
#' \code{boot_hr} computes bootstrap hazard ratios (HR) of two survival models (HR of model 1 vs. model 2)
#'
#' @param surv_model1 first survival model.
#' @param surv_model2 second survival model.
#' @param B number of bootstrap samples.
#' @param times time horizon the extrapolation of the survival model is done over.
#' @return
#' dataframe of hazard ratio statistics (2.5% percnetile, median, 97.5% percentile, time points)
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

#' Determine which rows the upper and lower values of each interval are (in the survival data set).
#'
#' \code{find_interval_limits} determines which rows the upper and lower values of each interval are (in the survival data set).
#'
#' @param start_time start time.
#' @param start_time survival times.
#' @return
#' matrix of limits.
#' @export
find_interval_limits <- function(start_time,
                                 surv_time){

  if (max(surv_time) > max(start_time))
    stop("Intervals must span all survival times. Later interval missing.")

  interval <- 1
  end_time <- start_time[2]
  upper <- NULL

  for (i in seq_along(surv_time)){
    if (surv_time[i] >= end_time) {
      upper[interval] <- i - 1
      interval <- interval + 1
      end_time <- start_time[interval + 1]
    }
  }

  cbind(lower = c(1, upper + 1),
        upper = c(upper, length(surv_time)))
}

#' Function to create at-risk table.
#'
#' \code{create_at_risk_table} creates at-risk table.
#'
#' @param survival_time vector of survival times.
#' @param Years vector of time intervals.
#' @param AtRisk vector of number of patients at risk.
#' @return
#' at-risk table (dataframe).
#' @export
create_at_risk_table <- function(survival_time, Years, AtRisk) {
  require(dplyr)
  time_interval <- c(Years, max(Years)+diff(Years)[1])
  df_nrisk <- as.data.frame(find_interval_limits(start_time = time_interval,
                                                 surv_time = survival_time))
  colnames(df_nrisk) <- c("Lower", "Upper")
  df_nrisk$Interval <- 1:nrow(df_nrisk)
  df_nrisk$Years <- Years[1:nrow(df_nrisk)]
  df_nrisk$AtRisk <- AtRisk[1:nrow(df_nrisk)]
  df_nrisk <- df_nrisk %>% dplyr::select(Interval, Years, Lower, Upper, AtRisk)
  return(df_nrisk)
}

#' Fit parametric mixture cure survival models for health economic evaluations.
#'
#' \code{fit.models.cure} fits parametric mixture cure survival models for health economic evaluations.
#'
#' @param formula a formula specifying the model to be used, in the form Surv(time,event)~treatment[+covariates] for flexsurv.
#' @param data A data frame containing the data to be used for the analysis. This must contain data for the 'event' variable. In case there is no censoring, then event is a column of 1s.
#' @param distr a (vector of) string(s) containing the name(s) of the model(s) to be fitted. Available options are: flexsurv: "exponential","gamma","genf","gengamma","gompertz","weibull", "weibullPH","loglogistic","lognormal" INLA: "exponential","weibull","lognormal","loglogistic" hmc: "exponential","gamma","genf","gengamma","gompertz","weibull","weibullPH", "loglogistic","lognormal".
#' @return
#' A model object.
#' @export
fit.models.cure <- function (formula = NULL, data, distr = NULL, method = "mle",
                             ...)
{
  exArgs <- list(...)
  exArgs$formula <- formula
  exArgs$data = data
  if (is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  method <- tolower(method)
  if (!method %in% c("mle")) {
    stop("Methods available for use are 'mle'")
  }
  survHE:::check_distributions(method, distr)
  if (method == "mle") {
    res <- survHE:::format_output_fit.models(lapply(distr, function(x) runMLE.cure(x,
                                                                                   exArgs)), method, distr, formula, data)
  }
  return(res)
}

#' run MLE on cure models
#' @export
runMLE.cure <- function (x, exArgs)
{
  formula <- exArgs$formula
  data = exArgs$data
  availables <- survHE:::load_availables()
  d3 <- survHE:::manipulate_distributions(x)$distr3
  x <- survHE:::manipulate_distributions(x)$distr
  tic <- proc.time()
  if (x == "survspline") {
    if (exists("bhazard", where = exArgs)) {
      bhazard <- exArgs$bhazard
    }
    else {
      bhazard <- NULL
    }
    if (exists("weights", where = exArgs)) {
      weights <- exArgs$weights
    }
    else {
      weights <- NULL
    }
    if (exists("subset", where = exArgs)) {
      subset <- exArgs$subset
    }
    else {
      subset <- NULL
    }
    if (exists("knots", where = exArgs)) {
      knots <- exArgs$knots
    }
    else {
      knots <- NULL
    }
    if (exists("k", where = exArgs)) {
      k <- exArgs$k
    }
    else {
      k <- 0
    }
    if (exists("bknots", where = exArgs)) {
      bknots <- exArgs$bknots
    }
    else {
      bknots <- NULL
    }
    if (exists("scale", where = exArgs)) {
      scale <- exArgs$scale
    }
    else {
      scale <- "hazard"
    }
    if (exists("timescale", where = exArgs)) {
      timescale <- exArgs$scale
    }
    else {
      timescale <- "log"
    }
    model <- flexsurv::flexsurvspline(formula = formula,
                                      data = data, k = k, knots = knots, bknots = bknots,
                                      scale = scale, timescale = timescale)
  }
  else {
    model <- flexsurvcure(formula = formula, data = data, dist = x, mixture = T)
  }
  toc <- proc.time() - tic
  model_name <- d3
  list(model = model, aic = model$AIC, bic = -2 * model$loglik +
         model$npars * log(model$N), dic = NULL, time2run = toc[3],
       model_name = model_name)
}

#' Randomly draw parameter values of survival models from multivariate normal distribution.
#'
#' \code{model.rmvnorm} randomly draws parameter values of survival models from multivariate normal distributio.
#'
#' @param dist.v a character string specifying the name of the survival model.
#' @param d.data a vector of mean parameter estimates of the survival model.
#' @param vc.data variance-covariance matrix (a matrix) of parameter estimates of the survival model.
#' @param n_sim number of random samples to draw.
#' Default = 100.
#' @param seed seed for random number generation.
#' Default = 421.
#' @return
#' A matrix of drawn parameter values.
#' @export
model.rmvnorm <- function(dist.v, d.data, vc.data, n_sim, seed = 421) {

  set.seed(seed) #

  if (!dist.v %in% c("Exponential", "Weibull (AFT)", "Gamma", "log-Normal",
                     "log-Logistic", "Gompertz", "Exponential Cure", "Weibull (AFT) Cure", "Gamma Cure", "log-Normal Cure",
                     "log-Logistic Cure", "Gompertz Cure")) {
    return(paste0("Incorrect distribution name, select from: Exponential, Weibull (AFT), Gamma, log-Normal,
                 log-Logistic, Gompertz, Exponential Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure,
                 log-Logistic Cure, Gompertz Cure."))
  }

  # function returns transition probability
  ### GAMMA, WEIBULL, LOGLOGISTIC ###
  if(dist.v == "Gamma" | dist.v == "Weibull (AFT)" | dist.v == "log-Logistic"){

    # 1) Transform MSM model est
    logPar1 <- log(d.data[1]) # log estimate paramter 1
    logPar2 <- log(d.data[2]) # log estimate paramter 2
    if(length(d.data)>2){     # covariate coefficients
      B       <- d.data[3:length(d.data)]
      op.par  <- c(logPar1, logPar2, B)
    }else {
      op.par  <- c(logPar1, logPar2)
    }
    # 2) Run MSM output in rmvnorm - multivariate normal distribution with mean = op.par and sigma = vc.data (variance cvariance matrix)
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution sample # o simulations

    #3) Transform rmvnorm MSM model est
    t.par1  <- exp(mvnorm.res[,1])  # transform multivariate normal estimates
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### CURE GAMMA, WEIBULL, LOGLOGISTIC ###
  } else if(dist.v == "Gamma Cure" | dist.v == "Weibull (AFT) Cure" | dist.v == "log-Logistic Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1]))       # theta estimate
    logPar2 <- log(d.data[2])                       # log estimate paramter 2
    logPar3 <- log(d.data[3])                       # log estimate paramter 3
    if(length(d.data)>3){
      B       <- d.data[4:length(d.data)]
      op.par  <- c(Par1, logPar2, logPar3, B)
    }else {
      op.par  <- c(Par1, logPar2, logPar3)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution

    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- exp(mvnorm.res[,2])
    t.par3  <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta    <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }

    ### EXPONENTIAL ###
  } else if (dist.v == "Exponential"){

    logPar1 <- log(d.data[1])  # log rate estimate
    if(length(d.data)>1){
      B       <- d.data[2:length(d.data)]
      op.par  <- c(logPar1,  B)
    }else {
      op.par  <- c(logPar1)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- exp(mvnorm.res[,1])
    if(length(d.data)>1){
      Beta    <- mvnorm.res[,2:ncol(mvnorm.res)]
      results <- cbind(t.par1,  Beta)
    }else {
      results <- cbind(t.par1)
    }

    ### CURE EXPONENTIAL ###
  } else if (dist.v == "Exponential Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1])) # theta estimate
    logPar2 <- log(d.data[2])                 # log rate estimate
    if(length(d.data)>2){
      B       <- d.data[3:length(d.data)]
      op.par  <- c(Par1, logPar2, B)
    }else {
      op.par  <- c(Par1, logPar2)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "log-Normal"| dist.v == "Gompertz"){
    Par1    <- d.data[1]      # meanlog(log-Normal). shape(Gompertz)
    logPar2 <- log(d.data[2]) # transform sdlog(log-Normal) or rate(Gompertz) estimate
    if(length(d.data)>2){
      B       <- d.data[3:length(d.data)]
      op.par  <- c(Par1, logPar2, B)
    }else {
      op.par  <- c(Par1, logPar2)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- mvnorm.res[,1]
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### CURE LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "log-Normal Cure"| dist.v == "Gompertz Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1]))      # theta
    Par2    <- d.data[2]      # meanlog(log-Normal). shape(Gompertz)
    logPar3 <- log(d.data[3]) # transform sdlog(log-Normal) or rate(Gompertz) estimate
    if(length(d.data)>3){
      B       <- d.data[4:length(d.data)]
      op.par  <- c(Par1, Par2, logPar3, B)
    }else {
      op.par  <- c(Par1, Par2, logPar3)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- mvnorm.res[,2]
    t.par3  <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta    <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }
  } else {
    print("no distribution found")
    results <- NA
  }
  return(results)
}

#' Randomly draw parameter values of survival distributions from multivariate normal distribution.
#'
#' \code{model.rmvnorm} randomly draws parameter values of survival models from multivariate normal distribution.
#'
#' @param dist.v a character string specifying the name of the survival distribution.
#' @param d.data a vector of mean parameter estimates of the survival distribution.
#' @param vc.data variance-covariance matrix (a matrix) of parameter estimates of the survival distribution.
#' @param n_sim number of random samples to draw.
#' Default = 100.
#' @param seed seed for random number generation.
#' Default = 421.
#' @return
#' A matrix of drawn parameter values.
#' @export
model.rmvnorm <- function(dist.v, d.data, vc.data, n_sim, seed = 421) {

  set.seed(seed)

  if (!dist.v %in% c("Exponential", "Weibull (AFT)", "Gamma", "log-Normal",
                     "log-Logistic", "Gompertz", "Expoenntial Cure", "Weibull (AFT) Cure", "Gamma Cure", "log-Normal Cure",
                     "log-Logistic Cure", "Gompertz Cure")) {
    return(paste0("Incorrect distribution name, select from: Exponential, Weibull (AFT), Gamma, log-Normal,
                log-Logistic, Gompertz, Expoenntial Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure,
                log-Logistic Cure, Gompertz Cure."))
  }

  # function returns transition probability
  ### GAMMA, WEIBULL, LOGLOGISTIC ###
  if(dist.v == "Gamma" | dist.v == "Weibull (AFT)" | dist.v == "log-Logistic"){

    # 1) Transform MSM model est
    logPar1 <- log(d.data[1]) # log estimate paramter 1
    logPar2 <- log(d.data[2]) # log estimate paramter 2
    if(length(d.data)>2){     # covariate coefficients
      B       <- d.data[3:length(d.data)]
      op.par  <- c(logPar1, logPar2, B)
    }else {
      op.par  <- c(logPar1, logPar2)
    }
    # 2) Run MSM output in rmvnorm - multivariate normal distribution with mean = op.par and sigma = vc.data (variance cvariance matrix)
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution sample # o simulations

    #3) Transform rmvnorm MSM model est
    t.par1  <- exp(mvnorm.res[,1])  # transform multivariate normal estimates
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### CURE GAMMA, WEIBULL, LOGLOGISTIC ###
  } else if(dist.v == "Gamma Cure" | dist.v == "Weibull (AFT) Cure" | dist.v == "log-Logistic Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1]))       # theta estimate
    logPar2 <- log(d.data[2])                       # log estimate paramter 2
    logPar3 <- log(d.data[3])                       # log estimate paramter 3
    if(length(d.data)>3){
      B       <- d.data[4:length(d.data)]
      op.par  <- c(Par1, logPar2, logPar3, B)
    }else {
      op.par  <- c(Par1, logPar2, logPar3)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution

    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- exp(mvnorm.res[,2])
    t.par3  <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta    <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }

    ### EXPONENTIAL ###
  } else if (dist.v == "Exponential"){

    logPar1 <- log(d.data[1])  # log rate estimate
    if(length(d.data)>1){
      B       <- d.data[2:length(d.data)]
      op.par  <- c(logPar1,  B)
    }else {
      op.par  <- c(logPar1)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- exp(mvnorm.res[,1])
    if(length(d.data)>1){
      Beta    <- mvnorm.res[,2:ncol(mvnorm.res)]
      results <- cbind(t.par1,  Beta)
    }else {
      results <- cbind(t.par1)
    }

    ### CURE EXPONENTIAL ###
  } else if (dist.v == "Exponential Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1])) # theta estimate
    logPar2 <- log(d.data[2])                 # log rate estimate
    if(length(d.data)>2){
      B       <- d.data[3:length(d.data)]
      op.par  <- c(Par1, logPar2, B)
    }else {
      op.par  <- c(Par1, logPar2)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "log-Normal"| dist.v == "Gompertz"){
    Par1    <- d.data[1]      # meanlog(log-Normal). shape(Gompertz)
    logPar2 <- log(d.data[2]) # transform sdlog(log-Normal) or rate(Gompertz) estimate
    if(length(d.data)>2){
      B       <- d.data[3:length(d.data)]
      op.par  <- c(Par1, logPar2, B)
    }else {
      op.par  <- c(Par1, logPar2)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- mvnorm.res[,1]
    t.par2  <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta    <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }

    ### CURE LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "log-Normal Cure"| dist.v == "Gompertz Cure"){
    Par1    <- log(d.data[1]/(1 - d.data[1]))      # theta
    Par2    <- d.data[2]      # meanlog(log-Normal). shape(Gompertz)
    logPar3 <- log(d.data[3]) # transform sdlog(log-Normal) or rate(Gompertz) estimate
    if(length(d.data)>3){
      B       <- d.data[4:length(d.data)]
      op.par  <- c(Par1, Par2, logPar3, B)
    }else {
      op.par  <- c(Par1, Par2, logPar3)
    }

    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)

    t.par1  <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2  <- mvnorm.res[,2]
    t.par3  <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta    <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }
  } else {
    print("no distribution found")
    results <- NA
  }
  return(results)
}

#' Calculate survival probabilities given a survival distribution and its parameter values.
#'
#' \code{model.dist} calculates survival probabilities given a survival distribution and its parameter values.
#'
#' @param dist.v a character string specifying the name of the survival distribution.
#' @param d.data a vector of parameter values of the survival distribution.
#' @param t a vector of time points to calculate the survival probabilities at.
#' @return
#' A vector of survival probabilities.
#' @export
model.dist <- function(dist.v, d.data, t){

  # save parameters as pars and coefficients as beta

  v_par1 <- c("Exponential")
  v_par2 <- c("Gamma", "Weibull", "log-Normal", "log-Logistic", "Gompertz", "Exponential Cure")
  v_par3 <- c("Gamma Cure", "Weibull Cure", "log-Normal Cure", "log-Logistic Cure", "Gompertz Cure")

  if (dist.v %in% v_par1){
    pars <- d.data[1]
    beta <- if(length(d.data) > 1){d.data[2:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]
  } else if(dist.v %in% v_par2){
    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]
  } else if(dist.v %in% v_par3){
    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]
  }

  # Get cumulative survival probability for each individual
  ############################ Gamma ############################
  if(dist.v == "Gamma"){
    for (i in 1:length(pars)) {
      if(i == 2){
        # print("rate") # beta affects the rate

        if(length(beta) == 0){beta.raw = 0 # if no covariates
        }else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      # else{print("Gamma")}
    }
    p.res   <- pgamma(t,        shape = pars[1], rate = pred, lower.tail = F)
    # p.res.1 <- pgamma(t - step, shape = pars[1], rate = pred, lower.tail = F)

    ############################ Gamma Cure ############################
  } else if (dist.v == "Gamma Cure"){
    for (i in 1:length(pars)){
      if (i == 3){
        print("rate") # beta affects the rate

        if(length(beta) == 0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)           # fit$dlist$inv.transforms
      }
      # else{print("Gamma Cure")}
    }
    p.res   <- pmixsurv(pgamma, t,        theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)
    # p.res.1 <- pmixsurv(pgamma, t - step, theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)

    ############################ Exponential ############################
  } else if (dist.v == "Exponential") {
    for (i in 1:length(pars)) {
      # print("rate") # beta affects the rate

      if(length(beta) ==0){beta.raw = 0}
      else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

      pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
      pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
    }

    p.res   <- pexp(t,        rate = pred, lower.tail = F)
    # p.res.1 <- pexp(t - step, rate = pred, lower.tail = F)

    ############################ Exponential Cure ############################
  } else if (dist.v == "Exponential Cure") {
    for (i in 1:length(pars)) {
      if (i == 2){
        # print("rate") # Beta affects the rate

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      # else{print("Exponential Cure")}
    }
    p.res   <- pmixsurv(pexp, t,        theta = pars[1], rate = pred, lower.tail = F)
    # p.res.1 <- pmixsurv(pexp, t - step, theta = pars[1], rate = pred, lower.tail = F)

    ############################ Weibull ############################
  } else if (dist.v == "Weibull") {
    for (i in 1:length(pars)) {
      if(i == 2) {
        # print("scale") # Beta affects the scale

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      # else{print("Weibull (AFT)")}
    }
    p.res   <- pWeibull (AFT)(t,        shape = pars[1], scale = pred, lower.tail = F)
    # p.res.1 <- pWeibull (AFT)(t - step, shape = pars[1], scale = pred, lower.tail = F)

    ############################ Weibull Cure ############################
  } else if (dist.v == "Weibull (AFT) Cure") {
    for (i in 1:length(pars)) {
      if(i == 3) {
        # print("scale") # Beta affects the scale

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      # else{print("Weibull (AFT) Cure")}
    }
    p.res   <- pmixsurv(pWeibull (AFT), t,        theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)
    # p.res.1 <- pmixsurv(pWeibull (AFT), t - step, theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)

    ############################ log-Normal ############################
  } else if (dist.v == "log-Normal") {
    for (i in 1:length(pars)) {
      if(i == 1) {
        # print("meanlog")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- pars[i] + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- pred.raw              # fit$dlist$inv.transforms
      }
      # else{print("log-Normal")}
    }
    p.res   <- plnorm(t,        meanlog = pred, sdlog = pars[2], lower.tail = F)
    # p.res.1 <- plnorm(t - step, meanlog = pred, sdlog = pars[2], lower.tail = F)

    ############################ log-Normal Cure ############################
  } else if (dist.v == "log-Normal Cure") {
    for (i in 1:length(pars)) {
      if(i == 2) {
        # print("meanlog")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- pars[i] + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- pred.raw              # fit$dlist$inv.transforms
      }
      # else{print("log-Normal Cure")}
    }

    p.res   <- pmixsurv(plnorm, t       , theta = pars[1], meanlog = pred, sdlog = pars[3], lower.tail = F)
    # p.res.1 <- pmixsurv(plnorm, t - step, theta = pars[1], meanlog = pred, sdlog = pars[3], lower.tail = F)
    ############################ Gompertz ############################
  } else if (dist.v == "Gompertz") {
    for (i in 1:length(pars)) {
      if(i == 2) {
        # print("rate")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      # else{print("Gompertz")}
    }
    p.res   <- pgompertz(t,        shape = pars[1], rate = pred, lower.tail = F)
    # p.res.1 <- pgompertz(t - step, shape = pars[1], rate = pred, lower.tail = F)

    ############################ Gompertz Cure ############################
  } else if (dist.v == "Gompertz Cure") {
    for (i in 1:length(pars)) {
      if(i == 3) {
        # print("rate")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      # else{print("Gompertz Cure")}
    }
    p.res   <- pmixsurv(pgompertz, t,        theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)
    # p.res.1 <- pmixsurv(pgompertz, t - step, theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)

    ############################ log-Logistic ############################
  } else if (dist.v == "log-Logistic") {
    for (i in 1:length(pars)) {
      if(i == 2) {
        # print("scale")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      # else{print("log-Logistic")}
    }
    p.res   <- pllogis(t,        shape = pars[1], scale = pred, lower.tail = F)
    # p.res.1 <- pllogis(t - step, shape = pars[1], scale = pred, lower.tail = F)

    ############################ log-Logistic Cure ############################
  } else if (dist.v == "log-Logistic Cure") {
    for (i in 1:length(pars)) {
      if(i == 3) {
        # print("scale")

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      # else{print("log-Logistic Cure")}
    }

    p.res   <- pmixsurv(pllogis, t,        theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)
    # p.res.1 <- pmixsurv(pllogis, t - step, theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)

    ############################ No Distribution ############################
  } else {
    print("no distribution found")
    p.res <- NA
    # p.res.1 <- NA
  }

  return(p.res)   # return survival probabilities
}
