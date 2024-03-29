% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_functions.R
\name{partsurv}
\alias{partsurv}
\title{Fit partitioned survival model to survival data}
\usage{
partsurv(
  pfs_survHE = NULL,
  os_survHE = NULL,
  l_d.data = NULL,
  l_vc.data = NULL,
  par = FALSE,
  chol = FALSE,
  choose_PFS = NULL,
  choose_OS = NULL,
  time = times,
  v_names_states,
  PA = FALSE,
  n_sim = 100,
  seed = 421,
  warn = TRUE,
  dat.x = 0
)
}
\arguments{
\item{pfs_survHE}{survHE obj fitting PFS.}

\item{os_survHE}{survHE obj fitting OS.}

\item{l_d.data}{list of mean parameter estimates (list containing 2 numerical estimates, 1st being for PFS and 2nd being for OS).}

\item{l_vc.data}{list of variance-covariance matrices (or their Cholesky decomposition) of parameter estimates (list containing 2 matrices, 1st being for PFS and 2nd being for OS).}

\item{par}{set to TRUE if parameter mean estimates and their variance-covariance matrices are used instead of survHE objects.
Default = FALSE}

\item{chol}{set to TRUE if l_vc.data contains Cholesky decomposition of the variance-covariance matrices instead of the actual variance-covariance matrices.
Default = FALSE}

\item{choose_PFS}{chosen PFS distribution. Choose from: Exponential, Weibull (AFT), Gamma, log-Normal, log-Logistic, Gompertz, Exponential Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure, log-Logistic Cure, Gompertz Cure.}

\item{choose_OS}{chosen OS distribution. Choose from: Exponential, Weibull (AFT), Gamma, log-Normal, log-Logistic, Gompertz, Exponential Cure, Weibull (AFT) Cure, Gamma Cure, log-Normal Cure, log-Logistic Cure, Gompertz Cure.}

\item{time}{numeric vector of time to estimate probabilities.}

\item{v_names_states}{vector of state names.}

\item{PA}{run probabilistic analysis.
Default = FALSE.}

\item{n_sim}{number of PA simulations.
Default = 100.}

\item{seed}{seed for random number generation.
Default = 421.}

\item{warn}{prints a warning message whenever PFS > OS}
}
\value{
a list containing Markov trace, expected survival, survival probabilities, transition probabilities.
}
\description{
\code{partsurv} fits partitioned survival model to survival data.
}
