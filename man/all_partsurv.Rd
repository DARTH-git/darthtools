% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_functions.R
\name{all_partsurv}
\alias{all_partsurv}
\title{Fit partitioned survival model on all combinations of chosen PFS and OS parametric survival functions.}
\usage{
all_partsurv(
  pfs_survHE,
  os_survHE,
  choose_PFS,
  choose_OS,
  time = times,
  PA = FALSE,
  n_sim = 100,
  seed = 421
)
}
\arguments{
\item{pfs_survHE}{survHE obj fitting PFS.}

\item{os_survHE}{survHE obj fitting OS.}

\item{choose_PFS}{preferred PFS distribution.}

\item{choose_OS}{preferred OS distribution.}

\item{time}{numeric vector of time to estimate probabilities.}

\item{PA}{run probabilistic analysis.
Default = FALSE.}

\item{n_sim}{number of PA simulations.
Default = 100.}

\item{seed}{seed for random number generation.
Default = 421.}
}
\value{
a list containing Markov trace, expected survival, survival probabilities, transition probabilities.
}
\description{
\code{all_partsurv} fits partitioned survival model on all combinations of chosen PFS and OS parametric survival functions.
}
