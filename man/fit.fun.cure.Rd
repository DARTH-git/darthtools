% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_functions.R
\name{fit.fun.cure}
\alias{fit.fun.cure}
\title{Fit multiple mixture cure models survival data}
\usage{
fit.fun.cure(time, status, data = data, extrapolate = FALSE, times)
}
\arguments{
\item{time}{numeric vector of time to estimate probabilities.}

\item{status}{numeric vector of event status.}

\item{data}{dataframe containing the time and status variables.}

\item{extrapolate}{extrapolate beyond model time horizon.}

\item{times}{time horizon the extrapolation is done over.}
}
\value{
a list containing all survival model objects.
}
\description{
\code{fit.fun.cure} fits multiple mixure cure models to survival data using flexsurv.
}
