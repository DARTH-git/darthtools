% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_transition_probabilities.R
\name{check_sum_of_transition_array}
\alias{check_sum_of_transition_array}
\title{Check if the sum of transition probabilities equal to one.}
\usage{
check_sum_of_transition_array(
  a_P,
  n_rows = NULL,
  n_states = NULL,
  n_cycles,
  err_stop = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{a_P}{A transition probability array/ matrix.}

\item{n_rows}{Number of rows (individuals), appropriate for microsimulation models.}

\item{n_states}{Number of health states in a Markov trace, appropriate for Markov models.}

\item{n_cycles}{Number of cycles.}

\item{err_stop}{Logical variable to stop model run if set up as TRUE.
Default = TRUE.}

\item{verbose}{Logical variable to indicate print out of messages.
Default = TRUE}
}
\value{
The transition probability array and the cohort trace matrix.
}
\description{
\code{check_sum_of_transition_array} checks if each of the rows of the
transition matrices sum to one.
}
