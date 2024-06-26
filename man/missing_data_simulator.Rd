% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simMFPCA.R
\name{missing_data_simulator}
\alias{missing_data_simulator}
\title{Function to simulate missingness in simulated data}
\usage{
missing_data_simulator(n_obs, dat, t = seq(0, 1, length.out = 51), seed = 51)
}
\arguments{
\item{n_obs}{number of observations to keep}

\item{dat}{original data}

\item{t}{time vector}

\item{seed}{random seed}
}
\value{
a data with artificial missing value
}
\description{
Function to simulate missingness in simulated data
}
