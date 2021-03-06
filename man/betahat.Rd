% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{betahat}
\alias{betahat}
\title{This function Local Linear Kernel Estimator vector of beta(t0) (p).}
\usage{
betahat(t0, X, y, t, h, mv)
}
\arguments{
\item{t0}{time point for evaluation}

\item{X}{is a covariates matrix (N*p),}

\item{y}{is the response vector (N),}

\item{t}{is the time vector (N)
For one variable, the longitudinals are vetorized together.}

\item{h}{is the seleted bandwidth.}

\item{mv}{is a vector (N) telling us how many repeated measurements for each individual}
}
\value{
Returns betahat(t0)
}
\description{
This function Local Linear Kernel Estimator vector of beta(t0) (p).
}
