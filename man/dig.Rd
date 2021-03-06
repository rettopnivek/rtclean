% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wald_functions.R
\name{dig}
\alias{dig}
\title{Density for the Inverse Gaussian Distribution}
\usage{
dig(t, kappa, xi, ln = F, summation = F)
}
\arguments{
\item{t}{vector of response times (t > 0).}

\item{kappa}{a scalar threshold value (kappa > 0).}

\item{xi}{a scalar rate of evidence accumulation (xi > 0).}

\item{ln}{logical; If TRUE, returns the log of the density.}

\item{summation}{logical; if TRUE, returns the sum of the logs of the
densities.}
}
\value{
If \code{ln} is FALSE, gives the likelihood, else gives the log-likelihood.
  If both \code{ln} and \code{summation} are TRUE, gives the sum of the log-likelihoods.
}
\description{
The density function for an inverse gaussian distribution,
parameterized in terms of a one-boundary wiener process.
}
\details{
The inverse gaussian, when parameterized in terms of brownian motion,
  has density \deqn{ f(t) = \kappa/\sqrt( 2*\pi*t^3 ) exp( (-.5/t)*(\kappa - \xi*t)^2 ) }
  where \eqn{t} is a response time, and \eqn{\kappa} is a threshold toward which evidence
  accumulates with average rate \eqn{\xi} and a fixed variance of 1.
}

