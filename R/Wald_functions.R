#
#
#

# Index
# Lookup - 01:  dig
# Lookup - 02:  dsig
# Lookup - 03:  dmsig
# Lookup - 04:  rtclean

# Useful functions for package creation
# library(devtools)
# library(roxygen2)

# Lookup - 01
#' Density for the Inverse Gaussian Distribution
#'
#' The density function for an inverse gaussian distribution,
#' parameterized in terms of a one-boundary wiener process.
#'
#' @param t vector of response times (t > 0).
#' @param kappa a scalar threshold value (kappa > 0).
#' @param xi a scalar rate of evidence accumulation (xi > 0).
#' @param ln logical; If TRUE, returns the log of the density.
#' @param summation logical; if TRUE, returns the sum of the logs of the
#'   densities.
#' @details The inverse gaussian, when parameterized in terms of brownian motion,
#'   has density \deqn{ f(t) = \kappa/\sqrt( 2*\pi*t^3 ) exp( (-.5/t)*(\kappa - \xi*t)^2 ) }
#'   where \eqn{t} is a response time, and \eqn{\kappa} is a threshold toward which evidence
#'   accumulates with average rate \eqn{\xi} and a fixed variance of 1.
#' @return If \code{ln} is FALSE, gives the likelihood, else gives the log-likelihood.
#'   If both \code{ln} and \code{summation} are TRUE, gives the sum of the log-likelihoods.
#' @export

dig <- function( t, kappa, xi, ln = F, summation = F ) {

  if ( length(kappa) > 1 | length(xi) > 1 ) stop('parameters must be scalars')

  # Initialize output
  out = rep( -Inf, length(t) )

  # Check for inadmissable values
  if ( (kappa <= 0) | (xi <= 0 ) ) {
    if (!ln) {
      out = exp(out)
    } else {
      if (summation) {
        out = sum(out)
      }
    }
  } else {

    # Response times can't be negative
    t[ t < 0 ] = 0

    # Calculate the log of the density for the inverse gaussian
    out = log( kappa/sqrt( 2*pi*t^3) ) +
      ( (-.5/t) * (kappa - t*xi)^2 )

    # Check for inadmissable inputs
    out[ is.na(out) ] = -Inf
    out[ t == 0 ] = -Inf
  }

  if (!ln) {
    out = exp( out ); # If indicated, return the density
  } else {
    # If indicated, return the sum of the log-likelihoods
    if (summation) {
      out = sum( out )
      if ( is.na( out ) ) out = -Inf
    }
  }

  return( out )
}

# Lookup - 02
#' Density for the Shifted Inverse Gaussian
#'
#' The density function for a shifted inverse gaussian distribution,
#' parameterized in terms of a one-boundary wiener process.
#'
#' @param t vector of response times (t > 0).
#' @param kappa a scalar threshold value (kappa > 0).
#' @param xi a scalar rate of evidence accumulation (xi > 0).
#' @param tau a scalar residual latency by which to shift the response times.
#'   ( 0 <= tau < min(t) ).
#' @param ln logical; If TRUE, returns the log of the density.
#' @param summation if TRUE, returns the sum of the logs of the
#'   densities.
#' @details The shifted inverse gaussian, when parameterized in terms of brownian motion,
#'   has density \deqn{ f(t - \tau) = \kappa/\sqrt( 2*\pi*(t-\tau)^3 ) exp( (-.5/(t-\tau))*(\kappa - \xi*(t-\tau))^2 ) }
#'   where \eqn{t} is a response time, \eqn{\kappa} is a threshold toward which evidence
#'   accumulates with average rate \eqn{\xi} and a fixed variance of 1, and \eqn{\tau} is a
#'   residual latency subtracted from the response time.
#' @return If \code{ln} is FALSE, gives the likelihood, else gives the log-likelihood.
#'   If both \code{ln} and \code{summation} are TRUE, gives the sum of the log-likelihoods.
#' @export

dsig <- function( t, kappa, xi, tau, ln = F, summation = F ) {

  if ( length(kappa) > 1 | length(xi) > 1 ) stop('parameters must be scalars')

  # Initialize output
  out = rep( -Inf, length(t) )

  # Check for inadmissable values
  if ( (kappa <= 0) | (xi <= 0 ) | (tau < 0) ) {
    if (!ln) {
      out = exp(out)
    } else {
      if (summation) {
        out = sum(out)
      }
    }
  } else {

    # Shift response time
    st = t - tau
    st[ st < 0 ] = 0

    # Calculate the log of the density for the inverse gaussian
    out = log( kappa/sqrt( 2*pi*st^3) ) +
      ( (-.5/st) * (kappa - st*xi)^2 )

    # Check for inadmissable inputs
    out[ is.na(out) ] = -Inf
    out[ t <= tau ] = -Inf
  }

  if (!ln) {
    out = exp( out ); # If indicated, return the density
  } else {
    # If indicated, return the sum of the log-likelihoods
    if (summation) {
      out = sum( out )
      if ( is.na( out ) ) out = -Inf
    }
  }

  return( out )
}

# Lookup - 03
#' Density for the Mixture of an Inverse Gaussian and a Uniform
#' Distribution
#'
#' The density function for the mixture of a shifted inverse
#' gaussian and a uniform distribution. The inverse gaussian
#' distribution is parameterized in terms of a one-boundary
#' wiener process.
#'
#' @param t vector of response times (t > 0).
#' @param kappa a scalar threshold value (kappa > 0).
#' @param xi a scalar rate of evidence accumulation (xi > 0).
#' @param tau a scalar residual latency by which to shift the response times.
#'   ( 0 <= tau < min(t) ).
#' @param lambda a scalar mixture probability for the inverse gaussian.
#' @param alpha a scalar for the lower boundary to use in the uniform density.
#' @param beta a scalar for the uppber boundary to use in the uniform density.
#' @param ln logical; If TRUE, returns the log of the density.
#' @param summation if TRUE, returns the sum of the logs of the
#'   densities.
#' @param sep logical; if TRUE, returns a matrix with separate columns for the weighted
#'   likelihoods of the inverse gaussian and the uniform distribution.
#' @return If \code{sep} is TRUE, gives a matrix whose columns are the separate weighted
#'   densities for the inverse gaussian and uniform distribution, respectively. If
#'   \code{ln} is FALSE, gives the likelihood, else gives the log-likelihood.
#'   If both \code{ln} and \code{summation} are TRUE, gives the sum of the log-likelihoods.
#' @export

dmsig <- function( t, kappa, xi, tau, lambda,
                   alpha = min(t), beta = max(t), ln = F,
                   summation = F, sep = F ) {

  # Calculate density for inverse gaussian distribution
  D = dsig( t, kappa, xi, tau, ln = T )
  # Calculate density for uniform distribution
  U = dunif( t, alpha, beta, log = 1 )

  # Initialize output
  out = rep(0,length(D))

  if (sep) {
    # If there are no inadmissable values
    if ( (kappa > 0) & (xi > 0 ) & (tau >= 0) &
         (lambda >= 0 & lambda <= 1 ) ) {
      out = cbind( log( lambda ) + D, log( 1 - lambda ) + U )
      out = exp( out )
      colnames(out) = c('IG','U')
    }
  } else {
    # If there are no inadmissable values
    if ( (kappa > 0) & (xi > 0 ) & (tau >= 0) &
         (lambda >= 0 & lambda <= 1 ) ) {
      out = exp( log( lambda ) + D ) + exp( log( 1 - lambda ) + U )
    }
    if (ln) {
      out = log( out );
      if (summation) {
        out = sum( out )
        if ( is.na( out ) ) out = -Inf;
      }
    }
  }

  return( out )
}

#' Function for Trimming Outlier Response Times
#'
#' Forthcoming
#'
#' @param rt a vector of response times.
#' @param exclude if \code{exclude} equals 1, any response time that has a relative
#'   probability less that 0.5 under the inverse gaussian distribution is trimmed. If
#'   \code{exclude} equals 2, response times are probabilistically trimmed using the
#'   relative probabilities under the inverse gaussian and uniform distribution.
#' @param lowerBoundary a vector defining the lower boundary for the range over which
#'   the starting values are generated.
#' @param upperBoundary a vector defining the upper boundary for the range over which
#'   the starting values are generated.
#' @param method determines the estimation method used by R's \code{optim} function.
#' @param maxit the maximum number of iterations to run in the optimization routine.
#' @param nRep the number of times to repeat the optimization routine with new starting
#'   values (higher values lead to a greater chance of avoiding local maxima at the cost
#'   of longer computational times).
#' @param alpha the width of the confidence intervals for parameter estimates.
#' @param cutoff a cutoff for excluding observations based on their relative
#'   probability under the inverse gaussian (default is .05).
#' @return \code{rtclean} returns an object of class "cleanrt".
#'
#' The function \code{summary} prints a comparison between the descriptive statistics
#' for the original and trimmed data sets, \code{anova} reports the Akaike Information
#' Criterion (AIC) values (corrected for small samples) and their relative weights for
#' the set of 3 models fit to the data.
#'
#' An object of class "cleanrt" is a list containing the following components:
#' \describe{
#'   \item{\code{rtNew}}{ a vector with the newly trimmed response times.}
#'   \item{\code{rtOld}}{a vector of the untrimmed response times.}
#'   \item{\code{M1}}{the \code{optim} results for the inverse gaussian model.}
#'   \item{\code{M2}}{the \code{optim} results for the shifted inverse gaussian model.}
#'   \item{\code{M3}}{the \code{optim} results for the mixture model.}
#'   \item{\code{M3}}{a list with the AICc values and relative weights.}
#'   \item{\code{ProbIG}}{a vector with the relative probabilities of each observation
#'     under the shifted inverse gaussian distribution.}
#'   \item{\code{exclude_1}}{a logical vector for the response times selected for
#'     exclusion based on having a relative probability of less than 0.5 for the
#'     shifted inverse gaussian.}
#'   \item{\code{exclude_2}}{a logical vector for the response times selected for
#'     exclusion based on probabilistic sampling using the relative probabilities
#'     for the inverse gaussian and uniform distributions.}
#' }

#' @export
rtclean = function( rt, exclude = 1, lowerBoundary=c(.5,.5,0,.8),
                    upperBoundary=c(5,5,.3,.95),
                    method = 'BFGS', maxit = 5000,
                    nRep = 5, alpha = .95, cutoff = -2.9957 ) {

  # Determine the smallest response time
  minVal = min( rt );
  # Set the range of starting values
  startRange = list( lowerBoundary, upperBoundary )

  # Fit the inverse gaussian
  type = 1
  model_1 = try( mleWrapper( rt, type,
                        mle_f, st_f, startRange, gr = gr_f,
                        method = method, maxit = maxit,
                        nRep = nRep, alpha = alpha ), silent = T )

  # Fit the shifted inverse gaussian
  type = 2
  model_2 = try( mleWrapper( rt, type,
                        mle_f, st_f, startRange, gr = gr_f,
                        method = method, maxit = maxit,
                        nRep = nRep, alpha = alpha ), silent = T )

  # Fit the mixture model
  # Note that we double the number of repetitions since the mixture model
  # can be unstable
  type = 3
  model_3 = try( mleWrapper( rt, type,
                        mle_f, st_f, startRange, gr = gr_f,
                        method = method, maxit = maxit,
                        nRep = nRep*2, alpha = alpha ), silent = T )

  # Check whether models were successfully estimated
  errorMessage = NULL
  if ( !is.list(model_1) ) {
    print( 'Warning: Inverse gaussian failed to fit' )
    errorMessage = list( errorMessage, IG = model_1 )
    model_1 = list( NULL )
  }
  if ( !is.list(model_2) ) {
    print( 'Warning: Shifted inverse gaussian failed to fit' )
    errorMessage = list( errorMessage, sIG = model_2 )
    model_2 = list( NULL )
  }
  if ( !is.list(model_2) ) {
    print( 'Warning: Mixture model failed to fit' )
    errorMessage = list( errorMessage, Mix = model_3 )
    model_3 = list( NULL )
  }

  # Compare the models using the AICc
  allAIC = c(model_1$AICc,model_2$AICc,model_3$AICc)
  names(allAIC) = c('M1','M2','M3')
  deltaAIC = allAIC - min(allAIC)
  relativeL = exp( -.5*deltaAIC )
  AICweights = relativeL/sum( relativeL )

  # Calculate the likelihoods for the separate distributions in
  # the mixture
  if ( length( model_3 ) > 0 ) {

    M3coef = model_3$par
    Obslik = dmsig( rt, M3coef[1], M3coef[2], M3coef[3], M3coef[4],
                    sep = T )
    ProbIG = Obslik[,1]/rowSums(Obslik)
    relativeEvidence = log(Obslik[,1]) - log(Obslik[,2])

    exclude_1 = relativeEvidence <= cutoff
    exclude_2 = apply( Obslik, 1,
                       function(x) sample( c(F,T), 1,
                                           prob = x[1:2]/sum(x[1:2]) ) )
    if (exclude == 1) {
      rtNew = rt[ !exclude_1 ]
    }
    if (exclude == 2) {
      rtNew = rt[ !exclude_2 ]
    }
  }

  output = list(
      rtNew = rtNew,
      rtOld = rt,
      M1 = model_1,
      M2 = model_2,
      M3 = model_3,
      AICc = list( AICc = allAIC, weights = AICweights ),
      ProbIG = ProbIG,
      exclude_1 = exclude_1,
      exclude_2 = exclude_2
      )

  # Create new class for output
  class(output) = append( class(output), "cleanrt" )

  invisible( output )
}

# Print method for new class
print.cleanrt = function( object ) {
  rng = range( object$rtNew )
  rmv = length(object$rtNew)/length(object$rtOld)

  out = paste( "Range: ", round(rng[1],3), " to ",
               round(rng[2],3),"; ", 100*round(1-rmv,2),
               "% removed.", "\n", sep = "" )

  cat( out )
}

# Summary method for new class
summary.cleanrt = function( object ) {

  rtOld = object$rtOld
  rtNew = object$rtNew
  rmv = length(rtNew)/length(rtOld)

  cat('Descriptives:','\n')
  cat('------------------------------------','\n')
  qntOld = quantile( rtOld )
  qntNew = quantile( rtNew )
  tmp = rbind( c(qntOld,mean(rtOld),sd(rtOld)),
               c(qntNew,mean(rtNew),sd(rtNew)) )
  rownames(tmp) = c('Original:','Trimmed:')
  colnames(tmp) = c('Min.','1st Q.','Median','3rd Q.','Max.','Mean','SD')
  printCoefmat( round(tmp,3) )
  cat( paste(100*round(1-rmv,2),'% removed',sep=''), '\n' )
  cat('------------------------------------','\n')

  cat('Parameter estimates:','\n')
  cat('------------------------------------','\n')
  tmp = rbind( c(object$M1$par,NA,NA),
               c(object$M2$par,NA),
               c(object$M3$par) )
  colnames(tmp) = c('kappa','xi','tau','lambda')
  rownames(tmp) = c('IG','sIG','msIG+U')
  printCoefmat( tmp, digits=2, na.print = ' ' )
  cat('------------------------------------','\n')
}

# Coef method for new class
# Extracts coefficients, standard errors, and confidence intervals
coef.cleanrt = function( object ) {
  list( par = list(
    M1 = object$M1$par,
    M2 = object$M2$par,
    M3 = object$M3$par ),
    SE = list(
      M1 = object$M1$SE,
      M2 = object$M2$SE,
      M3 = object$M3$SE ),
    CI = list(
      M1 = object$M1$CI,
      M2 = object$M2$CI,
      M3 = object$M3$CI )
  )
}

# ANOVA method for new class
# Displays AIC and relative weights
anova.cleanrt = function( object ) {

  aic = object$AICc

  cat('Akaike Information Criterion','\n')
  cat('(corrected for small samples)','\n')
  cat('------------------------------------','\n')

  values = round(aic[[1]])
  weights = round(aic[[2]],2)

  tmp = format( rbind( c(' ','IG','sIG','msIG+U' ),
    c('Values:',as.character(values)),
    c('Weights',as.character(weights)) ), justify = 'centre')

  cat( tmp[1,], '\n' )
  cat( tmp[2,], '\n' )
  cat( tmp[3,], '\n' )

  sel = which( aic[[2]] ==  max(aic[[2]]) )
  if (sel==1) model = 'Inverse gaussian'
  if (sel==2) model = 'Shifted inverse gaussian'
  if (sel==3) model = 'Mixture model'

  cat('------------------------------------','\n')
  cat( paste(model, ' is preferred.',sep=''),'\n' )
}
