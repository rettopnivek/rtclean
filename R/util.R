#-------------------#
# Utility functions #
#-------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01:  logistic
# Lookup - 02:  logit
# Lookup - 03:  igDerivative
# Lookup - 04:  sigDerivative
# Lookup - 05:  msigDerivative

# Lookup - 01
logistic = function(x) {
  # Purpose:
  # Forthcoming
  # Arguments:
  # Forthcoming

  1/(1+exp(-x))
}

# Lookup - 02
logit = function(p) {
  # Purpose:
  # Forthcoming
  # Arguments:
  # Forthcoming

  log( p/(1-p) )
}

# Lookup - 03
igDerivative = function( t, kappa, xi, summation = T, tran = F ) {
  # Purpose:
  # A function that calculates the first partial derivaties
  # for the log-likelihood of the shifted inverse gaussian
  # Arguments:
  # t         - A vector of response times
  # kappa     - A threshold value (kappa > 0)
  # xi        - A drift rate (xi > 0)
  # summation - If true, sums the derivatives across observations
  #             for each parameter
  # tran      - If true, transforms the inputted parameters
  # Returns:
  # A vector of 2 values, each parameter's partial derivative for
  # the sum of the log-likelihoods. If 'summation' is set to false
  # returns a matrix with the derivative for each individual
  # log-likelihood.

  # Select longest vector
  sz = c( length(t), length(kappa),
          length(xi) )
  N = max( sz )

  if (tran) {
    orig_kappa = kappa
    orig_xi = xi
    kappa = exp( kappa )
    xi = exp( xi )
  }

  # Check for inadmissable values
  chk = (t<0)
  if (sz[2]==1) chk = chk + rep( kappa <= 0, N ) else
    chk = chk + (kappa <= 0)
  if (sz[3]==1) chk = chk + rep( xi <= 0, N ) else
    chk = chk + (xi <= 0)

  # Initialize output
  out = matrix( NA, N, 2 )

  # Set values of t < 0 to 0
  t[ t < 0 ] = 0
  chk = chk + (t==0)

  ### Derivative for kappa ###

  # f = log( ( kappa/sqrt( 2*pi*t^3) )*exp( (-.5/t) * (kappa - t*xi)^2 ) )
  # f = log( kappa ) - log( sqrt( 2*pi*t^3) ) + ( (-.5/t)*( kappa - (t*xi) )^2 )
  # Let K1 = log( sqrt( 2*pi*t^3) ); K2 = -.5/t; K3 = t*xi
  # h = kappa - K3
  # f = log( kappa ) - K1 + K2*( h )^2
  # h' = 1
  # f' = (1/kappa) + 2*K2*( h )*h'
  # f' = (1/kappa) + 2*K2*( kappa - K3 )
  # To check in Wolfram Alpha: log( x ) - a + b*( x - c )^2

  # Define constants
  K2 = -.5/t
  K3 =  t*xi

  # Calculate derivative
  der_kappa = (1/kappa) + 2*K2*( kappa - K3 )
  der_kappa[ chk > 0 ] = NA
  if (tran) der_kappa = der_kappa*kappa

  ### Derivative for xi ###

  # log( ( kappa/sqrt( 2*pi*t^3) )*exp( (-.5/t) * (kappa - t*xi)^2 ) )
  # [ log( kappa/sqrt( 2*pi*t^3) ) ] + ( [-.5/t]*( [ kappa ] - [t]*xi )^2 )
  # h = kappa - t*xi
  # f = K1 + K2*( h )^2
  # h' = -t
  # f' = 2*K2*( h )*h'
  # f' = 2*K2*( kappa - t*xi )*(-t)
  # f' = -2*K2*t*( kappa - t*xi )
  # To check in Wolfram Alpha: a + b*( c - d*x )^2

  # Define constants
  K2 = -.5/t

  # Calculate derivative
  der_xi = -2*K2*t*(kappa-t*xi)
  der_xi[ chk > 0 ] = NA
  if (tran) der_xi = der_xi*xi

  # Store results into matrix
  out[,1] = der_kappa
  out[,2] = der_xi
  colnames( out ) = c( 'dkappa','dxi' )

  if (summation) out = colSums( out )

  out
}

# Lookup - 04
sigDerivative = function( t, kappa, xi, tau, summation = T, tran = F ) {
  # Purpose:
  # A function that calculates the first partial derivaties
  # for the log-likelihood of the shifted inverse gaussian
  # Arguments:
  # t         - A vector of response times
  # kappa     - A threshold value (kappa > 0)
  # xi        - A drift rate (xi > 0)
  # tau       - A residual latency (0 <= tau <= min(t))
  # summation - If true, sums the derivatives across observations
  #             for each parameter
  # tran      - If true, transforms the inputted parameters
  # Returns:
  # A vector of 3 values, each parameter's partial derivative for
  # the sum of the log-likelihoods. If 'summation' is set to false
  # returns a matrix with the derivative for each individual
  # log-likelihood.

  # Select longest vector
  sz = c( length(t), length(kappa),
          length(xi), length(tau) )
  N = max( sz )

  if (tran) {
    minVal = min( t )
    kappa = exp( kappa )
    xi = exp( xi )
    orig_tau = tau
    tau = logistic(orig_tau)*minVal
  }

  # Check for inadmissable values
  chk = (t<0)
  if (sz[2]==1) chk = chk + rep( kappa <= 0, N ) else
    chk = chk + (kappa <= 0)
  if (sz[3]==1) chk = chk + rep( xi <= 0, N ) else
    chk = chk + (xi <= 0)
  if (sz[4]==1) chk = chk + rep( tau < 0, N ) else
    chk = chk + (tau < 0)

  # Initialize output
  out = matrix( NA, N, 3 )

  # Shift response times
  ts = t - tau

  # Set values of ts < 0 to 0
  ts[ ts < 0 ] = 0
  chk = chk + (ts==0)

  ### Derivative for kappa ###

  # f = log( ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 ) )
  # f = log( kappa ) - log( sqrt( 2*pi*ts^3) ) + ( (-.5/ts)*( kappa - (ts*xi) )^2 )
  # Let K1 = log( sqrt( 2*pi*ts^3) ); K2 = -.5/ts; K3 = ts*xi
  # h = kappa - K3
  # f = log( kappa ) - K1 + K2*( h )^2
  # h' = 1
  # f' = (1/kappa) + 2*K2*( h )*h'
  # f' = (1/kappa) + 2*K2*( kappa - K3 )
  # To check in Wolfram Alpha: log( x ) - a + b*( x - c )^2

  # Define constants
  K2 = -.5/ts
  K3 =  ts*xi

  # Calculate derivative
  der_kappa = (1/kappa) + 2*K2*( kappa - K3 )
  if (tran) der_kappa = der_kappa*kappa
  der_kappa[ chk > 0 ] = NA

  ### Derivative for xi ###

  # log( ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 ) )
  # [ log( kappa/sqrt( 2*pi*ts^3) ) ] + ( [-.5/ts]*( [ kappa ] - [ts]*xi )^2 )
  # h = kappa - ts*xi
  # f = K1 + K2*( h )^2
  # h' = -ts
  # f' = 2*K2*( h )*h'
  # f' = 2*K2*( kappa - ts*xi )*(-ts)
  # f' = -2*K2*ts*( kappa - ts*xi )
  # To check in Wolfram Alpha: a + b*( c - d*x )^2

  # Define constants
  K2 = -.5/ts

  # Calculate derivative
  der_xi = -2*K2*ts*(kappa-ts*xi)
  if (tran) der_xi = der_xi*xi
  der_xi[ chk > 0 ] = NA

  ### Derivative for tau ###

  # l = log( ( kappa/sqrt( 2*pi*(t-tau)^3) )*exp( (-.5/(t-tau)) * (kappa - (t-tau)*xi)^2 ) )
  # f = log( kappa/sqrt( 2*pi*(t-tau)^3) )
  # g = (-.5/(t-tau))*(kappa - (t-tau)*xi)^2
  # l = f + g
  # Part 1
  # f = log( [ kappa/sqrt( 2*pi ) ] ) + log( (t-tau)^(-3/2) )
  # Let K1 = kappa/sqrt( 2*pi )
  # f = log(K1) + log( (t-tau)^(-3/2) )
  # f' = 3/( 2*(t - tau) )
  # To check in Wolfram Alpha: log(b) + log( (a-x)^(-3/2) )
  # Part 2
  # A) Focus on squared term
  # (kappa - (t-tau)*xi)^2
  # ( (kappa - t*xi) + tau*xi )^2
  # Let K2 = kappa - t*xi
  # Multiply out
  # K2^2 + 2*K2*xi*tau + (xi^2)*tau^2
  # B) Include numerator from outside term
  # -.5*K2^2 - K2*xi*tau - .5*(xi^2)*tau^2
  # Let K3 = -.5*K2^2; K4 = K2*xi; K5 = .5*(xi^2)
  # g = (K3 - K4*tau - K5*tau^2)/(t-tau)
  # g' = ( K3 + K5*tau^2 - t*(K4 + 2*K5*tau) )/(t-tau)^2
  # To check in Wolfram Alpha: (b - c*x - d*x^2)/(a-x)
  # l' = f' + g'

  # Define constants
  # K1 =  kappa/sqrt( 2*pi )
  K2 = kappa - t*xi
  K3 = -.5*K2^2 # b
  K4 = K2*xi # c
  K5 = .5*(xi^2) # d

  # Simplify equation
  fd = 3/( 2*(t - tau) )
  gd = ( K3 + K5*tau^2 - t*(K4 + 2*K5*tau) )/(t-tau)^2

  # Calculate derivative
  der_tau = fd + gd
  if (tran) der_tau = der_tau*minVal*exp(orig_tau)/( exp(orig_tau) + 1 )^2
  der_tau[ chk > 0 ] = NA

  # Store results into matrix
  out[,1] = der_kappa
  out[,2] = der_xi
  out[,3] = der_tau
  colnames( out ) = c( 'dkappa','dxi','dtau' )

  if (summation) out = colSums( out )

  out
}

# Lookup - 05
msigDerivative = function( t, kappa, xi, tau, lambda,
                           alpha = min(t), beta = max(t),
                           summation = T, tran = F ) {
  # Purpose:
  # A function that calculates the first partial derivaties
  # for the log-likelihood of the mixture of the shifted inverse
  # gaussian and the uniform distribution
  # Arguments:
  # t         - A vector of response times
  # kappa     - A threshold value (kappa > 0)
  # xi        - A drift rate (xi > 0)
  # tau       - A residual latency (tau > 0)
  # lambda    - A mixture probability (0 <= lambda <= 1)
  # alpha     - The lower boundary for the uniform distribution
  # beta      - The upper boundary for the uniform distribution
  # summation - If true, sums the derivatives across observations
  #             for each parameter
  # tran      - If true, transforms the inputted parameters
  # Returns:
  # A vector of 4 values, each parameter's partial derivative for
  # the sum of the log-likelihoods. If 'summation' is set to false
  # returns a matrix with the derivative for each individual
  # log-likelihood.

  # Select longest vector
  sz = c( length(t), length(kappa),
          length(xi), length(tau), length(lambda) )
  N = max( sz )

  # Transform parameters
  if (tran) {
    kappa = exp(kappa)
    xi = exp(xi)
    tau = exp(tau)
    orig_lambda = lambda
    lambda = logistic(orig_lambda)
  }

  # Check for inadmissable values
  chk = (t<0)
  if (sz[2]==1) chk = chk + rep( kappa <= 0, N ) else
    chk = chk + (kappa <= 0)
  if (sz[3]==1) chk = chk + rep( xi <= 0, N ) else
    chk = chk + (xi <= 0)
  if (sz[4]==1) chk = chk + rep( tau < 0, N ) else
    chk = chk + (tau < 0)
  if (sz[5]==1) chk = chk + rep( ( (lambda < 0) | (lambda > 1) ), N ) else
    chk = chk + ( (lambda < 0) | (lambda > 1) )

  # Initialize output
  out = matrix( NA, N, 4 )

  if ( sum(chk) == 0 ) {

    # Shift response times
    ts = t - tau
    # If less than 0, set to 0
    ts[ ts < 0 ] = 0

    # Constant for uniform distribution
    U = dunif( t, alpha, beta )

    ### Derivative for kappa ###
    # f = log( lambda*( ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 ) ) + (1-lambda)*(1/(beta - alpha)) )
    # Let K1 = (1-lambda)*(1/beta-alpha); K2 = (lambda/sqrt( 2*pi*ts^3))
    # Let K3 = -.5/ts; K4 = ts*xi
    # f = log( kappa*K2*exp( K3*(kappa - K4)^2 ) + K1 )

    # Define constants
    K1 = (1-lambda)*U; K2 = (lambda/sqrt( 2*pi*ts^3));
    K3 = -.5/ts; K4 = ts*xi;

    # Calculate derivative
    num = ( K2*exp(K3*(-K4 + kappa)^2) +
              2*K2*K3*exp(K3*(-K4 + kappa)^2)*kappa*(-K4 +
                                                       kappa))
    denom = (K1 + K2*exp(K3*(-K4 + kappa)^2)*kappa)
    der_kappa = num/denom;
    if (tran) der_kappa = der_kappa*kappa

    ### Derivative for xi ###
    # f = log( lambda*( ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 ) ) + (1-lambda)*(1/(beta - alpha)) )
    # Let K1 = (1-lambda)*(1/beta-alpha); K2 = lambda*kappa/sqrt( 2*pi*ts^3 )
    # Let K3 = -.5/ts;
    # f = log( K2*exp( K3*(kappa - ts*xi)^2 ) + K1 )

    # Define constants
    K1 = (1-lambda)*U;
    K2 = lambda*kappa/sqrt( 2*pi*ts^3 );
    K3 = -.5/ts;

    # Calculate derivative
    num = -2*K2*K3*ts*(kappa - ts*xi)*exp( K3*(kappa - ts*xi)^2 )
    denom = K1 + K2*exp( K3*(kappa - ts*xi)^2 )
    der_xi = num/denom
    if (tran) der_xi = der_xi*xi

    ### Derivative for tau ###

    # f = log( lambda*( ( kappa/sqrt( 2*pi*(t - tau)^3) )*exp( (-.5/(t - tau)) * (kappa - (t - tau)*xi)^2 ) ) + (1-lambda)*(1/(beta - alpha)) )
    # Let K1 = (1 - lambda)/(beta-alpha); K2 = lambda*kappa/sqrt( 2*pi );
    # Let K3 = kappa - t*xi;
    # f = log( K2*( (t - tau)^(-3/2) )*exp( ( -.5/(t - tau) )*( K3 + tau*xi )^2 ) + K1 )

    # Define constants
    K1 = (1 - lambda)*U;
    K2 = lambda*kappa/sqrt( 2*pi );
    K3 = kappa - t*xi;

    # Calculate derivative
    num = -.5*K2*( t*(2*K3*xi + 2*tau*xi^2 - 3 ) + K3^2 - tau*(tau*xi^2 - 3) )
    denom = ( (t-tau)^2 )*( K1*( (t - tau)^(3/2) )*exp(.5*( (K3 + tau*xi)^2 )/(t-tau) ) + K2 )
    der_tau = num/denom
    if (tran) der_tau = der_tau*tau

    ### Derivative for lambda ###

    # f = log( lambda*( ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 ) ) + (1-lambda)*(1/(beta - alpha)) )
    # Let K2 = (1/beta-alpha); K1 = ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 )
    # f = log( lambda*K1 + (1-lambda)*K2 )

    # Define constants
    K1 = ( kappa/sqrt( 2*pi*ts^3) )*exp( (-.5/ts) * (kappa - ts*xi)^2 );
    K2 = U;

    # Calculate derivative
    der_lambda = (K1 - K2)/(K1*lambda + K2*(-lambda) + K2)
    der_lambda[ ts == 0 ] = 1/(lambda-1)
    if (tran) der_lambda = der_lambda*exp(orig_lambda)/(exp(orig_lambda) + 1)^2

    # Output
    out[,1] = der_kappa
    out[,2] = der_xi
    out[,3] = der_tau
    out[,4] = der_lambda

    # Correct for inadmissable values
    out[ ts == 0, 1:3 ] = 0

  }

  colnames( out ) = c( 'dkappa','dxi','dtau','dlambda')

  if (summation) out = colSums( out )

  out
}

# Lookup - 06
AIC_c = function( logLik, K, n ) {
  # Purpose:
  # Calculates Akaike's Information Criterion corrected for small
  # samples.
  # Arguments:
  # logLik - The summed set of log-likelihoods for a data set
  # K      - The number of free parameters in the model
  # n      - The sample size (i.e. the number of observations)
  # Returns:
  # Forthcoming
  # References:
  # Forthcoming

  -2*logLik + 2*K + ( 2*K*(K+1) )/(n-K-1)

}

mleSE = function(fit) {

  fisher_info = solve(-1*fit$hessian)
  sigma_hat = sqrt(diag(fisher_info))

  return( sigma_hat )
}

# BIC = function( logLik, K, N ) {
# }

mleWrapper = function( rt,
                       type,
                       mle_f,
                       st_f,
                       startRange,
                       gr_f = NULL,
                       method = 'BFGS',
                       emStop = 20,
                       maxit = 5000,
                       nRep = 5,
                       best = T,
                       alpha = .95,
                       progress = F) {

  # Record starting time
  startTime = Sys.time()

  # Variables to store results
  logLikSum = rep(NA,nRep)
  start_val = st_f(type,startRange[[1]],startRange[[2]]) # Determine initial length
  parNames = names(start_val)
  start_val = matrix(NA,nRep,length(start_val))
  allEst = matrix(NA,nRep,ncol(start_val))
  colnames(allEst) = colnames(parNames)
  aic = vector("numeric",nRep)
  results = c()
  mleVar = c()
  seChk = rep(NA,nRep)
  # Initialize lists
  for (i in 1:nRep) {
    results=c(results,list(NULL))
    mleVar = c(mleVar,list(NULL))
  }
  inc = 1

  # Loop through repetitions
  for (nr in 1:nRep) {

    # Generate starting values that don't produce invalid
    # log-likelihoods
    chk = -Inf
    em = 0
    while( chk == -Inf ) {
      startVal = st_f(type,startRange[[1]],startRange[[2]])
      chk = mle_f( startVal, rt, type )
      em = em + 1
      if (em == emStop) stop("No feasible starting values were found",
                             call. = FALSE)
    }

    # Set default value in case estimation fails
    res = NULL

    # Run estimation routine within a 'try' statement
    # so that function won't stop if there's an error
    res = try( optim( startVal, mle_f,
                      rt = rt,
                      type = type,
                      gr = gr_f,
                      method = method,
                      hessian = T,
                      control = list( fnscale = -1, maxit = maxit ) ),
               silent = T )

    # Save results if no error occurred
    if (is.list(res)) {
      results[[inc]] = res
      start_val[inc,]=startVal
      logLikSum[inc]=res$value
      aic[inc] = AIC_c( res$value, length( startVal ), length( rt ) )

      # Extract standard errors and calculate confidence
      # intervals
      est = res$par
      SE = NULL
      SE = try( mleSE( res ), silent = T )
      if ( is.numeric(SE) ) {
        CI = NULL
        CI = try( rbind( est - qnorm(1 - (1-alpha)/2)*SE,
                         est + qnorm(1 - (1-alpha)/2)*SE ), silent = T )
        # Transform parameters
        if (type==1) CI = exp(CI)
        if (type==2) {
          CI[,1:2] = exp( CI[,1:2] )
          CI[,3] = logistic(CI[,3])*min(rt)
        }
        if (type==3) {
          CI[,1:3] = exp( CI[,1:3] )
          CI[,4] = logistic(CI[,4])
        }
        names(SE) = parNames
        colnames(CI) = parNames
        mleVar[[inc]] = list( SE = SE, CI = CI )
        seChk[inc] = 1;
      }
      # Transform estimates for the best-fitting parameters
      if (type==1) est = exp( est )
      if (type==2) {
        est[1:2] = exp(est[1:2])
        est[3] = logistic(est[3])*min(rt)
      }
      if (type==3) {
        est[1:3] = exp(est[1:3])
        est[4] = logistic(est[4])
      }
      allEst[inc,] = est

      inc = inc + 1
    }

  }

  # Remove any missing values
  chk = which(is.na(logLikSum))
  chk = unique( chk, which(is.na(seChk)) )
  if (length(chk)>0) {
    results = results[-chk]
    start_val = start_val[-chk,]
    logLikSum = logLikSum[-chk]
    allEst = allEst[-chk,]
    aic = aic[-chk]
    mleVar = mleVar[-chk]
  }

  if ( length(chk) == nRep ) stop("No viable estimates were found.",
                                  call. = FALSE)

  # Return results
  if (best) {

    sel = min( which( logLikSum == max( logLikSum ) ) )
    par = allEst[sel,]
    names(par) = parNames
    SE = mleVar[[sel]]$SE
    CI = mleVar[[sel]]$CI
    runTime = Sys.time() - startTime

    return( list(
      par = par,
      SE = SE,
      CI = CI,
      MLE = results[[sel]],
      AICc = aic[sel],
      startVal = start_val[sel,],
      Time = runTime
    ) )
  } else {
    return( list(
      MLE = results,
      logLikSum = logLikSum,
      AICc = aic,
      startVal = startVal,
      Time = runTime
    ) )
  }

}

# Likelihood function to pass to mleWrapper
mle_f = function(par,rt,type) {

  logLik = -Inf
  if (type==1) {
    par[1:2] = exp(par[1:2])
    logLik = dig( rt, par[1], par[2], ln = T, summation = T )
  }
  if (type==2) {
    par[1:2] = exp(par[1:2])
    par[3] = logistic(par[3])*min(rt)
    logLik = dsig( rt, par[1], par[2], par[3], ln = T, summation = T )
  }
  if (type==3) {
    par[1:3] = exp(par[1:3])
    par[4] = logistic(par[4])
    logLik = dmsig( rt, par[1], par[2], par[3], par[4], ln = T, summation = T )
  }

  return( logLik )
}

# Gradient function to pass to mleWrapper
gr_f = function(par,rt,type) {

  if (type==1) return( igDerivative( rt, par[1], par[2], tran = T ) )
  if (type==2) return( sigDerivative( rt, par[1], par[2], par[3], tran = T ) )
  if (type==3) return( msigDerivative( rt, par[1], par[2], par[3], par[4], tran = T ) )

}

# Starving values function to pass to mleWrapper
st_f = function(type,lowerBoundary,upperBoundary) {

  # Inverse gaussian
  if (type==1) {
    par = runif(2,lowerBoundary[1:2],upperBoundary[1:2])
    par = log( par )
    names(par) = c('kappa','xi')
  }
  if (type==2) {
    par = runif(3,lowerBoundary[1:3],upperBoundary[1:3])
    par[1:2] = log( par[1:2] )
    par[3] = logit( par[3] )
    names(par) = c('kappa','xi','tau')
  }
  if (type==3) {
    par = runif(4,lowerBoundary[1:4],upperBoundary[1:4])
    par[1:3] = log( par[1:3] )
    par[4] = logit( par[4] )
    names(par) = c('kappa','xi','tau','lambda')
  }

  return( par )
}
