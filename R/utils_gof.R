#' Calculate the estimate of W_{n}(u) process over a grid of pit values.
#'
#' @param S score matrix with n rows and p columns.
#' @param FI Fisher information matrix with p rows and p columns.
#' @param pit a numeric vector of PIT values.
#'
#' @return a matrix with n rows and n columns
#'
#' @noRd
calculateWnuhat = function(S, FI, pit){

  # Find the number of observations
  n       <- nrow(S)

  # Compute the values of indicator function over pit values
  ind     <- outer(pit, pit, '<=')

  # Estimate the Psi vector by covariance of indicator and score
  Psi_hat <- ( (n-1) * cov(ind, S) ) / n

  # Compute the covariance function over the grid (pit)
  Mat     <- ( ind - S %*% solve( FI ) %*% t(Psi_hat) )
  colnames(Mat) <- 1:n
  rownames(Mat) <- 1:n

  # Return the matrix
  return(Mat)

}


#' Calculate the estimate of W_{n}(u) process on a grid of equally spaced values over (0,1) interval
#'
#' @param S score matrix with n rows and p columns.
#' @param FI Fisher information matrix with p rows and p columns.
#' @param pit a numeric vector of PIT values.
#' @param M number of equally spaced values over (0,1) interval.
#'
#' @return a matrix with n rows and M columns
#'
#' @noRd
calculateWnuhat_manualGrid = function(S, FI, pit, M){

  # Find the number of observations
  n       <- nrow(S)

  # Create a grid points over (0,1) interval. Note that we need to add a small number to the edges of interval.
  # The covariance function does not define on the edges. We used 1e-5 for the epsilon.
  epsilon <- 1e-5
  gridpts <- seq(0 + epsilon, 1 - epsilon, length = M)

  # Compute the values of indicator function over pit values
  ind     <- outer(pit, gridpts, '<=')

  # Estimate the Psi vector by covariance of indicator and score
  Psi_hat <- ( (n-1) * cov(ind, S) ) / n

  # Compute the estimate of W_{n}(u) process over the grid (pit)
  Mat     <- ( ind - S %*% solve( FI ) %*% t(Psi_hat) )
  colnames(Mat) <- paste0('u', 1:M)
  rownames(Mat) <- 1:n

  return(Mat)

}


#' Compute Eigenvalues of the covariance matrix
#'
#' @param S score matrix with n rows and p columns.
#' @param FI Fisher information matrix with p rows and p columns.
#' @param pit a numeric vector of PIT values.
#' @param me the goodness-of-fit statistic, Cramer-von-Mises or Anderson-Darling.
#'
#' @return a numeric vector of Eigenvalues.
#'
#' @noRd
getEigenValues = function(S, FI, pit, me){

  # Find the number of observations
  n       <- nrow(S)

  # Find the number of parameters
  p       <- ncol(S)

  # Compute the estimate of W_{n}(u) process over a grid of pit values.
  Mat     <- calculateWnuhat(S, FI, pit)

  # Compute the covariance of the estimate of W_{n}(u) process and adjust for the sample size.
  W       <- var(Mat)
  #W       <- ( (n-1) * W ) / (n)
  W       <- ( (n-1) * W ) / (n-p-1)


  # Compute b vector to adjust W matrix
  pit <- sort(pit)
  l   <- length(pit)
  b   <- numeric()
  for( j in 1:l ){

    if( j == 1){
      b[j] <- pit[2]
    }

    if( j == 2){
      b[j] <- pit[3] - pit[1]
    }

    if( (j>2) & (j<=(n-1)) ){
      b[j] <- pit[j+1] - pit[j-1]
    }

    if( j == n){
      b[j] <- 1 - pit[n-1]
    }

  }
  b   <- b / 2
  b   <- b / sum(b)

  # Create diagonal matrix with b vector elements
  Q <- diag( sqrt(b) )

  W <- Q %*% W %*% Q


  # Compute the Eigenvalues of the covariance matrix depending on the goodness-of-fit statistic
  if( me == 'cvm' ){
    #ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / length(pit)
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
    return(ev)
  }

  if( me == 'ad'){
    adj.value <- sqrt( outer( pit * (1- pit), pit * (1- pit) ) )
    W       <- W / adj.value
    #ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / length(pit)
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
    return(ev)
  }

}


#' Compute Eigenvalues of the covariance matrix
#'
#' @param S score matrix with n rows and p columns.
#' @param FI Fisher information matrix with p rows and p columns.
#' @param pit a numeric vector of PIT values.
#' @param M number of equally spaced values over (0,1) interval.
#' @param me the goodness-of-fit statistic, Cramer-von-Mises or Anderson-Darling.
#'
#' @return a numeric vector of Eigenvalues.
#'
#' @noRd
getEigenValues_manualGrid = function(S, FI, pit, M, me){

  # Find the number of observations
  n       <- nrow(S)

  # Find the number of parameters
  p       <- ncol(S)

  # Compute the estimate of W_{n}(u) process over a grid of values equally spaced over (0,1).
  Mat     <- calculateWnuhat_manualGrid(S, FI, pit, M)

  # Compute the covariance of the estimate of W_{n}(u) process and adjust for the sample size.
  W       <- var(Mat)
  W       <- ( (n-1) * W ) / (n-p-1)


  # Compute b vector to adjust W matrix
  pit <- sort(pit)
  l   <- length(pit)
  b   <- numeric()
  for( j in 1:l ){

    if( j == 1){
      b[j] <- pit[2]
    }

    if( j == 2){
      b[j] <- pit[3] - pit[1]
    }

    if( (j>2) & (j<=(n-1)) ){
      b[j] <- pit[j+1] - pit[j-1]
    }

    if( j == n){
      b[j] <- 1 - pit[n-1]
    }

  }
  b   <- b / 2
  b   <- b / sum(b)

  # Create diagonal matrix with b vector elements
  Q <- diag( sqrt(b) )

  W <- Q %*% W %*% Q

  # Compute the Eigenvalues of the covariance matrix depending on the goodness-of-fit statistic
  if( me == 'cvm' ){
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
    return(ev)
  }

  if( me == 'ad'){
    s <- 1:M
    s <- s / (M + 1)
    adj.value <- sqrt( outer( s*(1-s) , s*(1-s) ) )
    W <- W / adj.value
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
    return(ev)
  }

}


#' Compute the Cramer-von-Mises statistics
#'
#' @param x a numeric vector of pit values
#'
#' @return a numeric value, CvM statistics
#'
#' @noRd
getCvMStatistic = function(x){

  n   <- length(x)
  Z   <- sort(x)
  a   <- seq(from = 1, to = 2*n-1, by = 2)
  a   <- a/(2*n)
  res <- sum( (Z - a)^2 ) + 1/(12*n)
  return(res)

}


#' Compute the Anderson-Darling statistics
#'
#' @param x a numeric vector of pit values
#'
#' @return a numeric value, AD statistics
#'
#' @noRd
getADStatistic = function(x){

  n   <- length(x)
  Z   <- sort(x)
  a   <- seq(from = 1, to = n, by = 1)
  a   <- (2 * a) - 1
  S   <- sum( a * ( log(Z) + log( 1 - rev(Z) ) ) )
  res <- (-S/n) - n
  return(res)

}


#' Calculate p-value for the goodness-of-fit test.
#'
#' @description The calculation of p-value is done by the aid of \code{CompQuadForm} package and \code{Farebrother} function.
#' It computes Pr(Q > u) where Q is a sum of squared Chi-quared variables weighted by non-zero Eigenvalues.
#'
#' @param u  a numeric value, Cramer-von-Mises statistic or Anderson-Darling statistic.
#' @param eigen a numeric vector of Eigenvalues.
#'
#' @return p-value
#'
#' @noRd
getpvalue = function(u, eigen){

  # set_1 is from i=1 to J1
  eigen <- eigen[eigen >= 1e-15]
  indx  <- which( eigen >= eigen / 2000)
  set_1 <- eigen[indx]

  # set_2 is from i=J1+1 to J1+J2
  uc      <- eigen[1] / 2000
  indx    <- which( (eigen >= 1e-15) & (eigen <= uc) )
  set_2   <- eigen[indx]

  # Compute the lower bound for the probability of Pr(Q > u)
  LB <- getLowerBoundForpvalue(statistic = u, lambda = eigen)

  # Compute the upper bound for the probability of Pr(Q > u)
  UB <- getUpperBoundForpvalue(statistic = u, lambda = eigen)

  if( LB >= 1e-7){
    # Compute p-value by imhof method
    pvalue <- CompQuadForm::imhof(q = u - sum(set_2), lambda = set_1)$Qq

    if( (pvalue < LB) & (pvalue > UB) ){
       # Imhof method failed to produce a valid pvalue. Compute p-value by farebrother method instead
       pvalue <- CompQuadForm::farebrother(q = u - sum(set_2), lambda = set_1)$Qq
    }
  }


  if( (1e-10 <= LB) & (LB < 1e-7) ){
    # Compute p-value by farebrother method
    pvalue <- CompQuadForm::farebrother(q = u - sum(set_2), lambda = set_1)$Qq
  }


  if( LB < 1e-10 ){
    # Compute p-value by farebrother method
    pvalue <- CompQuadForm::farebrother(q = u - sum(set_2), lambda = set_1)$Qq
  }

  # Check if the computed p-value by CompQuadForm package falls between lower and upper bound.
  if( (pvalue >= LB) & (pvalue <= UB) ){
    return(pvalue)
  }else{
    message(paste0('CompQuadForm failed to generate a valid p-value. The p-value lies between ', LB, ' and ', UB))
    return(pvalue)
  }

}


#' Compute lower bound for p-value
#'
#' @param statistic a numeric value, Cramer-von-Mises statistics or Anderson-Darling statistics
#' @param lambda a numeric vector containing Eigenvalues
#'
#' @return a numeric value
#'
#' @noRd
getLowerBoundForpvalue = function(statistic, lambda){

  term1 <- integrate(f = integrandForLowerBound, lower = 0, upper = statistic/lambda[1], ST = statistic, EV = lambda)$value
  term2 <- pchisq(q = statistic/lambda[1], df = 1, lower.tail = FALSE)
  return(term1 + term2)

}


#' Compute upper bound for p-value.
#'
#' @param statistic a numeric value, Cramer-von-Mises statistics or Anderson-Darling statistics.
#' @param lambda a numeric vector containing Eigenvalues.
#' @param tol the tolerance used to find the solution. The default value is 1e-10.
#' @param max.iter the maximum number of iteration to find the solution. The default value is 50.
#'
#' @return a numeric value
#'
#' @noRd
getUpperBoundForpvalue = function(statistic, lambda, tol = 1e-10, max.iter = 50){

  # Check to see if the root is at t = 0 and return the root and upper bound for pvalue.
  if( sum(lambda) >= statistic ){
    root <- 0
    fval <- exp(-root * statistic - 0.5 * sum(log(1-2*root*lambda)))
    return(UB_pvalue = fval)
  }

  # Define the interval to search for the root using bisection method.
  t_LB <- 0
  t_UB <- 1/(2*max(lambda)) - 1e-6

  # Evaluate values at the edges of interval
  f.LB  <- sum( lambda / (1-2*lambda*t_LB) ) - statistic
  f.UB  <- sum( lambda / (1-2*lambda*t_UB) ) - statistic

  # m counts the number of iteration and if max.iter reached, the process stops.
  m    <- 0

  # Work on solution as long as t_UB and t_LB are not close enough.
  while( abs(t_UB - t_LB) > tol) {

    # Increase iteration
    m <- m + 1

    # Check if the maximum iteration reached
    if (m > max.iter) {
      # maybe we should return CompQuadForm value if the maximum reached? Check with Richard
      message('maximum iteration reached.')
      warning(paste0('last root was:', root, ' and upper boud of ', fval))
      break
    }

    # Find the center of interval [t_LB,t_UB] and calculate the value of function at this middle point.
    t_mid <- (t_LB + t_UB)/2
    f.mid <- sum( lambda / (1-2*lambda*t_mid) ) - statistic

    #
    if (f.UB * f.mid < 0) {
      t_LB <- t_mid
      f.LB <- f.mid
    }else {
      t_UB <- t_mid
      f.UB <- f.mid
    }

  }

  root <- (t_LB + t_UB)/2
  fval <- exp(-root * statistic - 0.5 * sum(log(1-2*root*lambda)))

  return(UB_pvalue = fval)

}



#' Integrand function to compute lower bound
#'
#' @param t function argument
#' @param ST a numeric value, Cramer-von-Mises statistics or Anderson-Darling statistics.
#' @param EV a numeric vector containing Eigenvalues.
#'
#' @return a numeric value
#'
#' @noRd
integrandForLowerBound = function(t, ST, EV){
  a   <- EV[2] / EV[1]
  b   <- ST / EV[1]
  res <- pchisq(q = (b-t)/a, df = 1, lower.tail = FALSE) * dchisq(x = t, df = 1)
  return(res)
}
