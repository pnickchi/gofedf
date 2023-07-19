calculateWnuhat = function(S, FI, pit){

  n       <- nrow(S)
  ind     <- outer(pit, pit, '<=')
  Psi_hat <- ( (n-1) * cov(ind, S) ) / n

  Mat     <- ( ind - S %*% solve( FI ) %*% t(Psi_hat) )
  colnames(Mat) <- 1:n
  rownames(Mat) <- 1:n

  return(Mat)

}


getEigenValues = function(S, FI, pit, me){

  n       <- nrow(S)
  Mat     <- calculateWnuhat(S, FI, pit)
  W       <- var(Mat)
  W       <- ( (n-1) * W ) / n

  if( me == 'cvm' ){
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / length(pit)
    return(ev)
  }

  if( me == 'ad'){
    adj.value <- sqrt( outer( pit * (1- pit), pit * (1- pit) ) )
    W <- W / adj.value
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / length(pit)
    return(ev)
  }

}


calculateWnuhat_manualGrid = function(S, FI, pit, M){

  n       <- nrow(S)
  epsilon <- 1e-5
  gridpts <- seq(0 + epsilon, 1 - epsilon, length = M)
  ind     <- outer(pit, gridpts, '<=')
  Psi_hat <- ( (n-1) * cov(ind, S) ) / n

  Mat     <- ( ind - S %*% solve( FI ) %*% t(Psi_hat) )
  colnames(Mat) <- paste0('u', 1:M)
  rownames(Mat) <- 1:n

  return(Mat)

}


getEigenValues_manualGrid = function(S, FI, pit, M, me){

  n       <- nrow(S)
  Mat     <- calculateWnuhat_manualGrid(S, FI, pit, M)
  W       <- var(Mat)

  if( me == 'cvm' ){
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / M
    return(ev)
  }

  if( me == 'ad'){
    s <- 1:M
    s <- s / (M + 1)
    adj.value <- sqrt( outer( s*(1-s) , s*(1-s) ) )
    W <- W / adj.value
    ev      <- eigen(W, symmetric = TRUE, only.values = TRUE)$values / M
    return(ev)
  }

}


getCvMStatistic = function(x){

  n <- length(x)
  Z <- sort(x)
  a <- seq(from = 1, to = 2*n-1, by = 2)
  a <- a/(2*n)
  res <- sum( (Z - a)^2 ) + 1/(12*n)
  return(res)

}


getADStatistic = function(x){

  n <- length(x)
  Z <- sort(x)
  a <- seq(from = 1, to = n, by = 1)
  a <- (2 * a) - 1
  S <- sum( a * ( log(Z) + log( 1 - rev(Z) ) ) )
  res <- (-S/n) - n
  return(res)

}


getpvalue = function(u, eigen){

  LB <- getLowerBoundForpvalue(statistic = u, lambda = eigen)

  UB <- getUpperBoundForpvalue(statistic = u, lambda = eigen)

  pvalue <- CompQuadForm::farebrother(q = u, lambda = eigen)$Qq

  if( (pvalue >= LB) & (pvalue <= UB) ){
    return(pvalue)
  }else{
    warning(paste0('CompQuadForm failed to generate a valid p-value. The p-value lies between ', LB, ' and ', UB))
    return(pvalue)
  }

}


getLowerBoundForpvalue = function(statistic, lambda){

  term1 <- integrate(f = integrandForLowerBound, lower = 0, upper = statistic/lambda[1], ST = statistic, EV = lambda)$value
  term2 <- pchisq(q = statistic/lambda[1], df = 1, lower.tail = FALSE)
  return(term1 + term2)

}

integrandForLowerBound = function(t, ST, EV){
  # ST is the statistic from the sample
  # EV is the vector of Eigen values
  a   <- EV[2] / EV[1]
  b   <- ST / EV[1]
  res <- pchisq(q = (b-t)/a, df = 1, lower.tail = FALSE) * dchisq(x = t, df = 1)
  return(res)
}

getUpperBoundForpvalue = function(statistic, lambda, tol = 1e-10, max.iter = 50){

  # Check to see if the root is at t = 0 and return the root and upper bound for pvalue.
  if( sum(lambda) >= statistic ){
    root <- 0
    fval <- exp(-root * statistic - 0.5 * sum(log(1-2*root*lambda)))
    return(UB_pvalue = fval)
  }

  # Define the interval to search for the root
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
      warning('maximum iteration reached.')
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




