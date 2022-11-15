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

  Mat     <- calculateWnuhat(S, FI, pit)
  W       <- var(Mat)

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
  colnames(Mat) <- 1:n
  rownames(Mat) <- 1:n

  return(Mat)

}


getEigenValues_manualGrid = function(S, FI, pit, M, me){

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

  pvalue  <- CompQuadForm::imhof(q = u, lambda = eigen)$Qq

  if( pvalue < 0 ){
    pvalue <- 2 * pnorm( - sqrt( u / max(eigen) ) )
    warning('CompQuadForm generated a negative pvalue. The pvalue replaced by a lower bound.')
  }

  if( pvalue > 1 ){
    pvalue <- 1
    warning('CompQuadForm generated pvalue > 1. The pvalue replaced by an upper bound.')
  }

  return(pvalue)
}

