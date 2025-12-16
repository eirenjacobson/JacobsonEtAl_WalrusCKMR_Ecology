"ederiv2" <- function( Mfun, theta0, ...){
  ## Eigenderiv of lead evec *and* lead eval, for general matrix-valued function 'Mfun()' of 'theta0' (and perhaps other fixed parameters, in '...')
  
  # Lead evec, then its eval (called "lambda")
  elfun <- function( theta){
    MM <- Mfun( theta, ...)
    eigM <- eigen( MM)
    EE1 <- eigM$vectors[,1]
    EE1 <- EE1 / sum( EE1)
    return( c( EE1, lambda=eigM$values[1]))
  }
  
  M <- Mfun( theta0, ...)
  n <- nrow( M)
  stopifnot( n==ncol( M))
  
  D_M <- numderiv( Mfun, theta0, eps=1e-8, ...)[,,1]
  el <- elfun( theta0)
  e <- head( el, -1)
  lambda <- tail( el, 1)
  
  X <- cbind( rbind( M - lambda * diag( n), 1), c( -e, 0))
  D_e_direct = solve( X, c( -D_M %**% e, 0))
  
  D_e_check <- numderiv( elfun, theta0)
  returnList( e, lambda, D_e_check, D_e_direct);
}