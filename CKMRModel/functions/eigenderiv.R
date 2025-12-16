"eigenderiv" <- function( Mfun, theta0, ..., nits=10, quietly=TRUE){
  ## Iterate power-method to get lead evec of a matrix with lead eval==1
  ## It works, but NB isume==1 and D_isume==0, so a more efficient approach is to solve...
  
  require( vecless, quietly= TRUE)
  
  e1fun <- function( theta){
    MM <- Mfun( theta, ...)
    EE1 <- eigen( MM)$vectors[,1]
    EE1 <- EE1 / sum( EE1)
    return( EE1)
  }
  
  M <- Mfun( theta0, ...)
  n <- nrow( M)
  stopifnot( n==ncol( M))
  
  D_M <- numderiv( Mfun, theta0, eps=1e-8, ...)[,,1]
  e1 <- e <- e1fun( theta0)
  D_e <- 0*e # to start
  
  for( it in 1:nits){
    etemp[ i]:= M[i,j] %[j]% e[ j];
    D_etemp[ i]:= D_M[i, j] %[j]% e[ j] + M[ i, j] %[j]% D_e[ j];
    sume:= SUM_ %[i]% etemp[ i];
    D_sume:= SUM_ %[i]% D_etemp[ i];
    isume <- 1 / sume;
    D_isume <- -sqr( isume) * D_sume;
    if( !quietly){
      scatn( 'It %i: isume=%6.3f, D_isume=%6.3f', it, isume, D_isume)
    }
    e[ i]:= isume * etemp[ i];
    D_e[ i]:= D_isume * etemp[ i] + isume * D_etemp[ i];
  };
  
  # D_e[ i]:= D_M[i, j] %[j]% e[ j] + M[ i, j] %[j]% D_e[ j];
  # De = DM %*% e + M %*% De
  # => De = -(M-I)^-1 * (DM %*% e) but this is singular since M has one eval of +1
  # so need to impose sum-to-one constraint also
  
  X <- rbind( (M - diag( n))[-n,], 1)
  D_e_direct = solve( X, c( -D_M[-n,] %**% e, 0))
  
  D_e_check <- numderiv( e1fun, theta0)
  returnList( e1, e, D_e, D_e_check, D_e_direct);
}