"ldpois" <- function( o, e){
  ## sum log( Poisson prob of o given means e)
  # Might be more efficient than dpois( log=T)
  # coz no log-gamma, and only "a few" NZ obs needing logs
  
  o <- as.vector( o)
  e <- as.vector( e)
  # o %*% log( e+(o==0)*(e==0)) - sum( e) # do'em all version
  nz <- which( o>0)
  o[ nz] %*% log( e[ nz]) - sum( e)
}