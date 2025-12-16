# helps avoid out-of-range lookups (for efficiency)
"clippo" <- function( x, LIM1, LIM2){
  if( missing( LIM2)) 
    pmax( pmin( x, max( LIM1)), min( LIM1))
  else {
    stopifnot( LIM2 >= LIM1)
    pmax( pmin( x, LIM2), LIM1)
  }
}