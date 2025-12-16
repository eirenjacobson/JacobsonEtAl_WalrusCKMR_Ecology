"add_recap_kin" <- function( origP, me_first, me_again){
  ## quickin starts with unique samples; if recaps, add recap kin too
  # Recaps might be in 1st or/AND 2nd col
  # First col first, then expand, then second col
  
  m <- match( origP[,1], me_first, 0)
  if( any( m)){
    newP <- cbind( me_again[ m[ m>0]], origP[m>0,2])
    origP <- rbind( origP, newP)
  }
  
  m <- match( origP[,2], me_first, 0)
  if( any( m)){
    newP <- cbind( origP[m>0,1], me_again[ m[ m>0]])
    origP <- rbind( origP, newP)
  }
  return( origP)
}