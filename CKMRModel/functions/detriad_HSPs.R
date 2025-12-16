"detriad_HSPs" <- function( HSPs, B){ 
## Requires only Maternal or only Paternal sibs
## HSPs is 2-col matrix, ordered as (firstborn, secondborn)
## B is birth-years of each sample (not just those in HSPs)
## Returns possibly-shortened version of HSPs, with only "independent" pairs

  multisib <- table( c( HSPs))
  
  # If an HS appears 2 times, there will be 3 HSPs where only 2 are indept; etc
  # and there will be 2 other HSs in that group.
    
  # Weed out superfluous HSPs where individuals appear > 1 time. 
  # Only use the first and second appearances of any HSP-member in sib-pair.
  # This is simply horrible... buckle up!

  imulti <- which( multisib > 1)
  if( length( imulti)){ # else we only have pairs, no triads
    nms <- as.integer( names( multisib)[ imulti])
    dunyet <- rep( FALSE, length( imulti))
    repeat{
      i <- match( FALSE, dunyet, 0) # first triadite that's not yet processed
      if( !i){
    break # all dun
      }
      
      # Find all animals which share i's mother, and store as 'sibs'
      # Assumes sib-finding is perfect--- not meant for real data with false-negs!

      sibs <- integer() # start with empty family
      new_sibs <- nms[ i] # start testing with that animal
      while( length( new_sibs)){
        sibs <- c( sibs, new_sibs)
        checky <- which( (HSPs[,1] %in% new_sibs) | 
            (HSPs[,2] %in% new_sibs))
        new_sibs <- unique( c( HSPs[ checky,])) %except% sibs
      }

      # Sort by birth-order (just in case...)
      sibs <- sibs[ order( B[ sibs])]

      # Drop HSPs unless they are successive births within this family
      sibgap <- match( HSPs[,2], sibs, 0) - match( HSPs[,1], sibs, 0)
      # ... # 0 is OK (different family); 1 is OK (successive)
      HSPs <- HSPs[ sibgap < 2, ] 

      dunyet[ match( sibs, nms)] <- TRUE
    } # while not all dunyet
  } # if any triads to thin

return( HSPs)    
}
