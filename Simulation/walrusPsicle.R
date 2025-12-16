# breeding cycle -- affects females only!
# ignores M and juv (or any others with stock = 0, as set by walrusFounders)
psicle <- function(indiv = makeFounders(), moveMat) {
  
  library(fastmatch)
  if(nrow(moveMat) != ncol(moveMat)) warning("movement matrix is not square")
  if(nrow(moveMat) > length(unique(indiv[,7]))) warning("One or more stocks start empty")
  
  # only consider females within breeding psicle (stock == 1, 2, 3)
  who <- which(indiv$Stock != 0 & is.na(indiv$DeathY))
  
  oldStock <- indiv[who,7]
  newStock <- indiv[,7] # template, so only F are affected
  for(i in who) {
    newStock[i] <- sample(1:nrow(moveMat), 1, TRUE,
                          prob = moveMat[indiv[i,7],])
  }
  indiv[,7] <- newStock
  indiv[,7] <- as.integer(indiv[,7])
  
  # calves of f who move into stage 2 must become independent
  mtb <- which(indiv$Stock == 2 & is.na(indiv$DeathY)) # mums to be
  pdy <- which(indiv$AgeLast %in% 0:3 & is.na(indiv$DeathY) & indiv$Mum != "founder") # possible dependent young

  dc <- fmatch(indiv$Me[mtb], indiv$Mum[pdy]) # length 302
  dc <- dc[!is.na(dc)]
  
  indiv$Independent[pdy][dc] <- 1
  
  # calves of moms who have died must be independent
  
  dm <- which(!is.na(indiv$DeathY))
  ic <- fmatch(indiv$Me[dm], indiv$Mum[pdy])
  
  indiv$Independent[pdy][ic] <- 1

  return(indiv)
}




