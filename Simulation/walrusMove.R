walrusMove <- function(indiv = makeFounders(), moveMat) {
    if(nrow(moveMat) != ncol(moveMat)) warning("movement matrix is not square")
    if(nrow(moveMat) > length(unique(indiv[,7]))) warning("One or more stocks start empty")

    oldStock <- indiv[,7]
    newStock <- matrix(data = NA, nrow = nrow(indiv))
    for(i in 1:length(newStock)) {
        newStock[i] <- sample(1:nrow(moveMat), 1, TRUE,
                              prob = moveMat[indiv[i,7],])
    }
    indiv[,7] <- newStock
    indiv[,7][!is.na(indiv[,6])] <- oldStock[!is.na(indiv[,6])]
    indiv[,7] <- as.integer(indiv[,7])
    return(indiv)
}