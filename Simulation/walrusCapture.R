#' identify genetic captures/samples in population
#'
#' Works in a manner similar to \code{mort()}, assigning a year to captured individuals and killing
#' them if sampling is fatal. Sex specific sampling is allowed.
#' 
#' EKJ modified for walrus to allow for sampling with replacement and 
#' self-recaptures (across years). Now returns a list with two elements:
#' indiv (as before) and samples (same structure as indiv but if an individual
#' is recaptured it will appear more than once)
#' EKJ modified to take a vector of ages for each sex
#'
#' @param indiv A matrix of individuals, as from makeFounders(), mate(), or mort().
#' @param n Number of captures (genetic samples)
#' @param year Capture year
#' @param fatal Is sampling fatal?
#' @param sex Sex specific sampling (either \code{"M"}, \code{"F"} or \code{NULL})
#' @param age Integer vector. The age class(es) at which sampling occurs
#' @param replace TRUE/FALSE whether animals should be sampled with replacement
#' @param self TRUE/FALSE whether or not self-recaptures are possible
#' @param maleages ages of males that can be sampled
#' @param femaleages ages of females that can be sampled
#' @export

walrusCapture <- function(indiv = makeFounders(), n = 1, year = "-1", fatal = TRUE, 
                          maleages = NULL, femaleages = NULL, #sex = NULL, age = NULL, 
                          replace = FALSE, self = FALSE) {
  
  if (!is.null(maleages)) {
    eligible.male <- indiv[,2] == "M" & indiv[,8] %in% maleages
  } else {eligible.male <- rep(FALSE, nrow(indiv))}
  
  if (!is.null(femaleages)) {
    eligible.female <- indiv[,2] == "F" & indiv[,8] %in% femaleages
  } else {eligible.female <- rep(FALSE, nrow(indiv))}
  
  # alive (and of the correct sex)
  # or considered dead
  is.alive <- is.na(indiv[,6]) & (eligible.female | eligible.male)
  is.dead  <- !is.alive
  
  # sample support
  n.alive <- sum(is.alive)
  
  # sample captures have to be alive
  # in case where recaptures are possible, actually want n here
  if (replace == FALSE) {n <- min(n, n.alive)} else {n <- n}
  
  if(n > 0) {
    
    # sample row location in indiv data.frame
    sample.loc <- sample.int(n.alive, size = n, replace = replace)
    
    # record capture
    indiv[is.alive,][sample.loc, 9] <- year
    
    if (self == TRUE){samples <- unique(indiv[is.alive,][sample.loc,])}
    
    # kill if capture sampling is fatal
    if (fatal) {
      indiv[is.alive,][sample.loc, 6] <- year
      indiv[is.alive,][sample.loc, 12] <- "LETHAL"; samples[,12] <- "LETHAL"
    } else {indiv[is.alive,][sample.loc, 12] <- "NONLETHAL"; samples[,12] <- "NONLETHAL"}
  }
  if (self == TRUE & exists("samples")){return(list(indiv = indiv, samples = samples))} else {
  return(list(indiv = indiv))}
}
