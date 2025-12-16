#' breeding based on mature females
#' 
#' EKJ modified so that probability of breeding = 1 when in state 2
#' all mature males are considered as potential fathers (ie. not filtered by "stock")
#' tho note we may need to change this when we add USA/Russia
#' New individuals enter in Stock 0
#'
#' Returns an individual matrix with added newborns at age 0. This is the second
#' mate method, where the number of offspring is derived from the number of mature females,
#' such that each mature female produces a number of offspring specified by a sampling distribution,
#' and fathers are randomly drawn from all mature males within the mother's stock. Note that in
#' this mating system, *maturity by age* is specified, rather than *fecundity by age*. A single
#' probability distribution sets the number of offspring for each female, but the probability
#' that an individidual female is mature may vary by age. The same maturity by age structure applies
#' for males. It is possible, in cases where the maturity-by-age slope is shallow, that an individual
#' may be 'mature' in one breeding season, but then 'not mature' the next season.
#' Another mate method exists, where the number of newborns is set as a proportion of the population
#' size, and mating occurs until the required number of offspring are generated (if possible, given
#' breeding constraints) - see mate().
#'
#' @param indiv A matrix of individuals, e.g., generated in makeFounders() or output from move()
#' @param batchSize Numeric. The probability that a
#'                  female will have a (single) offspring, given that she is mature, if
#'                  fecundityDist = "binomial".
#'                  Used for all maturity structures. Note that, if run within a loop that goes
#'                  [move -> altMate -> mort -> birthdays], 'produced' is not the same as 'enters
#'                  age-class 1', as some individuals will die at age 0. Note for walrus, given that
#'                  females are being psi-cled around, batchSize will == 1 because these are females
#'                  that have entered the pregnancy state in the prev year
#' @param fecundityDist One of "poisson", "truncPoisson", or "binomial". Sets the distribution of the
#'                      number of offspring per mature female. Defaults to "poisson".
#' @param osr Numeric vector with length two, c(male, female), giving the sex ratio at birth
#'            (recruitment). Used to assign sexes to new offspring.
#' @param year Intended to be used in a simulation loop - this will be the iteration number, and
#'             holds the 'birthyear' value to give to new recruits.
#' @param firstBreedF Integer variable. The age at first breeding for females, default 1. The minimum age
#'                   at which individuals can breed. Applies to potential mothers and potential
#'                   fathers. 'firstBreed', 'maturityCurve', 'maleCurve, and 'femaleCurve' are
#'                   all capable of specifying an age at first breeding, and 'firstBreed' takes
#'                   precedence.
#' @param firstBreedM Integer. Age of first breeding for males                   
#' @param type The type of maturity-age relationship to simulate. Must be one of "flat",
#'             "age", or "ageSex". If "flat", the probability of parenthood is the same for
#'             all age:sex combinations above firstBreed. If "age", the probability that an
#'             individual is sexually mature is age-specific, set in 'maturityCurve'. If "ageSex",
#'             the probability that an individual is sexually mature is age- and sex-specific,
#'             set for males in 'maleCurve' and for females in 'femaleCurve'.
#' @param maxClutch Numeric value giving the maximum clutch / litter / batch / whatever size.
#'                  Reduces larger clutches to this size, for each breeding female.
#' @param singlePaternity TRUE/FALSE value indicating whether all the offspring produced by
#'                        female in a year should have the same father. Default TRUE. If
#'                        FALSE, each offspring will have a randomly-drawn father from within
#'                        the mother's stock. Note that this can lead to rapid exhaustion of
#'                        fathers if exhaustFathers = TRUE.
#' @param exhaustFathers TRUE/FALSE value indicating whether fathers should become 'exhausted'
#'                       by one breeding attempt. If exhausted, an individual will only mate with
#'                       one female, though may father more than one offspring - see
#'                       'singlePaternity' and 'batchSize'.
#' @param maturityCurve Numeric vector describing the age-specific maturity curve. One
#'                      value per age, over all ages from 0:max(indiv[,8]). Used if "type"
#'                      = "age". Note that 'firstBreed' can interfere with 'maturityCurve'
#'                      by setting maturities to zero for some age classes. Recommended
#'                      usage is to set 'firstBreed' to zero whenever 'maturityCurve' is
#'                      specified.
#' @param maleCurve Numeric vector describing age-specific maturity for males. One value per
#'                  age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                  Note that 'firstBreed' can interfere with 'maleCurve' by setting
#'                  maturities to zero for some age classes. Recommended usage is to set
#'                  'firstBreed' to zero whenever 'maleCurve' is specified.
#' @param femaleCurve Numeric vector describing age-specific maturity for females. One value
#'                    per age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                    Note that 'firstBreed' can interfere with 'femaleCurve' by setting
#'                    maturities to zero for some age classes. Recommended usage is to set
#'                    'firstBreed' to zero whenever 'femaleCurve' is specified.
#' @seealso [fishSim::mate()]
#' @export

walrusMate <- function(indiv = makeFounders(), batchSize = 0.5, fecundityDist = "poisson",
                    osr = c(0.5,0.5), year = "-1", firstBreedF = 7, firstBreedM = 15, type = "flat", maxClutch = Inf,
                    singlePaternity = TRUE, exhaustFathers = FALSE,
                    maturityCurve, maleCurve, femaleCurve) {
  if (!(type %in% c("flat", "age", "ageSex"))) {
    stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
  }
  if (!(fecundityDist %in% c("poisson", "truncPoisson", "binomial"))) {
    stop("'fecundityDist' must be one of 'poisson', 'truncPoisson', or 'binomial'.")
  }
  
  mothers <- subset(indiv, indiv[,2] == "F" & indiv[,8] >= firstBreedF & is.na(indiv[,6]) & indiv[,7] == 2)
  if(nrow(mothers) == 0) warning("There are no mature females in the population")
  fathers <- subset(indiv, indiv[,2] == "M" & indiv[,8] >= firstBreedM & is.na(indiv[,6]))
  if(nrow(fathers) == 0) warning("There are no mature males in the population")

  
  if (type == "ageSex") {
    
    #mothers <- mothers[runif(nrow(mothers)) < femaleCurve[mothers[,8]+1] , ,drop = FALSE]
    ## trims 'mothers' to just those that pass a random maturity test.
    
    #fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[,8]+1] , ,drop = FALSE]
    ## trims 'fathers' to just those that pass a random maturity test.
    ##        clutch <- rpois(n = nrow(mothers), lambda = batchSize)
   
    if(fecundityDist == "binomial") {clutch <- rbinom(nrow(mothers), 1, prob = batchSize)}
    
    mothers <- subset(mothers, clutch > 0)
    clutch <- clutch[clutch>0]
  }
  
  clutch[clutch > maxClutch] <- maxClutch ## delimits clutch sizes to not exceed maxClutch
  ## sprog.m <- matrix(data = NA, nrow = 0, ncol = 9) ## left empty if no-one breeds.
  #sprog.m <- walrusFounders(pop = 0) ## left empty if no-one breeds.
  
  for (s in 2) { ## s for 'stock'.
    mothersInStock <- mothers[mothers[,7] == s , , drop = FALSE]
    clutchInStock <- clutch[mothers[,7] == s]
    fathersInStock <- fathers
    if(nrow(fathersInStock) == 0) {
      warning (paste("There were no mature males in stock ",
                     s, ", so ", nrow(mothersInStock),
                     " mature females did not produce offspring",
                     sep = ""))
      ## sprog.stock <- matrix(data = NA, nrow = 0, ncol = 9)
      sprog.stock <- makeFounders(pop = 0)
    } else if(nrow(fathersInStock > 0)) {
      ## sprog.stock <- matrix(data = NA, nrow = sum(clutchInStock), ncol = 9)
      n.sprogs <- sum(clutchInStock)
      sprog.stock <- data.frame(Me = character(n.sprogs), Sex = character(n.sprogs),
                                Dad = character(n.sprogs), Mum = character(n.sprogs),
                                BirthY = integer(n.sprogs), DeathY = integer(n.sprogs),
                                Stock = integer(n.sprogs), AgeLast = integer(n.sprogs),
                                SampY = integer(n.sprogs), 
                                Region = character(n.sprogs), Independent = integer(n.sprogs),
                                Lethality = integer(n.sprogs))
      ticker <- 1
      for (m in 1:nrow(mothersInStock)) { ## m for 'mothers'
        if(nrow(fathersInStock) == 0) {
          warning(paste("All fathers in stock ", s, " are exhausted.", sep = ""))
        } else {
          ## assign mother, plus region of mother, and set as non-independent
          sprog.stock[ticker:(ticker+clutchInStock[m]-1), 4] <- mothersInStock[m,1]
          sprog.stock[ticker:(ticker+clutchInStock[m]-1), 10] <- mothersInStock[m,10]
          sprog.stock[ticker:(ticker+clutchInStock[m]-1), 10] <- 0
          if (singlePaternity == TRUE) {
            sprog.stock[ticker:(ticker+clutchInStock[m]-1), 3] <-
              fathersInStock[sample(1:nrow(fathersInStock), 1), 1]
          } else if (singlePaternity == FALSE) {
            if(nrow(fathersInStock) >= clutchInStock[m]) {
              sprog.stock[ticker:(ticker+clutchInStock[m]-1), 3] <-
                fathersInStock[sample(1:nrow(fathersInStock), clutchInStock[m]), 1]
            } else {
              sprog.stock[ticker:(ticker+nrow(fathersInStock)-1),3] <-
                fathersInStock[ ,1]
            }
          } ## Assign father(s).
          ## Note the potential conflict here with 'exhaustFathers' - but maybe not
          ## of concern, because there can't be many species with multiple paternity
          ## that also exhaust fathers after one mating attempt.
          if(exhaustFathers == TRUE) {
            fathersInStock <- fathersInStock[!fathersInStock[,1] %in% sprog.stock[,3],
                                             ,drop = FALSE]
            ## remove the used fathers from 'fathersInStock'
          }
          ticker <- ticker+clutchInStock[m]  ## increment ticker
        }
      }
    }
    sprog.stock <- sprog.stock[!is.na(sprog.stock[,3]), , drop = FALSE]
    
    sprog.stock[,1] <- uuid(n = nrow(sprog.stock), drop_hyphens = TRUE)
    sprog.stock[,2] <- sample(c("M", "F"), nrow(sprog.stock), TRUE, prob = osr)
    sprog.stock[,5] <- c(rep(year, nrow(sprog.stock)))
    sprog.stock[,6] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[,7] <- as.integer(rep(0), nrow(sprog.stock))
    sprog.stock[,8] <- c(rep(0, nrow(sprog.stock)))
    sprog.stock[,9] <- c(rep(NA, nrow(sprog.stock)))
    sprog.stock[,12] <- c(rep(NA, nrow(sprog.stock))) 
    ## sprog.stock <- dfify(sprog.stock)
    
    #sprog.m <- rbind(sprog.m, sprog.stock)
  }
  #names(sprog.m) <- names(indiv)
  indiv <- rbind(indiv, sprog.stock)
  
  return(indiv)
}