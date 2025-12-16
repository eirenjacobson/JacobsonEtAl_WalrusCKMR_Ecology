#' make a founding population
#' 
#' EKJ edit for walrus to allow Year 0 to be set 
#'
#' Returns a pop-by-8 data.frame, with each row being an individual in the
#' founder population.
#' [,1] is a unique (uuid) identifier for each animal.
#' [,2] is sex, values "M" or "F".
#' [,3] is "founder" for all animals.
#' [,4] is "founder" for all animals.
#' [,5] is the animal's birth year. Implicitly assumes that 'makeFounders' occurs at the
#'      very start of year 1, just after the 'birthdays' step of year 0.
#' [,6] is NA for all animals. Holds the death year in most cases.
#' [,7] is the stock membership for each animal. FOr walrus, stock means breeding phase.
#' [,8] is the age of each animal (in 'breeding seasons') at the beginning of year 1,
#'      given that birthdays occur at the very end.
#' [,9] is NA for all animals (will hold sample year)
#' [,10] is designed to hold region, but may not be used
#' [,11] is independence of animals (1 for founders bc they don't have moms)
#' [,12] is lethality of sampling (NA for now)
#' makeFounders() will throw a warning if osr, stocks, or survCurv do not sum to 1. It is
#' not strictly necessary that they sum to 1 (proportionality within each class is sufficient),
#' but error-checking and readability is easiest if they do sum to 1.
#'
#' @param pop The size of the founder population.
#' @param osr A numeric vector describing the sex ratio, c([male], [female]).
#' @param stocks A numeric vector describing the probability that an individual
#'               is in each stock.
#' @param maxAge Numeric. The max age to which an animal may survive.
#' @param survCurv Numeric vector. Describes the probability within the founder cohort of belonging
#'                 to each age-class for age=-classes 1:maxAge. Cannot be blank.
#' @param year1 default to 1, but can be set to a real year (e.g., 2000)
#' @param senescence default to FALSE
#' @seealso [fishSim::make_archive()]
#' @export
#' 
library(ids)

walrusFounders <- function(pop = 1000, osr = c(0.5,0.5), stocks = c(0.3,0.3,0.4),
                         maxAge = 20, survCurv = 0.7^(1:maxAge)/sum(0.7^(1:maxAge)),
                         year1 = 1, senescence = FALSE, regions = c(0.5, 0.5),
                         independenceAge = 6) {
  
  if(sum(osr) != 1) warning("osr does not sum to 1")
#  if(sum(stocks) != 1) warning("stocks do not sum to 1")
#  if(sum(survCurv) != 1) warning("survCurv does not sum to 1")
  if(length(survCurv) != maxAge) warning("survCurv and maxAge imply different maximum ages")
  
  ### indiv <- matrix(data = NA ,nrow = pop, ncol = 9)
  indiv <- data.frame(Me = character(pop), Sex = character(pop), Dad = character(pop),
                      Mum = character(pop), BirthY = integer(pop), DeathY = integer(pop),
                      Stock = integer(pop), AgeLast = integer(pop), SampY = integer(pop),
                      Region = character(pop), Independent = integer(pop), Lethality = integer(pop))
  indiv[,1] <- uuid(n = pop, drop_hyphens=TRUE)  ## uuid IDs for each animal.
  # assign sexes by probability
  indiv[,2] <- sample(c("M", "F"), pop, TRUE, prob = osr)
  indiv[,3] <- c(rep("founder", pop)) ## has no father
  indiv[,4] <- c(rep("founder", pop)) ## has no mother
  indiv[,6] <- as.integer(c(rep(NA, pop))) ## founders are not yet dead
  # assign ages and birth years
  indiv[,8] <- sample.int(maxAge, pop, TRUE, prob = survCurv)
  indiv[,5] <- year1 - indiv[,8] ## back-infer birth year from age
  indiv[,9] <- as.integer(c(rep(NA, pop))) ## founders are not yet sampled
  # assign psicle state
  indiv[,10] <- sample(c("USA","Russia"), pop, TRUE, prob = regions)
  indiv[,11] <- 1 # all founders are independent bc they don't have moms
  indiv[,12] <- as.integer(c(rep(NA, pop))) # lethality of sampling 
  indiv[,7] <- as.integer(sample(1:length(stocks), pop, TRUE, prob = stocks))
  # put all males and juveniles into Stock 0
  indiv[which(indiv$Sex == "M"),]$Stock <- 0
  indiv[which(indiv$AgeLast<4),]$Stock <- 0
  # put all 4-yr-old F into Stock 3
  indiv[which(indiv$Sex == "F" & indiv$AgeLast == 4),]$Stock <- 3
  # five-year-olds can only be in stock 1 or 3
  fiveyearoldf <- which(indiv$Sex == "F" & indiv$AgeLast == 5)
  indiv$Stock[fiveyearoldf] <- sample(x = c(1,3), size = length(fiveyearoldf), 
                                      prob = c(0.5, 0.5), replace = TRUE)
    # if senescence == TRUE, put all 30-plus-yr old F into Stock 0
  if(senescence == TRUE){
    indiv[which(indiv$Sex == "F" & indiv$AgeLast >= 30),]$Stock <- 0
  }
  return(indiv)
}