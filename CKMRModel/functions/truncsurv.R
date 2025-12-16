"truncsurv" <- function( AMAX=30,  phi){
## If adults survive AMAX years beyond maturity, with a surv of phi
## up til then, at which point it's 0: what would the equiv surv of
## a population that had no cut-off death rate?

  # truncsurv ( 1, 0.9) # check: yep, should be ~0.5
  # truncsurv( 30, 0.9622) # 0.9458
  N <- rep(0, AMAX)
  N <- phi ^ (0:AMAX)
  surv <- phi + 0*N
  surv[ length( surv)] <- 0
  N %**% surv / sum(N)
}

"truncsurvsenesce" <- function(phi1, phi2, ASEN = 24, AMAX = 30){
  ## If adults survive AMAX years beyond maturity, with a surv of phi1
  ## from AMAT until ASEN-1 and phi2 from ASEN to AMAX
  ## at which point it's 0: what would the equiv surv of
  ## a population that had no cut-off death rate?
  
#  AMAX <- 30 # yrs post maturity
#  ASEN <- 24 # yrs post maturity
#  phi1 <- 0.99
#  phi2 <- 0.55
  N <- rep(0, AMAX+1)
  N[1:(ASEN+1)] <- phi1 ^ (0:ASEN)
  N[(ASEN+2) : (AMAX+1)] <- phi1^ASEN * phi2^(1:(AMAX-ASEN))
  surv <- c(rep(phi1, ASEN), rep(phi2, (AMAX - ASEN + 1)))
  surv[ length( surv)] <- 0
  N %**% surv / sum(N)
}