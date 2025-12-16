"effective_juve_surv" <- function( s1up){
## Takes juve surv vector over all ages that count as "juvenile"
## Finds the constant-at-age surv that would give same average
## survival for a randomly-chosen juve
## There might be other possible definitions, not sure...

  # First, numbers-at-age:
  # ... put 1 animal at first age, scale the rest;
  # surv of final juve age irrel for this part
  nju <- c( 1, cumprod( head( s1up, -1)))

  # Now an average of age-specific survs:
  avsurv <- nju %**% s1up / sum( nju)
return( avsurv)  
}
