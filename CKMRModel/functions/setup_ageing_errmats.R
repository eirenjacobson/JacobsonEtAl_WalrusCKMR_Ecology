"setup_ageing_errmats" <- function( nlocal=sys.parent()) mlocal({
## Callable only from add_data()
  r"--{
  Assumes PRE-KNOWN matrix of ageing "error" (uncertainty):
  Pr_agerr_a
  
  MAXAGERR is greatest poss +/- deviation between truth & actual;
  limiting the range of error is much more efficient than using a full-sized AMAX*AMAX matrix
  
  A is *true* age here; E is estimated
  }--"

  # Ageing error in fact assumed to be zero for 2024 results,
  # but code is all ready for it (except prelim deconvolution).
  
  # Distro of age in samples (meant for ageing-error stuff)
  Pr_e_Y_sampled <- offarray( table( A, Y),
      dimseq=list( A=AGES, Y=SYEARS)
  ) # not normed yet
  recip_denom <- 1/sumover( Pr_e_Y_sampled, 'A')
  Pr_e_Y_sampled <- autoloop(
      A= AGES, Y= SYEARS,
    Pr_e_Y_sampled[ A, Y] * recip_denom[ Y]
  ) # normed to sum to 1 per year
  
  # Now we *should* deconvolve this, to get a priori est of
  # Pr[ TRUE-age | Y & was-sampled]
  # If no ageing error, there's nothing to do ;)
  if( MAXAGERR>0){
warning( "Skipping deconvolution--- forgiveable for Design..?")
  }
  Pr_Atru_Y_sampled <- Pr_e_Y_sampled

  agerr_indices <- list( a= AGES, e= AGES, y= SYEARS) 
  
  # this is the numerator of Bayes' theorem
  Pr_a_ey <- autoloop( 
    indices= agerr_indices,
    (abs( a-e) <= MAXAGERR) *
      Pr_agerr_a[ clippo( e-a, -MAXAGERR, MAXAGERR), a] *
      Pr_Atru_Y_sampled[ a, y]
  )
  
  # ... must sum to 1
  cumpr <- sumover( Pr_a_ey, 'a')
  # this is the denominator of Bayes' theorem
  recip_cumpr <- 1/cumpr # minimize divisions by doing once only (faster)
  Pr_a_ey <- autoloop(
    indices= agerr_indices,
    Pr_a_ey[a,e,y] * recip_cumpr[ e,y]
  )
  
  # Tedium to deal with ageing error in first sample for SelfP
  # All adult ages (Amat up) are condensed to a single class
  self_Pr_a_ey <- Pr_a_ey[ SELF1AGES, SELF1AGES, ]
  
  # Easy when e<AMAT:
  self_Pr_a_ey[ SLICE=AMAT,,] <- 
      sumover( Pr_a_ey[ AMAT:AMAX,SELF1AGES,], 1)
  
  # When e=AMAT, gotta do weighted sum...
  glurp <- autoloop( 
      A=AGES, Y=SYEARS,
      SUMOVER=list( E=AMAT:AMAX),
    Pr_a_ey[A,E,Y] * Pr_e_Y_sampled[ E, Y]
  )
  
  self_Pr_a_ey[ , SLICE=AMAT,] <- glurp[ SELF1AGES,] # OK for juves but...
  self_Pr_a_ey[ SLICE=AMAT, SLICE=AMAT,] <- 
      sumover( glurp[ AMAT:AMAX,], 1)
      
  r"--{ 
  That's technically valid, but slightly inconsistent with the (invalid) assumption that Devstage is always right even in the presence of ageing error, which is assumed for the 2nd sample. Getting the latter right is too horrible to contemplate; will wait for full age-structured model
  }--"
})
