"popdyn_walrus" <- function( nlocal=sys.parent()) mlocal({
  ## Population dynamics
  Nfad_y <- autoloop( Y=PDYEARS, {
    Nfad_ystart * exp( RoI * (Y-YSTART))
  })
  recip_Nfad_y <- 1/Nfad_y # I prefer calcing reciprocal once...
  
  # Juvenile version, assuming QEQM:
  # "Juves" are things included in Self comps (and probably as potential XHSPs and potential O's in MOPs, though those _could_ be something different I guess, at least if we had exact age)
  # 24/12: fixed minor (fortunately) bug that omitted one of 
  # the "+RoI"'s below
  Nfju_y <- Nfad_y * 
    ((1+RoI-fadsurv) / (1+RoI-fjusurv)) *
    ((fjusurv ^ -(AMAT-AMIN)) - 1)
  Nf_yd <- offarray( 0, dimseq=list( PDYEARS, DEVSTAGES))
  Nf_yd[,'JuveF'] <- Nfju_y
  Nf_yd[,'AdF'] <- Nfad_y
  
  ## Breeding psicle
  # col= FROM, row= TO (Mark's version) so that T %*% oldprob = newprob
  # Birth happens when Mum is in state 2
  Tbreedf <- rbind( 
    c( 0, psicle[2], psicle[3]),
    c( 1, 0, 0),
    c( 0, 1-psicle[2], 1-psicle[3])
  )
  
  ## Ppn ad fem calving in any given year 
  
  # ... neglects ramping-up as maturity is reached
  # No longer going the eigenroute; get stationary eivec by...
  # ... solving eigeneq with eigenval 1
  M <- Tbreedf - diag( 3)
  M[3,] <- 1
  Pbreedf <- solve( M, c( 0,0,1)) # eqm distro of states
  # Whole (adult) popln, at any one time
  ppn_breedy <- Pbreedf[2] # tho *really* should be age-wted somehow I guess
  recip_ppn_breedy <- 1 / ppn_breedy  
  
  # ... and Do Something for TMB version! 
  # IIRC there's an explicit soln given form of Tbreedf. And if not,
  # then there's other ways, but they're harder.
  
  # Powers of the breed trans mat
  Pr_Bt_B0T <- offarray( 0, 
                         dimseq=list( Bt=1:3, B0=1:3, T=0:MAX_popgap))
  Pr_Bt_B0T[,,SLICE=0] <- diag( 3)
  for( T in 1:MAX_popgap){
    Pr_Bt_B0T[,,T] <- Tbreedf %**% Pr_Bt_B0T[,,SLICE=T-1]
  }
  
  ## Fec at age: start 2(?) yrs prior to AMAT in might-get-preg
  # max( FECAGES) serves as "mini plus group" here--- assumed const after
  # Construction is a bit ugly, coz gotta get offset right
  
  fec_a <- offarray( 
    c( Pr_Bt_B0T[ SLICE=2,SLICE=3,FECAGES-(AMAT-2)]),
    dimseq=list( FECAGES))
  fec_a <- fec_a * recip_ppn_breedy
  fec_a[ max( FECAGES)] <- 1 # force asymptote!
  
  ## Distro of rebreed prob, given mum still alive
  # This is not QUITE the same as Pr_w, since we forcibly
  # wean animals at MAX_WEANAGE, but mum might still not give birth
  # [sic] "This only applies at MAX_WEANAGE"
  Pr_breedyagain_Db <- Pr_Bt_B0T[ SLICE=2, SLICE=2, ]

  ## ...e emoved a LOT of code for conditioning on nursingness

  ## Survival of adults; used in XHSP code
  Pr_fadsurv_t <- autoloop(
      t= 0:MAX_Tsep, # I *think* MAX_Tsep is limit...
    fadsurv ^ t
  )
  
  r"---{    
  Age-specific version (ie from given starting age, not necessarily adult). Only needed for starting ages up to AMAT; for animals known to be adult (eg parents), it's the same as Pr_fadsurv_t. Allows different survs for JUVEF and ADF. Currently used in ideal_selfP; should be used in ideal_MOP too.
  }---"
  
  # Create then loop is easiest here...  
  Pr_fsurv_ta <- offarray( 1, dimseq=list( t= 0:MAX_Dsyears, a= SELF1AGES))
  
  for( tt in 1 %upto% MAX_Dsyears){
    Pr_fsurv_ta[ tt,] <- Pr_fsurv_ta[ tt-1,] * 
        ifelse( SELF1AGES + tt-1 >= AMAT, fadsurv, fjusurv)
  }

  ## Ageing error setup already done 
  ## (assuming it does not depend on pop dyn--- OK in this case)  
  ## Movement stuff done separately
})
