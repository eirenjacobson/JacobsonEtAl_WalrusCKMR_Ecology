"movement_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"--{
  I've only "tested" this for the homogenous case where everything moves around completely at random. It gets that one right, but..!
  
  Now added option when there's only one region, to avoid breaking later code. Fingers crossed..!
}--"
  
  if( n_REGIONS==1){
    # Will deffo be in same place..!
    Pr_rr_t <- offarray( 1, dimseq=list( 
      r1= REGIONS, r2= REGIONS, dt= 0:MAX_Tsep)
    )
    Pr_r <- 1
  } else { # 2 regions...
    # romega is the prob of Roaming from current state (are you Russianing or Americking?) to the other  
    # MVB's convention (more typical of general maths): 
    # Pr[new] = TM %**% Pr[ old]
    # IE row= to, col= from
    Tmove <- rbind( 
      c( 1-romega[1], romega[2]),
      c( romega[1], 1-romega[2])
    )
    # Equilib distro: No longer going the eigenroute; get stationary eivec by...
    M <- Tmove - diag( 2)
    M[2,] <- 1
    Pr_r <- solve( M, c( 0,1))
    
    r"---{
    Prob that *one* indiv is in r1 at start, and r2 at end, of tsep years
    Pr[r1,r2,dt] = Pr[r1] * Pr[r2|r1,t] = Pr[r2] * Pr[r1|r2,-t]
    Pr[r2|r1,t] is t-th power of transition matrix
    We need T^t %**% diag( Pr_region)--- RH bit there requires thought!
    No longer going the eigenway
    Direct loop easiest
    }---"  
    
    Pr_rr_t <- offarray( 0, dimseq=list( 
      r1= REGIONS, r2= REGIONS, dt= 0:MAX_Tsep)
    )
    
    Tmove_pow <- diag( 2)
    for( tsep in 0:MAX_Tsep){
      # NB next line could be efficientized in C--- RH thing is diag
      Pr_rr_t[ ,,tsep] <- Tmove_pow %**% diag( Pr_r)
      Tmove_pow <- Tmove %**% Tmove_pow
    }
  } # if 2 regions
  
  # Whether 1 or 2 regions
  names( Pr_r) <- REGIONS
  recip_Pr_r <- 1 / Pr_r
})