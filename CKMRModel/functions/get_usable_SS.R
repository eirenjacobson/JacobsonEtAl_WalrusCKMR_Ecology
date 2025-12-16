"get_usable_SS" <- function( m_SYEL, Devstage_A){
  r"--{
  Extract usable SS for different types of comp, eg female-only, live, seg-agged, etc. Starts from full array of samps by [Sex, Year, E(stimated age), Lethality]. Returns a list suitable for passing to make_n_comps()
  
  Even this isn't the full story, coz make_n_comps() will omit certain comparisons (ie will zero the number of comps in certain categories) on the basis that those results would not be useful (eg mother & child sampled in same year).

  For most design calcs this is enough, but for a few things (eg effects of not using CK) you might want to tweak some SSes before make_n_comps (or afterwards).
  }--"
  
  ddd <- dimseq( m_SYEL)
  SEXES <- ddd[[1]]
  SYEARS <- ddd[[2]]
  AGES <- ddd[[3]]
  LETHALITY <- ddd[[4]]
  
  # Devstage_A is the most convenient way to summarize age biol
  all_ages <- dimseq( Devstage_A)
  AMAT <- all_ages[ match( 'AdF', Devstage_A)]
  AMAX <- tail( all_ages, 1)
  DEVSTAGES <- unique( Devstage_A) # JuveF comes first... ?calves?
  # SELF1AGES: all (used) juve ages, plus first adult age (for aggregating)
  SELF1AGES <- all_ages[1] %upto% AMAT
  
  mF_YEL <- m_SYEL[ SLICE='F',,,]  # potential Mothers
  m_YE <- sumover(m_SYEL, cq( S, L)) # pot Offspring/Sib: either sex, dead OK

  ## SelfP sample sizes
  # Only Females (hence F at start)
  # First samp must be non-lethal

  # Gotta partially aggregate by age: 
  # 1st samp: all adults -> AMAT  
  # 2nd samp: split by Devstage
  mselfF_YEL <- mF_YEL[,SELF1AGES,] # correct dim...
  # ... put all adults into AMAT
  mselfF_YEL[,AMAT,] <- sumover( mF_YEL[,AMAT:AMAX,], 'A')
  
  # First sample for self: must be nonlethal and female
  mliveF_YE <- mselfF_YEL[,, SLICE='NONLETHAL']

  # Second sample for self: lethality status is irrelevant
  mF_YD <- autoloop(
      Y=SYEARS, D=DEVSTAGES, 
      SUMOVER=list( A=SELF1AGES, L=LETHALITY),
    mselfF_YEL[ Y, A, L] * (Devstage_A[ A]==D)
  )
  
returnList( 
    mF_YEL, # mums
    m_YE,   # offspring, sibs
    mliveF_YE, # self 1st time
    mF_YD   # self 2nd time 
  )
}
