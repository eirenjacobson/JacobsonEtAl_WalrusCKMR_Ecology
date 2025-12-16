"make_n_comps" <- function( 
  mlist, # should contain m_YE etc--- see below
  FIRST_PDYEAR,
  AMAT,
  MAX_BGAP_HSPS,
  MAX_WEANAGE
){
  r"--{
  mlist presumbaly from calling prepare_usable_SS, or hand-tweaking afterwards
  
  Sometimes we will not use results of certain kin-comps. 
  EG POP when caught in same year (to avoid unweaned offspring being sampled near mother); 
  HSP when birth-gap is long enuf that GGP becomes likely. 
  These were _not_ necessarily filtered out when MOPs and HSPs were made earlier.
  }--"
  
  # Create these things (and ONLY these things):
  extract.named( mlist[ cq( m_YE, mF_YEL, mliveF_YE,  mF_YD)])
  
  SYEARS <- dimseq( m_YE, 1) # slightly more efficient than dimseq()[[1]]
  AGES <- dimseq( m_YE, 2)
  LETHALITY <- dimseq(mF_YEL, 3) 
  
  SELF1AGES <- dimseq( mliveF_YE, 2)
  DEVSTAGES <- dimseq( mF_YD, 2)
  
  # _Could_ do consistency checks on all input dimensions
  # but leave that to R (it's really up to the caller!!!)
  
  # All numbers-of-pairwise-comparisons take the form:
  # possible-zeroing * number-first * number-second 

  # Lethal-potoff OK: for potpar, condition on Lethality
  n_comp_MOP_EYEYL <- autoloop(
    aj= AGES, yj = SYEARS,
    ac = AGES, yc = SYEARS, 
    lc = LETHALITY, {
      bad <- yc - ac
      bju <- yj - aj
      
      (bju >= FIRST_PDYEAR) *
        (bju >= (bad+AMAT)) *
        ((yj != yc) | (aj > MAX_WEANAGE)) * # no risk offspring & mum nearby
        m_YE[yj, aj] * mF_YEL[yc, ac, lc]
    })

  # HSP: lethal OK
  n_comp_XmHSP_EYEY <- autoloop(
    a1= AGES, y1 = SYEARS,
    a2 = AGES, y2 = SYEARS, {
      b1 <- y1 - a1
      b2 <- y2 - a2
      
      (b1 >= FIRST_PDYEAR) * (b2 >= FIRST_PDYEAR) *
        (b2>b1) * # don't double-count
        ((b2-b1) <= MAX_BGAP_HSPS) * # GGP risk beyond this
        m_YE[y1,a1] * m_YE[y2,a2]
    })
  
  # Self-recap: 
  # female only, 1st nonlethal, 2nd optional
  # Assumes ageing error does NOT change assessment of "Devstage"
  n_comp_selfP_EYDY <- autoloop(
    a1= SELF1AGES, y1= SYEARS,
    d2= DEVSTAGES, y2= SYEARS,{
      (y2 > y1) * # duhhh
        mliveF_YE[ y1, a1] * mF_YD[ y2, d2]
    })
  
  # Return all n_comp_<blah>
  returnList( mget( ls( environment(), pattern='^n_comp_')))
}
