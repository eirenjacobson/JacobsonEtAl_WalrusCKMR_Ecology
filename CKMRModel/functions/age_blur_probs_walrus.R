"age_blur_probs_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{
  "A" is true age: unknown, but that's what affects probs
  "E" stands for Estimated age (now forced to have same range as actual ages).
  
  Gonna special-case MAXAGERR==0 coz code should be much faster (and that's the base case anyway).
  
  If MAXAGERR>0, this should be re-written to loop over AGERR=-MAXAGERR:+MAXAGERR only, since that will much quicker than running over the full age range; for most measured ages, most _real_ ages are not possible (you can't really be 20 if DNAage says 30).
  }---"
  
  #     abs( B[ HSPs[,1]] - B[ HSPs[,2]]) <= MAX_BGAP_HSPS,] 
  
  eyHSP <- list( 
    e1= AGES, y1= SYEARS,
    e2= AGES, y2= SYEARS)
  
  if( MAXAGERR==0){ # e1==a1, e2==a2
    Pr_XmHSP_EYEY <- autoloop( 
        indices= eyHSP, {
      b1 <- y1-e1
      b1clip <- pmax( b1, FIRST_PDYEAR)
      b2 <- y2-e2
      b2clip <- pmax( b2, FIRST_PDYEAR)

      (b1 >= b1clip) * (b2 >= b2clip) *
        Pr_XmHSP_bb[b1clip, b2clip] 
    })
  } else {
    Pr_XmHSP_EYEY <- autoloop( 
        indices= eyHSP, 
        SUMOVER= list( a1= AGES, a2= AGES), {
      b1 <- y1-a1
      b1clip <- pmax( b1, FIRST_PDYEAR) # Anti OOB
      b2 <- y2-a2
      b2clip <- pmax( b2, FIRST_PDYEAR)

      (b1 >= b1clip) * (b2 >= b2clip) *
        Pr_XmHSP_bb[ b1clip, b2clip] * 
        Pr_a_ey[ a1, e1, y1] * Pr_a_ey[ a2, e2, y2]
    })
  }
  
  ## MOPs: 
  r"---{
  NB for efficiency, birth-year of Parent is pushed up to min( PCOHORTS), since any parent born prior to that will have reached asymp fec by the time the offspring is born
  }---"

  eyeylMOP <- list(
    ej= AGES, yj= SYEARS,
    ec= AGES, yc= SYEARS, lc=LETHALITY) 
  
  if( MAXAGERR==0){ # ej==aj, ec==ac
    Pr_MOP_EYEYL <- autoloop(
        indices= eyeylMOP, {
      bj <- yj-ej
      bjclip <- pmax( bj, FIRST_PDYEAR)

      bc <- yc-ec
      bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)

      (bj==bjclip) * (bc <= LAST_PCOHORT) *
        Pr_MOP_bybyl[bjclip, yj, bcclip, yc, lc]
    })
  } else { # AGE ERROR
    Pr_MOP_EYEYL <- autoloop(
        indices= eyeylMOP,
        SUMOVER= list( aj= AGES, ac= AGES), {
      bj <- yj-aj
      bjclip <- pmax( bj, FIRST_PDYEAR)

      bc <- yc-ac
      bcclip <- clippo( bc, FIRST_PCOHORT, LAST_PCOHORT)

      (bj==bjclip) * (bc <= LAST_PCOHORT) *
        Pr_MOP_bybyl[bjclip, yj, bcclip, yc, lc] *
        Pr_a_ey[ aj, ej, yj] * Pr_a_ey[ ac, ec, yc]
    })
  }
  
  ## Self recap: assumes no Devstage error in 2nd animal
  # abs( y2-y1) below since prob may as well be computed for impossibles---
  # but gotta avoid OOB---
  # easier than deciding which are possible!
  # unless you specifically ask it to return several things.
  # Semicolon is optional, but arguably clearer for separating
  # statements within autoloop()

  eySelf <- list( 
    e1= SELF1AGES, y1= SYEARS,
    d2= DEVSTAGES, y2= SYEARS)

  if( MAXAGERR==0){ # ej==aj, ec==ac
    Pr_selfP_EYDY <- Pr_selfP_AYDY
  } else {
    Pr_selfP_EYDY <- autoloop( 
        indices=eySelf, 
        SUMOVER= list( a1= SELF1AGES, a2= SELF1AGES), {
      d2 <- ifelse( a2 >= AMAT, "AdF", "JuveF");
      Pr_selfP_AYDY[ pmin( a1, AMAT), y1, d2, y2] * 
        Pr_a_ey[ a1, e1, y1]
    })
  }  
})
