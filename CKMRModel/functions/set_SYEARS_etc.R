"set_SYEARS_etc" <- function( nlocal=sys.parent()) mlocal({
  # mlocal so that it creates vars
  
  PDYEARS <- FIRST_PDYEAR %upto% max( SYEARS)
  MAX_Dsyears <- length( SYEARS) # 0 up to last one
  
  DEVSTAGES <- cq( JuveF, AdF)
  Devstage_A <- autoloop( A=AGES,
      ifelse( AGES < AMAT, 'JuveF', 'AdF')
    )
  
  # PCOHORTS: Only need to consider pot-parents born within this
  # Earlier-born have same prob as first one, for all pot offs 
  # later-born can't be parents of any sample
  PCOHORTS <- (FIRST_PDYEAR - FEC_ASYMP_AGE) %upto%
      (max( SYEARS) - AMAT) 
  
})
