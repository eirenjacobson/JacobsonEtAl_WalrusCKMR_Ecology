"make_n_comps_1region" <-
function( 
    m_YA,
    m_FLive_YAN,
    m_FLive_YAself,
    m_FLive_YD,
    FIRST_PDYEAR,
    AMAT,
    MAX_BGAP_HSPS
){
  ## Abstracted from add_data, for re-use in pure designs
  SYEARS <- dimseq( m_YA)[[1]]
  AGES <- dimseq( m_YA)[[2]]
  NURSINGNESS <- dimseq( m_FLive_YAN)[[3]]
  
  SELF1AGES <- dimseq( m_FLive_YAself)[[2]]
  DEVSTAGES <- dimseq( m_FLive_YD)[[2]]
  
  # _Could_ do consistency checks on all input dimensions
  # but leave that to R (it's really up to the caller!!!)
  
  # Half-way house between region-aware and not...  
  # cos none of the arrays are regionized!
  # So, keep "region" dimensions but make 'em all the same
  REGIONS <-'Arctic' 
  sole_samp_region <- REGIONS[1]
  
  # Cheekily, I am really using A (exact age) but calling it E
  n_comp_MOP_REYREYN <- autoloop(
    rj = REGIONS, aj= AGES, yj = SYEARS,
    rc = REGIONS, ac = AGES, yc = SYEARS, nc=NURSINGNESS, {
      bad <- yc - ac
      bju <- yj - aj
      
      (rj==sole_samp_region) * (rc==sole_samp_region) *
        (bju >= FIRST_PDYEAR) *
        (bju >= (bad+AMAT)) *
        (yj != yc) * # should only apply to pre-max-weanage, but...
        m_YA[yj, aj] * m_FLive_YAN[yc, ac, nc] 
    })
  n_comp_MOP_REYREYN <- sumover( n_comp_MOP_REYN, 'nc')
  
  # At least for now, nursy means not spatial
  n_comp_MOP_EYEYN <- sumover( n_comp_MOP_REYN, cq( rj, rc))
  
  n_comp_XmHSP_REYREY <- autoloop(
    r1 = REGIONS, a1= AGES, y1 = SYEARS,
    r2 = REGIONS, a2 = AGES, y2 = SYEARS, {
      b1 <- y1 - a1
      b2 <- y2 - a2
      
      (r1==sole_samp_region) * (r2==sole_samp_region) *
        (b1 >= FIRST_PDYEAR) * (b2 >= FIRST_PDYEAR) *
        (b2>b1) * # don't double-count
        ((b2-b1) <= MAX_BGAP_HSPS) * # GGP risk beyond this
        m_YA[y1,a1] * m_YA[y2,a2]
    })
  
  n_comp_selfP_RYARDY <- autoloop(
    r1= REGIONS, y1= SYEARS, a1= SELF1AGES,
    r2= REGIONS, d2= DEVSTAGES, y2= SYEARS,{
      (y2 > y1) * # duhhh
        m_FLive_YAself[ y1, a1] * m_FLive_YD[ y2, d2]
    })
  
  # Return all n_comp_<blah>
  returnList( mget( ls( environment(), pattern='^n_comp_')))
}
