"make_precision_table_completely" <- function( 
    stuff, Dstuff, Hbits, 
    data_env, 
    new_sslist
){
  ## Recalcs n_comps, then H, then Vpar, then Vstuff
  ## based on new_sslist
  
  extract.named( mget( 
    cq( AMAT, MAX_BGAP_HSPS, FIRST_PDYEAR), 
    data_env)
  )
  
  n_comps <- with( new_sslist, 
                   make_n_comps_1region(
                     m_YA,
                     m_FLive_YAN,
                     m_FLive_YAself,
                     m_FLive_YD,
                     FIRST_PDYEAR,
                     AMAT,
                     MAX_BGAP_HSPS
                   ))
  
  E_comp <- H_comp <- list()
  for( comptype in names( Hbits$DSP)){
    ncomp <- n_comps[[ sub( 'DSP', 'n_comp', comptype)]]
    H_comp[[ comptype]] <- finfo_onetype( Hbits$DSP[[ comptype]], 
                                          ncomp)
    E_comp[[ sub( 'DSP_', '', comptype)]] <- ncomp * 
      Hbits$Prkin[[ sub( 'DSP', 'Pr', comptype)]]
  }
  
  # Add up Hessian bits:
  H <- 0*H_comp[[1]]
  for( comptype in names( Hbits$DSP)){
    H <- H + H_comp[[ comptype]]
  }
  
  Vpar <- solve( H)
  Vstuff <- Dstuff0 %**% Vpar %**% t( Dstuff0)
  mpt <- make_precision_table( stuff0, sqrt( diag( Vstuff)), E_comp)
  attr( mpt, 'Vpar') <- Vpar
  return( mpt)
}