"adjust_future_sample_sel" <- function( data_env, upju=1, from=2024){
  ## Change all samp size arrays so that juves are sampled
  # upju times as hard, but annual samp sizes kept the same
  
  # Drag them in here...
  extract.named( as.list( denv)[ cq(
    m_YA,
    m_FLive_YAN,    
    m_FLive_YAself,
    m_FLive_YD,
    AMAT
  )])
  
  # This seq of ops is repetitive, but a bit hard to functionize
  # So I shall be repetitive
  mtot <- sumover( m_YA, 'A')
  m_YA <- autoloop(
    indices= dimseq( m_YA),
    m_YA[ Y, A] * 
      ifelse( (A<AMAT) & (Y >= 2024), upju, 1)
  )
  mtot_new <- sumover( m_YA, 'A')
  m_YA <- autoloop(
    indices= dimseq( m_YA),
    m_YA[Y,A] * 
      mtot[ Y] / mtot_new[ Y]
  )
  m_YA[ VECSUB=which( is.na( m_YA))] <- 0
  
  mtot <- sumover( m_FLive_YAN, 'A')
  m_FLive_YAN <- autoloop(
    indices= dimseq( m_FLive_YAN),
    m_FLive_YAN[ Y, A, N] * 
      ifelse( (A<AMAT) & (Y >= 2024), upju, 1)
  )
  mtot_new <- sumover( m_FLive_YAN, 'A')
  m_FLive_YAN <- autoloop(
    indices= dimseq( m_FLive_YAN),
    m_FLive_YAN[Y,A,N] * 
      mtot[ Y,N] / mtot_new[ Y,N]
  )
  m_FLive_YAN[ VECSUB=which( is.na( m_FLive_YAN))] <- 0
  
  
  mtot <- sumover( m_FLive_YAself, 'A')
  m_FLive_YAself <- autoloop(
    indices= dimseq( m_FLive_YAself),
    m_FLive_YAself[ Y, A] * 
      ifelse( (A<AMAT) & (Y >= 2024), upju, 1)
  )
  mtot_new <- sumover( m_FLive_YAself, 'A')
  m_FLive_YAself <- autoloop(
    indices= dimseq( m_FLive_YAself),
    m_FLive_YAself[Y,A] * 
      mtot[ Y] / mtot_new[ Y]
  )
  m_FLive_YAself[ VECSUB=which( is.na( m_FLive_YAself))] <- 0
  
  mtot <- sumover( m_FLive_YD, 'D')
  m_FLive_YD <- autoloop(
    indices= dimseq( m_FLive_YD),
    m_FLive_YD[ Y, D] * 
      ifelse( (D=='JuveF') & (Y >= 2024), upju, 1)
  )
  mtot_new <- sumover( m_FLive_YD, 'D')
  m_FLive_YD <- autoloop(
    indices= dimseq( m_FLive_YD),
    m_FLive_YD[Y,D] * 
      mtot[ Y] / mtot_new[ Y]
  )
  m_FLive_YD[ VECSUB=which( is.na( m_FLive_YD))] <- 0
  
  returnList( 
    m_YA,
    m_FLive_YAN,  
    m_FLive_YAself,
    m_FLive_YD
  )
}