"fec_check" <- function( EO){
## Either supply n_MOP_EYEYL directly, or extract it...
  if( !is.numeric( EO)){
    EO <- EO$n_MOP_EYEYL
  }
  
  AGES <- dimseq( EO, 1)
  AMIN <- AGES[1]
  AMAX <- tail( AGES, 1)
  SYEARS <- dimseq( EO, 2)
  LETHALITY <- dimseq( EO, 5)

  n_MOP_AL <- autoloop( 
      l=LETHALITY, ac_bj= AGES,
      SUMOVER= list( yc=SYEARS, aj=AGES, yj=SYEARS), {
    bj <- yj - aj;
    # ac_bj <- ac - (yc - bj); IE
    ac <- ac_bj + (yc-bj);
    (ac >= AMIN) * (ac <= AMAX) * # ensure in-range; 0 if not
      EO[ aj, yj, clippo( ac, AMIN, AMAX), yc, l]
  })
    
return( n_MOP_AL)
}
