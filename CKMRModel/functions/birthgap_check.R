"birthgap_check" <- function( EO, max_db=9){
## Either supply n_XmHSP_EYEY directly, or extract it...
  if( !is.numeric( EO)){
    EO <- EO$n_XmHSP_EYEY
  }
  
  AGES <- dimseq( EO, 1)
  AMIN <- AGES[1]
  AMAX <- tail( AGES, 1)
  SYEARS <- dimseq( EO, 2)
 
  # Birth gaps: anything above max_db goes into max_db
  DBRANGE <- 1:max_db

  # this is a bit inefficient but never mind
  n_HSP_Db <- autoloop( 
      db= DBRANGE,
      SUMOVER= list( a1=AGES, y1=SYEARS, a2=AGES, y2=SYEARS), {
    actual_db <- clippo( abs( (y1-a1) - (y2-a2)), DBRANGE);
    (actual_db == db) *
      EO[ a1, y1, a2, y2]
  })
    
return( n_HSP_Db)
}
