"make_precision_table" <- function( stuff, SE, Ekin){
  is_CVable <- grepl( '^N', names( stuff)) # abunds only (?!)
  
  EMOP_df <- as.data.frame( Ekin$MOP_EYEYL)
  EMOP_df$bj <- with( EMOP_df, yj-aj)
  
  EHSP_df <- as.data.frame( Ekin$XmHSP_EYEY)
  EHSP_df$b2 <- with( EHSP_df, y2-a2)
  
  ESelfP_df <- as.data.frame( Ekin$selfP_EYDY)
  ESelfP_df$b1 <- with( ESelfP_df, y1-a1)
  
  SE_or_CV <- SE
  SE_or_CV[ is_CVable] <- 100 * SE[ is_CVable] / stuff[ is_CVable]
  stuff <- c( stuff, 
              'EMOP25+'= with( EMOP_df, sum( response[ bj >= 2025])),
              'EHSP25+'= with( EHSP_df, sum( response[ b2 >= 2025])),
              'ESelfP25+'= with( ESelfP_df, sum( response[ y1 >= 2025])
              )
  )
  
  stuff[ stuff > 100] <- round( stuff[ stuff> 100])
  rbind( signif( stuff, 3), signif( c( SE_or_CV, NA, NA, NA), 2))
}
