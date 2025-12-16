"prepare_H_E" <- function( denv, Hbits){
  E <- list() # by pairwise comp type, eg MOP

  # H contains expected Hessian (d^2 lglk / dparams^2)
  # ie inverse covariance between parameters
  # Add it up from all lglk components
  
  # MOPs
  # Not bothering with Nursy status now (no subscript N)
  ncomp <- denv$n_comp_MOP_EYEYL 
  H <- finfo_onetype( Hbits$DSP$DSP_MOP_EYEYL, ncomp)
  E$n_MOP_EYEYL <- ncomp * Hbits$Prkin$Pr_MOP_EYEYL


  # Half sibs
  ncomp <- denv$n_comp_XmHSP_EYEY
  H <- H + finfo_onetype( Hbits$DSP$DSP_XmHSP_EYEY, ncomp)
  E$n_XmHSP_EYEY <- ncomp * Hbits$Prkin$Pr_XmHSP_EYEY

  # self
  ncomp <- denv$n_comp_selfP_EYDY
  H <- H + finfo_onetype( Hbits$DSP$DSP_selfP_EYDY, ncomp)
  E$n_selfP_EYDY <- ncomp * Hbits$Prkin$Pr_selfP_EYDY

return( list( H=H, E=E))
}
