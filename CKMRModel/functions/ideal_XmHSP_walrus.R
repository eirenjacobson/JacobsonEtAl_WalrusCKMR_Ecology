"ideal_XmHSP_walrus" <- function( nlocal=sys.parent()) mlocal({
## Non-spatial only here.
  r"--{ 
  Spatial version has the extra complexity that juve will stay with mum until weaning (by defn) so the two sibs are not moving independently. Code is quite a lot more complicated.
  
  This HSP code incorporates ppn_breedy explicitly, rather than computing a TRO
  }--"
  
  Pr_XmHSP_bb <- autoloop(
      b1= PDYEARS, b2= PDYEARS,
    (b1 < b2) *   # non-symmn
      Pr_fadsurv_t[ abs( b2-b1)] * # avoid OOB
      Pr_breedyagain_Db[ abs( b2-b1)] * 
      recip_Nfad_y[ pmax( b1, b2)] *
      recip_ppn_breedy *
      (1-PPN_2KP_FALSE_NEG) # new in 2025
    )
})
