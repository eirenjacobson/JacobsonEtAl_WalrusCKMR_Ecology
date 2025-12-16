"ideal_MOP_walrus" <- function( nlocal=sys.parent()) mlocal({
  r"---{
  Standard CK prob. If yc<bj, mother would have to survive post her sampling, and could have been immature then, so Pr_survival needs to take account of her age then.
  
  TRO is *not* adjusted to allow for ramp-up in fec around maturity. Small beer.
  
  PCOHORTS is potential-parent's birth-year. It's a wide range. We do not need the pop dyn to operate then, however; it is enough to know she was alive then.
  }---"
  
  Pr_MOP_byby_nonlethal <- autoloop( 
    bj= PDYEARS, yj = SYEARS,
    bc= PCOHORTS, yc= SYEARS, { 
      aj = yj - bj;
      (yc > bc) * # you can't be sampled before you exist!
        (yc - bc <= AMAX) * # coz then you die      
        (bc+AMAT <= bj) * # gotta be mature
        fec_a[ clippo( bj-bc, FECAGES)] * # cow's ERO (cf avg adult)
        recip_Nfad_y[ bj] *  # competition
          # ... how many more years must Mum live, starting from sampling?
          # NB if juv > 0 mum must have also survived its first year
        Pr_fsurv_ta[ pmax( 0, (bj+ifelse( aj>0, 1, 0))-yc), 
            clippo( yc-bc, SELF1AGES)]
    })
  
  Pr_MOP_byby_lethal <- autoloop( 
    bj= PDYEARS, yj = SYEARS,
    bc= PCOHORTS, yc= SYEARS, {
      aj = yj - bj;
      (yc > bc) * # you can't be sampled before you exist!
        (yc - bc <= AMAX) * # coz then you die      
        (bc+AMAT <= bj) * # gotta be mature
        fec_a[ clippo( bj-bc, FECAGES)] * # cow's ERO (cf avg adult)
        recip_Nfad_y[ bj] *   # competition
        (bj > yc)    # lethal sample must be strictly after birth year; simplified by RT
    })
  
  Pr_MOP_bybyl <- offarray( 0, dimseq=list( 
      bj= PDYEARS, yj= SYEARS, 
      bc= PCOHORTS, yc= SYEARS, lc= LETHALITY))
  Pr_MOP_bybyl[,,,,SLICE='LETHAL'] <- Pr_MOP_byby_lethal
  Pr_MOP_bybyl[,,,,SLICE='NONLETHAL'] <- Pr_MOP_byby_nonlethal
})
