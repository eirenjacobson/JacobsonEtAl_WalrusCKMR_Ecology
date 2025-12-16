"nested_make_n_comps_nonsparsely" <- function( nlocal=sys.parent()) mlocal({
  r"--{
  Sometimes we will not use results of certain kin-comps. 
  EG POP when caught in same year (to avoid unweaned offspring being sampled near mother); 
  HSP when birth-gap is long enuf that GGP becomes likely. 
  These were _not_ necessarily filtered out when MOPs and HSPs were made earlier.
  
  ALSO this version of make_n_comps definitely applies the nonsparsity patches.
  
  In theory, it should first remove all prior-to-last-capture SelfP members from all MOP and XmHSP comps, and adjust sample sizes for MOPs & XmHSPs accordingly. However, that's already been done in Eiren's two-part dataset structure. If rewriting from scratch, I would do things differently, starting from one file that includes ALL samples, but that's more effort and waaaay too much work for now.
  
  Next, it removes all O's in MOPs from XmHSP comparisons.
  
  Next, it counts triads etc in SelfPs and XmHSPs, removes the superfluous pairs, and adjusts overall sample sizes (somewhat crudely) for those two types of comparison to bring the expected total into line with the superfluity-reduced outcomes.
  }--"

  get_usable_SS( m_SYEL, Devstage_A) |> extract.named()
  # ... generate m_YE, mF_YEL, mliveF_YE,  mF_YD from m_SYEL
  
  # All numbers-of-pairwise-comparisons take the form:
  # possible-zeroing * number-first * number-second 

  ## MOPs
  # Lethal-potoff OK: for potpar, condition on Lethality
  n_comp_MOP_EYEYL <- autoloop(
    aj= AGES, yj = SYEARS,
    ac = AGES, yc = SYEARS, 
    lc = LETHALITY, {
      bad <- yc - ac
      bju <- yj - aj
      
      (bju >= FIRST_PDYEAR) *
        (bju >= (bad+AMAT)) *
        ((yj != yc) | (aj > MAX_WEANAGE)) * # no risk offspring & mum nearby
        m_YE[yj, aj] * mF_YEL[yc, ac, lc]
    })

  ## HSP
  # Discard known O's
  O2drop <- unique( MOPs[,2])
  
  # Adjust sample size first
  mYE2drop <- offarray( table( Y[ O2drop], A[ O2drop])) # mebbe >1 in some YE combos
  droppies <- whichoff( mYE2drop>0)
  m_YE[ MATSUB=droppies] <- m_YE[ MATSUB=droppies] - mYE2drop[ MATSUB=droppies]
  
  # Then drop any actual XmHSPs involving those O's
  m1 <- match( XmHSPs[,1], O2drop, 0)
  m2 <- match( XmHSPs[,2], O2drop, 0)
  XmHSPs <- XmHSPs[ m1==0 & m2==0, ]
  
  # Now we adjust for sib triads etc.
  orig_nrow_XmHSPs <- nrow( XmHSPs)  
  XmHSPs <- detriad_HSPs( XmHSPs, B) # fairly complicated code

  # Adjust sample sizes for HSP comps across-the-board, to match the reduced
  # total of (now independent) HSPs. 
  adjuss <- nrow( XmHSPs) / orig_nrow_XmHSPs
  m_YE <- m_YE * sqrt( adjuss) 

  # Finally, set #comps based on adjusted SS. NB Lethal is OK for HSP
  n_comp_XmHSP_EYEY <- autoloop(
    a1= AGES, y1 = SYEARS,
    a2 = AGES, y2 = SYEARS, {
      b1 <- y1 - a1
      b2 <- y2 - a2
      
      (b1 >= FIRST_PDYEAR) * (b2 >= FIRST_PDYEAR) *
        (b2>b1) * # don't double-count
        ((b2-b1) <= MAX_BGAP_HSPS) * # GGP risk beyond this
        m_YE[y1,a1] * m_YE[y2,a2]
    })
  # And... breathe.
  
  
  ## Self-recap: 
  # female only, 1st nonlethal, 2nd optional
  # Assumes ageing error does NOT change assessment of "Devstage"
  n_comp_selfP_EYDY <- autoloop(
    a1= SELF1AGES, y1= SYEARS,
    d2= DEVSTAGES, y2= SYEARS,{
      (y2 > y1) * # duhhh
        mliveF_YE[ y1, a1] * mF_YD[ y2, d2]
    })
  
  # Nothing to return: this is a nested function that 
  # modifies/creates things directly in its caller
})
