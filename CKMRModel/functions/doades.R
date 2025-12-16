"doades" <- structure( function( 
  m_SYEL,     # samp sizes
  Hbits,      # from get_Hbits
  reporteez,  # function that returns quants-of-interest
  lglk,       # prolly needed by reporteez, plus its env has "metadata"
  eps= 1e-5,  # for numvbderiv
  comp_wts= NULL
){
  ptru <- Hbits$PARS_FOR_H
  npars <- length( ptru)

  # Get param covar mat for this sampling scheme
  denv <- environment( lglk) # various useful things live here
  m_breakdown <- get_usable_SS( m_SYEL, denv$Devstage_A)
  
  # Other things required by make_n_comps:
  # IMO this is ugly code (why is it these 3 random-seeming quantities?!)
  # but getting it all generic and interdependent seems to 
  # decrease readability (at least in this case)
  ncomp_list <- make_n_comps( m_breakdown, 
      FIRST_PDYEAR= denv$FIRST_PDYEAR, 
      AMAT= denv$AMAT, 
      MAX_BGAP_HSPS= denv$MAX_BGAP_HSPS,
      MAX_WEANAGE = 5)
  
  H <- matrix( 0, npars, npars)
  
  # Loop over the different lglk components
  # (in this case MOP, HSP, SelfP) ... but I like writing generic code

  default_comp_wts <- rep( 1, length( Hbits$DSP))
  names( default_comp_wts) <- names( Hbits$DSP)
  
  if( !is.null( comp_wts)){ 
    # Allow user to specify just some of them, eg 
    # comp_wts= c( MOP=0)
    # Look for the right comp-type(s) inside name 
    # in the middle of the DSP name
    # dropping "s" for plural
    # I should know better than to write stuff like this...
    whicho <- match( sub( 's$', '', sprintf( '_%s', names( comp_wts))),
#        sub( '^[^_]*(_[^_s]+)_.*', '\\1', names( Hbits$DSP)), 0)
         sub( '^DSP(_.+)_.+', '\\1', names( Hbits$DSP)), 0)
    if( !all( whicho>0)){
stop( sprintf( "Don't understand which comp to weight: %s", 
        paste( names( comp_wts)[ whicho==0], collapse=', ')))}
    default_comp_wts[ whicho] <- comp_wts
  }
  comp_wts <- default_comp_wts
  
  for( comptype in names( Hbits$DSP)){
    this_ncomp <- ncomp_list[[ sub( 'DSP', 'n_comp', comptype)]]

stopifnot( my.all.equal( # dimensions had better match...
        unname( dimseq( Hbits$DSP[[ comptype]])[-1]), # 1st is params
        unname( dimseq( this_ncomp))
    ))
        
    this_H <- finfo_onetype( Hbits$DSP[[ comptype]], 
      this_ncomp)
    H <- H + this_H * comp_wts[ comptype]
  }
  # Other H-contributors, eg RE priors or general priors, should go here
  # H <- H + other_stuff
  
  Vpar <- solve( H)

  # Quants-of-interest: var via Delta-method
  # Just popdyn is quick (probs are slow) so don't parallel this
  stuff <- reporteez( ptru, lglk)
  
  Dstuff <- numvbderiv( reporteez, ptru, lglk=lglk, eps=1e-5)
  Vstuff <- Dstuff %**% Vpar %**% t( Dstuff)
  SEstuff <- sqrt( diag( Vstuff))

  # Up to the user to make sure reporteez() puts names on its results
  # but you don't absolutely HAVE to...
  names( SEstuff) <- names( stuff)
  dimnames( Vstuff) <- list( names( stuff), names( stuff))

  SEstuff@V <- Vstuff
  
return( SEstuff)
},
"doc" =
structure(c("doades    package:not-yet-a-package", "", "", "Design calculation for walrus", 
"", "", "DESCRIPTION", "", "This computes a set of SEs for user-specified quantities-of-interest, based on one specified sampling scheme, using pre-computed results ('Hbits') for one pop-dyn scenario. For Walrus.", 
"", "There's \"really\" three nested steps to CK Design: ", "", 
" - getting 'Hbits' (the magic derivatives), for one pop-dyn case;", 
" - turning that into a parameter covar mat, for one sampling-scheme;", 
" - getting \"interesting\" covar mat, for one particular set", 
" ", "'doades' combines the last two. This theoretically leads to slight inefficiency if computing different quants-of-int for the same pop-dyn and samp-scheme, coz the middle step gets done twice, but actually it is pretty quick.", 
"", "", "USAGE", "", "doades( m_SYEL, Hbits, reporteez, lglk, eps = 0.00001, comp_wts=NULL) ", 
"", "", "ARGUMENTS", "", " m_SYEL: 4D offarray of sample sizes by Sex, Year, E(stimated age), Lethality. These had better match the ranges already used in computing...", 
" ", " Hbits: ... which comes from 'despack::get_Hbits()', for one specific pop dyn scenario and (maximum) ranges of things like Age, Sampling-year, etc.", 
" ", " reporteez: Function taking \"true\" params and 'lglk' (see next), and returning a vector of \"interesting\" pop dyn summaries, eg abundance in some particular year-of-interest. The 'names()' of the result will be applied to the output of 'doades'. There's an example in the (non-runnable) EXAMPLES section.", 
" ", " lglk: You generated this long ago, from 'lglk_this_dataset <- add_data( lglk_walrus, <simfile>)'. Its environment contains useful constants, and you can run it with 'want=\"popdyn\"' to get just the pop-dyn summaries (which is presumably the first thing that your 'reporteez' function will do).", 
" ", " eps: for 'numvbderiv' numerical differentation. Not too big, not too small, is my advice. The anxious can try changing the default by a factor of 10 either way; it should make little difference.", 
" ", " comp_wts: How to weight the various lglk components, basically to explore the effect of (almost) turning off CKMR and/or IMR data. If non-NULL, it should be a numeric vector with names such as 'MOP\" or \"XmHSP\" or \"selfP\" (as seen inside names of 'Hbits$DSP'--- ie, the different types of pairwise comparison.  A value of '1' means \"treat this comp-type as-is\", and you only need to supply non-1 values, since 1 is the default. So, to explore eg turning off IMR, just set eg 'comp_wts=c( selfP=0.01)'--- which will leave the MOP and XmHSP weights at 1.", 
"", "", "VALUE", "", "A vector of standard errors, corresponding to the vector of things returned by 'reporteez()'. It will have an attribute 'V' containing the entire covariance matrix; usually just the SE's, ie 'ssqrt(diag(.))', are what's wanted, but You Never Do Know.", 
"", "", "SEE.ALSO", "", "'despack::get_Hbits'", "", "", "EXAMPLES", 
"", "## Don't run", "thingz_of_int <- function( p, L){", "  popd <- L( p)", 
"  Nfad_2015 <- popd$Nfad_y[ 2015]", "  Deaths_2025 <- popd$Nfad_y[ 2025] * ", 
"      (1-inv_logit( popd$lgt_fadsurv))", "  thingz <- unlist( returnList( # a handy way to name-and-get", 
"    Nfad_2015, Deaths_2025)", "return( thingz)", "}", "", "# ... now pass that into doades()", 
"", "## End don't run"), class = "docattr")
)
