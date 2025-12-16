"add_data" <- function( 
  lglk= NULL,
  indiv= NULL,
  qk= NULL,
  simfile= "WalrusSamples.RData",
  simdir= "./simulation",
  AMAX= 37, # and then you die
  FIRST_PDYEAR= 2000,
  SYEARS, # REMOVE DEFAULT
  YSTART= 2015, # could use median PDYEARS but clarity is best
  AMAT= 6,
  AMIN= 1, # calves not modelled 
  MAX_WEANAGE= AMAT-1,
  MAX_BGAP_HSPS= 2*AMAT + 2, # GGPs unlikely within this
  FEC_ASYMP_AGE= 15,
  MAXAGERR= 0, # bold!
  PPN_2KP_FALSE_NEG= 0, # possibly unwise to have such a low default... 0.15?
  nonsparse= FALSE, # if TRUE, tweak sample sizes and kin for better independence
  ... # just to allow ending all previous arglines with commas
){
  # Auto-trimming is good (remove calves)
  op <- options( offarray_table_warn_trim=FALSE)
  on.exit( options( op))
  
  if( ...length()){
stop( sprintf( "Unknown args: %s", 
    paste( ...names(), collapse=', ')))
  }
  
  # Some "constants", independent of sampling
  SEXES <- cq( F, M)
  LETHALITY <- cq(LETHAL, NONLETHAL)
  
  AGES <- AMIN:AMAX
  # For self-recap: at 1st cap time, all adult ages are equiv
  SELF1AGES <- AMIN:AMAT # AMAT serves as "plus group" for selfs

  r"--{  
  Assume ageing error does not change this range! IE we only give "plausible" estimated ages,  not an age of minus-3 etc!  Also that ageing error does NOT mess up ID of "Devstage".

  Age "A" in simulations is assumed to ALREADY INCLUDE ageing error, even tho in Nov 2024 there isn't any yet. To add ageing error in future,  there should be another column "Atru" in the sim output.
  
  From now on in these data-summary steps, A means "measured A",  andbut I use "E" in the final arrays to emphasize that it's Estimated--- also for consistency with all prob calcs.  The ideal_<blah>_prob() functions all start with A (meant as real), and then transform the probs via age_blur_probs_walrus() to allow for error.
  }--"
  
  # Sim data: might be either be in simfile, or previously loaded
  # or might be non-existent--- in which case, will make all-zero
  # dataset waaaay below here
  if( is.null( indiv) && !is.null( simfile)){
    if( !file.exists( simfile)){
      simfile <- file.path( simdir, simfile)
    }
    
    what <- load( simfile) # indiv, or samples, or ....
    if( 'indiv' %not.in% what){
      # work out what its name was...
      for( i in what){
        indiv <- get( i)
        if( (indiv %is.a% 'data.frame') &&
            all( indiv %has.name% cq( Me, Mum, Dad, SampY)))
          break
      }
    } # now we should have indiv
    # last_raw_indiv <<- indiv # naughty code, but useful...
  } # loading of "raw" data
stopifnot( exists( 'indiv', environment(), inherits=FALSE))  

  # MVB fix: very important!!!
  indiv <- within( indiv, {
    Dad[ Dad=='founder'] <- NA
    Mum[ Mum=='founder'] <- NA
  })

  samples <- indiv[!is.na(indiv$SampY),]

  # SYEARS: make sure there's no gaps
  syrange <- if( missing( SYEARS)){
    range( samples$SampY)
    } else {
      range( SYEARS)
    }
  SYEARS <- syrange[1] %upto% syrange[2]
  
  # Confusingly, the next function does _not_ actually create
  # SYEARS: we've just done that. But it used to... 
  # ... and it still does create other things ("etc")...
  set_SYEARS_etc()

  if( is.null( qk)){
    recaps <- which( duplicated( samples$Me)) # EKJ NOTE:...
    # ... this won't happen with sim data in main sim file
    # bc simulation overwrites capture years. AllSamples has this info.
    
    nrecaps <- length( recaps)
    nuniq <- nrow( samples) - nrecaps
    if( nrecaps>0){
      # Re-order to put them at the end; simpler later
      samples <- samples[ c( seq_len( nuniq) %except% recaps, recaps),]
    }
    # Normal quickin() wants entire sim input, including unsampled. But 
    # can do POPs & sibs (not GGPs) just from samples, as below
    # qk <- quickin( samples[ seq_len( nuniq),], max_gen = 2, 
    #     cousins. = FALSE, thiatics. = FALSE)
    qk <- quickin_only_popsib( samples[ seq_len( nuniq),])

    if( nrecaps>0){
      # Add in POPs & HSPs for the recap too
      me_again <- samples$Me[ -seq_len( nuniq)]
      me_first <- match( me_again, samples$Me)

      qk$POP <- add_recap_kin( qk$POP, me_first, me_again)
      qk$MatHSP <- add_recap_kin( qk$MatHSP, me_first, me_again)
    }

    # last_qk <<- qk
  }

  # Extract key variables...
  extract.named( 
    rename.els( 
      samples[ cq( Me, Sex, SampY, BirthY, Lethality)],
      BirthY='B', SampY='Y', Lethality = 'L')
  ) # Me, Sex, Y, B, L now available

  # (Estimated) age at sampling:
  A <- Y - B

  # Find *indices* of matching pairs
  MOPs <- matrix( match( qk$POP, Me), nrow( qk$POP), 2)
  XmHSPs <- matrix( match( qk$MatHSP, Me), nrow( qk$MatHSP), 2)

  # Make sure in order. Even if DNAge imperfect, can prolly tell
  # which-is-which. I think later code assumes the birth-order is correct
  swappo <- B[ MOPs[,1]] > B[ MOPs[,2]]
  MOPs[ swappo,] <- MOPs[ swappo, 2:1]

  swappo <- B[ XmHSPs[,1]] > B[ XmHSPs[,2]]
  XmHSPs[ swappo,] <- XmHSPs[ swappo, 2:1]

  # Offspring considered: only non-calfs born post 2000
  offposs <- (A>0) & (B >= FIRST_PDYEAR)
  # Parents: zap calves, and samples taken too late in the study to mature
  # (latter is optional, they'll get all-zero probs anyway).
  # In principle, calves too could be used as potential future mothers
  # (ie A==0 is possible)
  # but they are not part of the model in general, so zap em  
  parposs <- (A>0) & (B <= (max(PDYEARS) - AMAT))

  # Elim founders... shouldn't be any!!!
  print( sum( offposs & is.na( samples$Dad))) # 0, whew
  offposs <- offposs & !is.na( samples$Dad) # in case

  # Filter MOPs & XmHSPs based only on pop-dyn year-range *for now*
  MOPs <- MOPs[parposs[MOPs[,1]] & offposs[MOPs[,2]],]
  XmHSPs <- XmHSPs[offposs[ XmHSPs[,1]] & offposs[ XmHSPs[,2]],]

  if( FALSE){ # Debugging code to demonstrate lack of 6yo new mums
    # Parent's age at offspring's birth
    Ac_Bj <- A[ MOPs[,1]] - (Y[ MOPs[,1]] - B[ MOPs[,2]])
    table( sort(Ac_Bj))
  }

  # Oh yeah--- the parent in a MOP had better be a Mother not a Father!!!
  MOPs <- MOPs[ samples$Sex[ MOPs[,1]]=='F',]

  ## SAMPLE SIZES: 
  # m_subset_xyz is samp size of subset-type samples with covars x,y,z
  # if subset omitted, it's all samples.
  # Now no attempt to pre-decide "poss par" and/or "poss off"
  # those will be incorped at the n_comp... stage
  # NB 'dimseq' arg in case ends-of-ranges don't have any samples
  # to make sure arrays are full-size

  # Using E for (Estimated) age hereon, instead of A(ge).

  m_SYEL <- offarray( table( 
      S=Sex, Y, A, L), 
      dimseq=list( SEXES, SYEARS, AGES, LETHALITY))

  ## Self-recaps: done from separate file
  # NB that SS is _slightly_ different in the separate file
  # cos recaps are excluded from the original. But it's <2% so use orig
  # Not quite perfect; should use absolutely everything, then adjust as for XmHSPs
  # Only first & last captures are used, but there aren't many threecaps
  # in fact just 3 animals in the "base" simulation
  # Look for file in same folder as simfile
  allsimfile <- file.path( dirname( simfile),
                           'All' %&% basename( simfile))
  what <- load( allsimfile)        
  if( 'allsamples' %not.in% what){
    # work out what its name was...
    for( i in what){
      allsamples <- get( i)
      if(( allsamples %is.a% 'data.frame') &&
         all( allsamples %has.name% cq( Me, Mum, Dad, SampY)))
        break
    }
  } # now we should have allsamples

  allsamples <- rename.els( 
    allsamples[ cq( Me, Sex, SampY, BirthY, Lethality)],
    BirthY='B', SampY='Y', Lethality = 'L')
  allsamples$A <- with( allsamples, Y-B)
  # Females only, no calves
  allsamples <- allsamples %where% ((Sex=='F') & (A>0))

  # For Xtuple recaps, just use 1st and last...
  # this is not 100 %ideal; should follow HSP approach as above
  # and recalc sample sizes purely for SelfPs. But anyway.
  # MVB: I have run out of meaningful variable names here
  metab <- table( allsamples$Me)
  threecaps <- names( metab[ metab>2])
  if( length( threecaps)){
    scrungio <- which( allsamples$Me %in% threecaps)
    blurk <- allsamples$Me[ scrungio]
    firsto <- match( threecaps, blurk)
    lasto <- length( blurk) + 1 - match( threecaps, rev( blurk))
    middlo <- scrungio[ -c( firsto, lasto)]
    if( length( middlo)){
      allsamples <- allsamples[ -middlo,]
    }
  }

  recaps <- names( metab[ metab==2])
  scrungio <- which( allsamples$Me %in% recaps)
  blurk <- allsamples$Me[ scrungio]
  firsto <- match( recaps, blurk)
  lasto <- length( blurk) + 1 - match( recaps, rev( blurk))
  selfPs <- cbind( scrungio[ firsto], scrungio[ lasto])
  
  ## Numbers of comparisons, and of actual kin-pairs
  # Create offarrays n_comps_<blah>
  if( !nonsparse){ # old 2024 version
    # To avoid over-changing existing code: this sparsity-assuming version is 
    # a non-nested function, whose args gotta be passed in
    # explicitly, and whose results are returned as a list--- which 
    # extract.named then turns into variables right here
    extract.named( 
      make_n_comps(
        mlist= get_usable_SS( m_SYEL, Devstage_A), # for different kin-types
        FIRST_PDYEAR,
        AMAT,
        MAX_BGAP_HSPS,
        MAX_WEANAGE # added 16/1/2025: now needed coz allowing *some* MOP comps when yj==yc
      ))
  } else {
     # This modifies eg m_YE to account for nonsparsity, as well as creating 
     # n_comps_<blah>. It's a nested function
     # so doesn't need args or return-value--- and it's easier to follow!
     nested_make_n_comps_nonsparsely()
  }
  
  ## MOPs:
  n_MOP_EYEYL <- offarray(table(
    aj = A[MOPs[,2]],
    yj = Y[MOPs[,2]],
    ac = A[MOPs[,1]],
    yc = Y[MOPs[,1]],
    lc = L[MOPs[,1]]
  ), template = n_comp_MOP_EYEYL) # use these dims 
  
  # Zero-out unusable cases, even if they exist
  zap_MOP <- whichoff( n_comp_MOP_EYEYL==0 & n_MOP_EYEYL>0)
  n_MOP_EYEYL[ MATSUB = zap_MOP] <- 0

  ## HSPs... make it compatible with offposs etc
  r"--{
     First we apply false-neg filter. This is done *after* creating the full XmHSP list and applying any nonsparsity stuff, which is a bit weird coz you couldn't do it with real data. The reason is to avoid logical problems with triads where some pairwise outcomes are falsely lost. This is painful detail. Cross the real-data bridge when we come to it.
  }--"
  XmHSPs <- XmHSPs[ runif( nrow( XmHSPs)) > PPN_2KP_FALSE_NEG,]  
  n_XmHSP_EYEY <- offarray( table(
    a1= A[ XmHSPs[,1]],
    y1= Y[ XmHSPs[,1]],
    a2= A[ XmHSPs[,2]],
    y2= Y[ XmHSPs[,2]]
  ), template = n_comp_XmHSP_EYEY) # use these dims
  
  # Zero-out unusable cases, even if they exist
  zap_HSP <- whichoff( n_comp_XmHSP_EYEY==0 & n_XmHSP_EYEY>0)
  n_XmHSP_EYEY[ MATSUB = zap_HSP] <- 0
  
  ## selfPs: comps first.
  # Indices are different here
  
  # number of self-recaps. 
  # MUST use A etc from allsamples NOT general A! 
  # MVB 8/11: "NOOFF" must now be called "nooff", for some reason
  n_selfP_EYDY <- with( allsamples,  # CRUCIAL!
      offarray(table(
        aself1 = pmin( A[selfPs[,1]], AMAT),
        y1 = Y[selfPs[,1]],
        d2 = Devstage_A[ A[ selfPs[,2]], nooff=TRUE], # nooff stops warnings
        y2 = Y[selfPs[,2]]
      ), template = n_comp_selfP_EYDY # use these dims
  )) 
  
  # Zero-out unusable cases, even if they exist
  zap_selfP <- whichoff( n_comp_selfP_EYDY==0 & n_selfP_EYDY>0)
  n_selfP_EYDY[ MATSUB = zap_selfP] <- 0
  
  # Maybe tell the user that we've squished some data...
  nzaps <- do.on( cq( MOP, HSP, selfP), 
      nrow( get( sprintf( 'zap_%s', .))))
  if(  sum( nzaps)> 0){
    nzaps <- nzaps %except% 0
warning( "Pairs found but ignored cozza zero usable comps: " %&%
    paste( sprintf( '%s: %i', names( nzaps), nzaps), collapse=', '))
  }
  
  ## Ageing error "data"--- assumed exact for now,
  # but code allows generality
  expand_by_maxagerr <- function( ARANGE)
    max( min( ARANGE)-MAXAGERR, 1) %upto% (max( ARANGE) + MAXAGERR)
  
  # This gets a matrix of the right shape, but with hardwired no-error
  Pr_agerr_a <- offarray( 0, 
      first=c( -(1+MAXAGERR), 1),
      last=c( 1+MAXAGERR, AMAX))
  Pr_agerr_a[ 0,] <- 1 # no agerr for now, regardless of MAXAGERR

  setup_ageing_errmats() # for much later

  
  ## Return stuff
  env <- if( is.null( lglk)) .GlobalEnv else environment( lglk)
  env <- new.env( parent=env)
  e <- list2env( envir=env, returnList( 
    n_comp_MOP_EYEYL,
    n_MOP_EYEYL,
    n_comp_XmHSP_EYEY,
    n_XmHSP_EYEY,
    n_comp_selfP_EYDY,
    n_selfP_EYDY,
    
    m_SYEL, # not the breakdowns from prepare_usable_SS()--- though why not???
    
    zap_MOP,
    zap_HSP,
    zap_selfP,
    
    Pr_e_Y_sampled, # samp age compo (including error...)
    Pr_a_ey,
    self_Pr_a_ey, # agged adults (sort-of...)
    
    SEXES,
    LETHALITY,
    FIRST_PDYEAR,
    PDYEARS,
    MAX_Db= diff( range( PDYEARS)), # birth gap HSPs
    YSTART,
    
    SYEARS,
    MAX_Dsyears= diff( range( SYEARS)), # recap gap selfPs
    
    MAXAGERR,
    Pr_agerr_a,
    AMIN,
    AMAX,
    AMAT,
    FEC_ASYMP_AGE,
    FECAGES= (AMAT-1):FEC_ASYMP_AGE,
    MAX_WEANAGE,
    MAX_BGAP_HSPS,
    MAX_popgap = 50, # who knows..?
    AGES,
    SELF1AGES,
    DEVSTAGES,
    Devstage_A,
    
    PCOHORTS,
    FIRST_PCOHORT= PCOHORTS[1],
    LAST_PCOHORT= tail( PCOHORTS, 1),
    
    MAX_Tsep= 100, # excessive, but too hard to precalc for now..!
    
    PPN_2KP_FALSE_NEG, # needed by lglk_walrus
    
    nsamples = nrow(samples),
    nallsamples = nrow(allsamples)
  ))
  
  if( !is.null( lglk)){
    copy_lglk <- lglk # actually unnecessary, but clearer...
    environment( copy_lglk) <- e
    return( copy_lglk)
  } else {
    return( e)
  }
}
