#### Compare *totals* of kin-pairs by type
# There are more checkables than just totals, but it's a start

# 22/12: replicates are in the test folder
compcheck <- array( 0, 
    c( 3, 3, length( suffixes)),
    dimnames=list( 
      cq( MOP, XmHSP, SelfP),
      cq( Exp, Obs, P),
      suffixes
    )
  )

# Loop over sim datasets, loading...
sim_lglks <- list() # we will keep this for later

tempenv <- new.env()
outs <- matrix( 0, 2, length( suffixes), 
    dimnames= list( names( out), suffixes))
ptrus <- matrix( ptru, length( ptru), length( suffixes),
    dimnames= list( names( ptru), suffixes))
    
for( i in seq_along( suffixes)){
  # Lglk fun for *this* dataset:
  # ... and a place to keep the data summary 
  # Give them names
  suff <- suffixes[ i]
  sim_lglks[[ suff]] <- add_data( lglk_walrus,
      simfile = sprintf( 'test/WalrusSamples_%s_%s.RData', TEST_CASE, suff),
      YSTART = denv$YSTART) # presumably 2015
  denvi <- environment( sim_lglks[[ i]]) # where stuff lives

  compcheck[ 'MOP', 'Obs', i] <- sum(denvi$n_MOP_EYEYL)
  compcheck[ 'XmHSP', 'Obs', i] <- sum(denvi$n_XmHSP_EYEY)
  compcheck[ 'SelfP', 'Obs', i] <- sum(denvi$n_selfP_EYDY) 

  # Expected values: 
  # First, get "right" params for *this* simulation
  load( envir=tempenv, sprintf( './simulation/test/Nfad_RoI_%s_%s.RData', 
    TEST_CASE, suff))
  outs[,i] <- tempenv$out

  ptrui <- ptru
  ptrui[ 'log_Nfad_y0'] <- with( tempenv, 
    log( unname( out[ 'Nfad_2000'] * exp( 15*out['RoI']))))
  ptrui[ 'RoI'] <- tempenv$out[ 'RoI']
  ptrus[,i] <- ptrui
  
  probsi <- sim_lglks[[ suff]]( ptrui, want='probs')
  
  # Could use E for expecteds (calced in fit_walrus_CKMR_xxxx.r), 
  # but that's based on values of "out" (ie empirical param "ests" from 
  # complete info) for one specific simulation
  # Instead, Expecteds via ncomps * prob for this p'tic ptru
  # tho in fact all ptrui are *fairly* similar (~1% diff?)
  # as you would hope
  
  for( ikin in dimnames( compcheck)[[1]]){
    # Get full name of relevant array (incl subscripts)
    # Case-insensitive coz I was not consistent...
    wotizit <- ls( denvi, patt= sprintf( '(?i)n_comp_%s', ikin)) |>
        xsub( 'n_comp_', '')
    n_comp <- denvi[[ sprintf( 'n_comp_%s', wotizit)]]
    pr <- probsi[[ sprintf( 'Pr_%s', wotizit)]]
    Ekin <- n_comp * pr
    compcheck[ ikin, 'Exp', i] <- sum( Ekin)
  }
} # for suffix (ie sim dataset)

save( sim_lglks, file='sim_lglks.RData')

if( FALSE){ 
  # Expecteds straight from E
  compcheck[ 'MOP', 'Exp', ] <- sum(E$n_MOP_EYEYL) # same per sim
  compcheck[ 'XmHSP', 'Exp', ] <- sum(E$n_XmHSP_EYEY)
  compcheck[ 'SelfP', 'Exp', ] <- sum(E$n_selfP_EYDY) 
}

# All probs at once:
compcheck[ , 'P', ] <- ppois( compcheck[,'Obs',], compcheck[,'Exp',])

# Fisher's method for combining multiple independent P-values
pchisq( -2 * rowSums( log( compcheck[,'P',])), 2*length( suffixes))

# MOPs fine (on total, at least...), not so HSP or SelfP
# where there's about a 6% and 4% bias respectively in the totals
rowMeans( compcheck[,'Exp',] / compcheck[,'Obs',])

r"--{
  Is there a sim-specific pattern? 
  When using same Nfad_RoI value for each, perhaps yes... each seems either high or low, for all 3 comp types. 
  Not really obvious with new sim-specific Nfad_RoI
}--"
matplot( t( compcheck[,'Obs',] / compcheck[ ,'Exp',]), 
    ylim=c( 0.5, 1.5), xlim=c( 0, 10),type='p')
abline( h=1)

## 29/11: HSPs & Selfs look OK, but not MOPs
# Look at birth-gaps for XmHSPs:

birthgap_check( environment( sim_lglks[[1]])) # same as denv
# All sims (the t() means sims are rows, which compares better
bgsims <- t( do.on( sim_lglks, birthgap_check( environment( .))))
colMeans( bgsims)
birthgap_check( E)

# Fairest is to renorm the observed ones so the totals match:
colMeans( bgsims) * sum( birthgap_check( E)) / sum( colMeans( bgsims))
# ... very close

# MOPs: Look at distro by adult age-at-offspring-birth
# fec_check() takes nearly a whole second, so store it :)
fcE <- fec_check( E)
sum( fcE) # matches compcheck
sum( fcE[ SLICE='LETHAL',]) # no Lethals expected
fcO <- fec_check( denv)
sum( fcO) # also matches compcheck
fcEnl <- fcE[ SLICE='NONLETHAL',]
fcOnl <- fcO[ SLICE='NONLETHAL',]
ylimmo <- 1.05 * max( c( fcEnl, fcOnl))
plot( 1:37, fcEnl, ylim=c( 0, ylimmo), col='blue',
  xlab= "Mother's age at offspring's birth",
  main= "Reprod output by age class")
points( 1:37, fcOnl, col='orange')
# NB _no_ sim parents at age 6..?
# Fec in E swings heavily over early years of adulthood.
# 


#### Now check bias, based on lglk derivs at true pars
Dlglk0 <- esterr <- 0 * ptrus

# Parallelize for speed. Same code as "fit_walrus_ckmr_<xxxx>.r"
ncores <- parallel::detectCores( logical=FALSE)-2 # leave spare...
ncores <- pmin( ncores, 8) # for now...
CLUSTO <- makeCluster( ncores)
# For some reason, that gives problems 1st time (windows security)
# and then is OK 2nd time..?
registerDoParallel( CLUSTO, ncores, nocompile=TRUE)

# Subsidiary functions needed by lglk_walrus:
export_funs <- rownames( mvbutils::foodweb( 
    where=.GlobalEnv, prune='lglk_walrus', 
    plotting=FALSE, ancestors=FALSE, descendents=TRUE)$funmat
  )

for( i in seq_along( suffixes)){
  # For parallel derivs, bundle all necessary functions in the envir
  # To avoid "environmental pollution" (cozza reference semantics)
  # make a deep copy of all objects already in that envir
  
  # Arguably neater than same-effect code in 'fit_walrus_ckmr_xxxx'
  
  sim_lglki <- sim_lglks[[ i]]
  e <- environment( sim_lglki) # for brevity
  denvi <- list2env( mget( lsall( e), envir=e), parent=.GlobalEnv)
  list2env( mget( export_funs), envir=denvi) # ... and add functions
  environment( sim_lglki) <- denvi

  Dlglk0[ ,i] <- numvbderiv_parallel( 
      sim_lglki, ptrus[,i], want='just_lglk',
      eps=1e-5, 
      FOREACH_ARGS=list( 
        .packages= c( 'mvbutils', 'offarray')
      )
    )
  delta <- Vpar %**% Dlglk0[,i] # coz Vpar = -solve( H)
  # Point est would be close to ptru +/- delta; correct sign is "obvious" but
  # hard to guess!

  if( FALSE){ # code to check...
    L0 <- sim_lglks[[i]]( ptrus[,i], want='just_lglk')
    Lplus <- sim_lglks[[i]]( ptrus[,i] + delta, want='just_lglk')
    Lminus <- sim_lglks[[i]]( ptrus[,i] - delta, want='just_lglk')
    # ptru+delta is better than ptru-delta (and better than ptru)!
    # So, est would be at ptru+delta, but it should be at ptru if perfect
    # Thus the esterr calc below...
  }
  
  esterr[,i] <- (-delta)
} # bias

stopImplicitCluster()
stopCluster( CLUSTO)
rm( CLUSTO)

rowMeans( esterr) # absolute bias (on scale of parameters)
rowMeans( esterr) / sqrt( diag( Vpar)) # bias rel to SD
rowMeans( Dlglk0) / sqrt( diag( H)) # ditto

r"--{
We can only check overall bias on the *parameter* scale using the complete model ie using all kin-types, coz only then is the Hessian guaranteed pos-deft and therefore invertible. (Actually it might be pos-deft if selfP are left out, but in general this point is true.)

However, we can check for "asymptotic consistency" using Dlglk0 even just on individual kin-type components, because then H itself (not its inverse) is the relevant covariance matrix, so it doesn't need to be inverted. But H would need to be re-calculated using only the included kin-types. Alternatively, and perhaps more simply, if we have 100 realizations of Dlglk0, then we can easily use the empirical variance (per component) instead of re-calculating H.

100 is just a number I made up. But 10 is probably too few for this.
}--"

