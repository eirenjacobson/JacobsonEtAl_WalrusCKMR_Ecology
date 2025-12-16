#### Compare *totals* of kin-pairs by type
# There are more checkables than just totals, but it's a start

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

for( i in seq_along( suffixes)){
  # Lglk fun for *this* dataset:
  # ... and a place to keep the data summary
  sim_lglks[[i]] <- add_data( lglk_walrus,
      simfile = sprintf( 'test/WalrusSamples_D0_L1_S0_%s.RData', suffixes[i]),
      YSTART = denv$YSTART) # presumably 2015
  denvi <- environment( sim_lglks[[ i]]) # where stuff lives

  compcheck[ 'MOP', 'Obs', i] <- sum(denvi$n_MOP_EYEYL)
  compcheck[ 'XmHSP', 'Obs', i] <- sum(denvi$n_XmHSP_EYEY)
  compcheck[ 'SelfP', 'Obs', i] <- sum(denvi$n_selfP_EYDY) 
} # for suffix (ie sim dataset)

save( sim_lglks, file='sim_lglks.RData')

# Expecteds:
compcheck[ 'MOP', 'Exp', ] <- sum(E$n_MOP_EYEYL) # same per sim
compcheck[ 'XmHSP', 'Exp', ] <- sum(E$n_XmHSP_EYEY)
compcheck[ 'SelfP', 'Exp', ] <- sum(E$n_selfP_EYDY) 

# All probs at once:
compcheck[ , 'P', ] <- ppois( compcheck[,'Obs',], compcheck[,'Exp',])

## 29/11: HSPs & Selfs look OK, but not MOPs
# Look at birth-gaps for XmHSPs:

birthgap_check( E)
birthgap_check( environment( sim_lglks[[1]])) # same as denv
# All sims (the t() means sims are rows, which compares better
bgsims <- t( do.on( sim_lglks, birthgap_check( environment( .))))
colMeans( bgsims)
birthgap_check( E)
# ... Those look pretty good to me. The smallest & largest @gap 3 are
# P~=0.03, 0.97. The largest @gap 2 is 0.9998... but it's waaay higher
# than the next biggest.



# MOPs: Look at distro by adult age-at-offspring-birth
sum( fec_check( E)) # matches compcheck
sum( fec_check( E)[1,]) # no Lethals expected
sum( fec_check( denv)[2,]) # ditto in simulatin
plot( 1:37, fec_check(E)[SLICE=2,], col='blue', ylim = c(0,50))
points( 1:37, fec_check(denv)[SLICE=2,], col='orange')
# NB _no_ sim parents at age 6..?
# Fec in E swings heavily over early years of adulthood.



#### Now check lglk derivs at true pars:
Dlglk0 <- matrix( 0, length( suffixes), length( ptru))

for( i in seq_along( suffixes)){
  Dlglk0[ i,] <- numvbderiv( sim_lglks[[i]], ptru, want='just_lglk')
}


# MVB version of loading the Nfad (and any other) results into a matrix
# Though Nfad & RoI (and anything else) shouldn't vary much, should they..?
sim_summ <- NULL
etemp <- new.env()
for (i in seq_along(suffixes)){
  load( sprintf( './simulation/test/Nfad_RoI_D0_L1_S0_%s.RData', suffixes[ i]), 
      envir=etemp)
  # I happen to know it's called 'out'
  
  if( i==1) {
    sim_summ <- matrix( 0, length( suffixes), length( etemp$out))
    dimnames( sim_summ) <- list( suffixes, names( etemp$out))
  }
  
  sim_summ[i,] <- etemp$out
  
  #  suffix <- suffixes[i]
  #  load(paste0("./simulation/Nfad_RoI_", suffix, ".RData"))
}
save( sim_summ, file='sim_summ.RData')  


