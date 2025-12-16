"interesting_stuff" <- function( pars, lglk, ...){
  raw_stuff <- lglk( pars, want='popdyn', ...)
  
  # Only *scalar* interesting things are allowed
  # Vectors are not interesting. End of.
  
  # Make it easy to change the list of calculees, and 
  # ... avoid having to explicitly refer to raw-stuff everywhere
  e <- new.env( parent=list2env( raw_stuff)) 
  evalq( envir=e, {
    # Well, these are not *that* interesting, really...  

# MVB 8/11: offarray now keeps identity even for scalar extractions
# so either need to wrap the output in c(), or use [SLICE=blah]
# otherwise you get "incomparables" problem at the rel_ad_15_25 stage
    Nfad_2015 <- Nfad_y[ SLICE=2015] 
    Nfad_2020 <- Nfad_y[SLICE = 2020]
    Nfad_2025 <- Nfad_y[ SLICE=2025] 
    rel_ad_15_25 <- Nfad_y[ SLICE=2025]/Nfad_y[ SLICE=2015] - 1
    extract.named( real_params[ cq( fadsurv, fjusurv)])
    ppn_breedy <- ppn_breedy
  })
  
  # Must return a vector of interesting scalars!
  stuff <- unlist( as.list( e))
  
  # Sorting by names helps a bit... abunds first
  iscap <- substring( names( stuff), 1, 1) %in% LETTERS
  o <- order( !iscap, names( stuff))
  stuff <- stuff[ o]
return( stuff)  
}
