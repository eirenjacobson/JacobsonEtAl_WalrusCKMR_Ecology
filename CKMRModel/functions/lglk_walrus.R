"lglk_walrus" <- function( 
  params, 
  want= c( 'just_lglk', 'lglk', 'probs', 'popdyn', 'ALL')
){
## All the steps in lglk calculation, based on data/constants already in this
## function's environment and put there by add_data(). Can be run _without_ 
## data (just ranges) to just calc probs or just popdyn, based on params

r"--{
Default behaviour is to return lots of things, as a list, equivalent to 'want="ALL"'. If you just want a _scalar_ lglk value, set 'want="just_lglk". Otherwise, you can specify a subset of the options in 'want', in which case relevant variables will be returned as a list. For example, 'want=c("probs","popdyn")' will return numerous popdyn quantities (such as numbers-at-age-and-year) as well as all the _observable_ pairwise kinship prob arrays (ie dependent upon observable attributes of animals, such as _estimated_ age, but not _true_ age).

Those options are useful in different settings: "just_lglk" for fitting to data (or RTMB) where a scalar result is required; "probs" for 'get_Hbits()' (see how it is called in the main script); "popdyn" for 'interesting_stuff()'. If you don't need the probabilities or lglk value, then 'want="popdyn"' is fastest.
}--"
  
  if( MAXAGERR>0){
stop( "Ageing error not fully implemented yet")
    # That's only _slightly_ true now; AFAIK all the code is ready, you
    # just gotta supply Pr_agerr_a[,] and MAXAGERR
  }
  
  ## First needs to have its envir set, to know about DATA
  # Pretty lazy partial test!
stopifnot( exists( 'PDYEARS', environment( sys.function()), 
      inherits=FALSE))
  
  ## Unpack parameters
  # get em one (or more) at a time
  # mitigates counting errors!
  iparam <- 0
  next_param <- function( n=1) {
    iparam <<- iparam + n
    unname( params[ seq( to=iparam, length=n)])
  }
  
  Nfad_ystart <- exp( next_param())
  RoI <- next_param()
  
  lgt_fadsurv <- next_param()
  diff_lgt_fjusurv <- next_param() # 0 means same surv as adult
  fadsurv <- inv.logit( lgt_fadsurv)
  fjusurv <- inv.logit( lgt_fadsurv + diff_lgt_fjusurv)
  
  # breeding cycle
  psicle <- offarray( inv.logit( next_param( 2)),
      first=2, last=3)
  
  assign( 'last_params', params, environment( sys.function())) # debugging etc
  
  # Do what? Might return early...  
  want <- match.arg( want, several.ok=TRUE)
  if( 'ALL' %in% want){
    # everything in the default argument
    want <- eval( formals( sys.function())$want)
  }
  returnables <- list() # filled in below
  
  ## General setup (boring housekeeping) if any
  # Well, there isn't any; all done in popdyn, or 
  # previously in add_data
  # so don't do it:
  # setup_walrus()

  ## Pop dyn
  popdyn_walrus()
  if( 'popdyn' %in% want){
    returnables <- returnList(
      Nfad_y,
      Nfju_y,
      Pbreedf,
      ppn_breedy,
      recip_ppn_breedy,
      fec_a,
      
      real_params= unlist( returnList( 
        Nfad_ystart,
        RoI,
        fadsurv,
        fjusurv,
        psicle
      ))
    )
  }
  want <- want %except% 'popdyn'
  if( !length( want)){
return( returnables) # that's all folks
  }
  
 # movement_walrus() # even if only 1 region!
  
  ## Ideal probs (if we knew everything)
  # although these DO incorporate space
  ideal_MOP_walrus()
  ideal_XmHSP_walrus() # and ignore GGPs--- restrict comps instead
  ideal_selfP_walrus()
  
  ## Probs given observables
  age_blur_probs_walrus()
  
  # Housekeeping: useful to keep some stuff after function exits. 
  # Trust me, this works...
  ## WRONG VARS FOR WALRUS, OF COURSE!
  list2env( mget( cq( 
    fadsurv, psicle, Nfad_ystart, RoI, Nfad_y, fec_a)),
    envir=environment( sys.function()))
  
  if( 'probs' %in% want){
    # Only OBSERVABLE-conditioned ones, to avoid confusing despack
    returnables <- c( returnables, returnList(
      Pr_MOP_EYEYL,
      Pr_XmHSP_EYEY,
      Pr_selfP_EYDY
    ))
  }
  want <- want %except% 'probs'
  if( !length( want)){
return( returnables) # that's all folks  
  }
  
  # lglk (if there's any data, per kin-type)
  # Not an "A" to be seen below; 
  # all now with E(stimated age)
  # Can turn off any or all kin-types
  # by setting their n_KINP_subs to NULL
  # This allows checking Dlglk closeness-to-zero per kin-type

  lglk <- 0

  if( length( n_MOP_EYEYL)>0){ # use eg NULL for no-data  
    lglk <- lglk + ldpois( n_MOP_EYEYL, 
        n_comp_MOP_EYEYL * Pr_MOP_EYEYL) 
  }
  
  if( length( n_XmHSP_EYEY)>0){ # use eg NULL for no-data
    lglk <- lglk + ldpois( n_XmHSP_EYEY, 
        n_comp_XmHSP_EYEY * Pr_XmHSP_EYEY)
  }
  
  if( length( n_selfP_EYDY)>0){ # use eg NULL for no-data    
    lglk <- lglk + ldpois( n_selfP_EYDY, 
        n_comp_selfP_EYDY * Pr_selfP_EYDY)
  }

  if( identical( want, 'just_lglk')){
    returnables <- c( lglk) # remove 1x1 dimensions, return scalar
  } else {
    returnables$lglk <- c( lglk) # no dims please
  }
  
return( returnables)
}
