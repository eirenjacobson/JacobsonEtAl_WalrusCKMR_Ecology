#######
## Here's some code for Design calculations
## Assumes that the following exists, eg from 'fit_walrus_ckmR_24_11.R'
# ptru, Vpar, denv (which contains sample size info)
# and various functions

# Ready for import into Lyx--- which works pretty well
if( !dir.exists( 'results')){
  mkdir( 'results') # MVB 8/11--- tho might wanna put elsewhere anyway
}


# Need a SS scenario
denv <- environment( lglk_with_data)
m_SYEL_0 <- denv$m_SYEL

basic0 <- doades( 
  m_SYEL_0,          # could be anything
  Hbits,             # *assume* this p'tic Hbits matches this lglk_with_data
  interesting_stuff, # look at function code
  lglk_with_data
)

basic0_noCK <- doades( 
  m_SYEL_0,          # could be anything
  Hbits,             # *assume* this p'tic Hbits matches this lglk_with_data
  interesting_stuff, # look at function code
  lglk_with_data,
  comp_wts= c( MOP=0.01, XmHSP=0.01)
)


# Would be better to use a filename that reflects something about this particular Design... pop-dyn scenario label, SS label
write.csv( basic0, file='results/basic0.csv')
