"tsep" <- function( yj, yc, wyj){
  ## Tsep = #(indept movement-years) 
  ## between juve (potoff) sampled at yj and cow at yc
  ## given juve was weaned in year wyj (not age)
  
  r"--{
  ## For MOP (and perhaps HSP), we'll need Prob that an M and her O are in certain places at certain times. Key insight: this is the _same_ as prob that *one* individual will be in those two places given a time-gap equal to the amount of movement-years there's been. The latter is 's affected by when weaning happened (possibly after both are sampled).

  Let wyj be weaning year. The if( FALSE) code below is the long-winded way to figure "effective years of movement". More succinctly, that's abs( yc-wyj) + abs( yj-wyj) *unless* wyj exceeds both yj & yc--- in which case, it's just abs( yj-yc). I believe this is correct (checked pictorially... but I am not the most thorough) but need a diagram here!
  
  This was originally a (large) precalced array in R (and prolly TMB) but makes more sense as an on-the-spot function; in R it will get yj,yc,wyj "in parallel" from autoloop
  }--"
  
  ifelse(  (wyj >= yj) & ( wyj >= yc),
           abs( yc-yj),
           abs(yc-wyj) + abs(yj-wyj)
  )
}