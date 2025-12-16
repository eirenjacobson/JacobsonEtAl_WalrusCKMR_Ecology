
library(kinsimmer)
library(mvbutils)
library(offarray)


# USER INPUTS #
maxAge <- 44
years <- 2000:2028 # birth years of potential offspring
syears <- 2013:2028 # years in which samples were collected
Amat <- 7
# END USER INPUTS #

## MVB has special circle of hell for tidyverse stuff!
if( !exists( 'library_tidybollox', mode='function')){
  library_tidybollox <- library
}
library_tidybollox(dplyr)
# library_tidybollox(tidyr) # not needed; pipes in base R now

library(fastmatch)

simdir <- './simulation'
if( !dir.exists( simdir)){
  # ... MVB is not running this from current folder
  simdir <- 'd:/github/walrusCKMR/simulation'
}

source( file.path( simdir, 'ckmr_funs.R'))

walrusQuickin <- try( kinsimmer::quickin) # MVB checking stuff
if( walrusQuickin %is.a% 'try-error'){
  source( file.path( simdir, 'walrusQuickin.R'))
  # creates walrusQuickin()
}


# load simulated data generated from walrusSim 
load( file.path( simdir, "WalrusSamples.RData")) # indiv

#load((file.path(simdir, "AllWalrusSamples.RData"))) # for self-recap

# MVB fix: very important!!!
samples <- within( samples, {
  Dad[ Dad=='founder'] <- NA
  Mum[ Mum=='founder'] <- NA
})


# MVB: get MOPs and MatHSPs. 

qk <- walrusQuickin(samples, max_gen = 1, cousins. = FALSE, thiatics. = FALSE) 
# save( qk, file='postquickin.rda')

# convert to Bravo format
# MVB: I've replaced %>% with |> since pipes are now base-R
Samps <- samples |>
  select(Me, Sex, SampY, AgeLast) |>
  rename(Y = SampY, A = AgeLast)
  
# MVB: FYI non-dplyr alternative with mvbutils
# NB mvbutils::rename.els unfortunately doesn't work with pipes 
# (cozza bad design on my part)
# Samps <- rename.els( Samps[ cq( Me, Sex, SampY, AgeLast)], 
#    AgeLast='A', SampY='Y')

# plonk those 4 things (Me, Y, etc) into the global env
extract.named(Samps)

# find *indices* of matching pairs

if( FALSE){ # NB this is inefficient
  MOPs <- array(rep(NA, nrow(qk$POP)*2), dim = c(nrow(qk$POP), 2))
  for(i in 1:nrow(MOPs)){
    MOPs[i, 1] <- which(Samps$Me == qk$POP[i, 1])
    MOPs[i, 2] <- which(Samps$Me == qk$POP[i, 2])
  }

  HSPs <- array(rep(NA, nrow(qk$MatHSP)*2), dim = c(nrow(qk$MatHSP), 2))
  for (i in 1:nrow(HSPs)){
    HSPs[i, 1] <- which(Samps$Me == qk$MatHSP[i, 1])
    HSPs[i, 2] <- which(Samps$Me == qk$MatHSP[i, 2])
  }
} else {
  # MVB version: MOPs should have same shape as qk$POP

  MOPs <- matrix( match( qk$POP, Me), 
    nrow( qk$POP), 2)

  HSPs <- matrix( match( qk$MatHSP, Me), 
    nrow( qk$MatHSP), 2)
}

# calculate birth years of samples
B <- Y - A

offposs <- B >= years[1] # only consider offspring since year 2000
parposs <- B < (max(years) - Amat) # must have matured by end of sampling

r"--{
MVB: should also check that none of the offposs is a founder, coz those won't show up in a POP. In fact, the date-ranges make this very unlikely or perhaps impossible (there were none when I checked) but let's make sure...
}--"

#print( sum( offposs & is.na( samples$Dad))) # 0, whew
#offposs <- offposs & !is.na( samples$Dad) # in case

MOPs <- MOPs[parposs[MOPs[,1]] & offposs[MOPs[,2]],]
HSPs <- HSPs[offposs[ HSPs[,1]] & offposs[ HSPs[,2]],]

# first year in which offspring could have been born
y0 <- years[1]


# Region placeholder
# MVB: I think we want 2 regions here even if sim only uses one
Rad_range <- Rju_range <- cq( Russia, USA)

Aad_range <- 0:maxAge # 1 to max age (of any sampled animals)
Aju_range <- 0:(syears[length(syears)]-years[1])
Bju_range <- years # birth years of potential offspring
Yad_range <- syears # sampling years (of any samples) 

# m_... is samp size
# array of Y by A for possible parents
# EKJ NOTE: I fixed this array so it contains all possible combinations (not
# just those that appear in the data). This may need to happen for other arrays at some point?
m_ad_YA_temp <- offarray(x = rep(0, length(years)*length(Aad_range)), 
                         dimseq = list(Y = years, A = Aad_range))
m_ad_YA <- offarray(table(Y=Y[parposs], A=A[parposs]), template = m_ad_YA_temp)

# array of Y by A for possible offspring
m_ju_YA <- offarray(table(Y=Y[offposs], A = A[offposs]))

m_ju_A <- offarray(table(A = A[offposs]))
m_ad_A <- offarray(table(A = A[parposs]))

# for half-siblings
m_ad_Y <- offarray( table( Y)) # for noage
m_ju_B <- offarray(table(B = B[offposs]))

# decide which years and ages we will use comparisons from


# number of comparisons
r"--{
MVB: added region check; FOR NOW all samples come ONLY from USA, but historically there are also Rus samples that must be built in. Though (for a while at least) we could just pretend they are USA, given baseline assumption about movement
}--"

n_comp_MOP_RAY <- autoloop(
  rj = Rju_range, aj=Aju_range, yj = Yad_range, 
  rc = Rad_range, ac = Aad_range, yc = Yad_range,{
    bad <- yc - ac # birth year of adult/cow
    bju <- yj - aj # birth year of juvenile
    (rj=='USA') * (rc=='USA') *
    (bju >= (bad+Amat)) *
    m_ju_YA[yj, aj] * m_ad_YA[yc, ac] 
  })

Region <- rep( 'USA', nrow(MOPs))
# number of pairs
n_MOP_RAY <- offarray(table(
  rj = Region,
  aj = A[MOPs[,2]],
  yj = Y[MOPs[,2]],
  rc = Region,
  ac = A[MOPs[,1]],
  yc = Y[MOPs[,1]]
), template = n_comp_MOP_RAY) # template ensures dims match n_comp_MOP_RAY

if( TRUE){ 
  # MVB: We don't neeeed this, though mebbe helpful for display
  "boring_dfize" <-
    function( n_MOP, n_comp_MOP) {
      MOP_df <- as.data.frame( n_comp_MOP, name_of_response='n_comp_MOP')
      temp <- as.data.frame( n_MOP, name_of_response='n_MOP')
      MOP_df <- cbind( MOP_df, n_MOP=temp$n_MOP) # really ought to check that the index rows match...
      MOP_df <- MOP_df %where% (n_comp_MOP > 0)
      return( MOP_df)
    }

  MOP_df <- boring_dfize( n_MOP_RAY, n_comp_MOP_RAY)
}

# And for HSPs...
n_comp_HSP_B1B2 <- autoloop( B1=Bju_range, B2=Bju_range, {
  # NB *exclude* double-count and same-cohort and back-to-back years
  m_ju_B[ B1] * m_ju_B[ B2] * (B2>(B1+1))
})

n_HSP_B1B2 <- offarray( table( B1=B[ HSPs[,1]], B2=B[ HSPs[,2]]),
                        template=n_comp_HSP_B1B2)
# 
# n_comp_POP_noage_BY <- autoloop( Bju=Bju_range, Yad=Yad_range, {
#   m_ju_B[ Bju] * m_ad_Y[ Yad] 
# })
# 
# # Now tot up number of MOPs seen, in the same way
# n_POP_noage_BY <- offarray( table( Bju=B[ MOPs[,2]], Yad=Y[ MOPs[,1]]),
#                             template=n_comp_POP_noage_BY)

# MVB: since we are perforce using mvbutils, might as well use it fully
# out <- list(n_comp_MOP_RAY = n_comp_MOP_RAY, n_MOP_RAY = n_MOP_RAY)
mop_out <- returnList( n_comp_MOP_RAY, n_MOP_RAY) # auto namer
save(mop_out, file = "MOPs.RData")

hsp_out <- returnList(n_comp_HSP_B1B2, n_HSP_B1B2)
save(hsp_out, file = "HSPs.RData")
