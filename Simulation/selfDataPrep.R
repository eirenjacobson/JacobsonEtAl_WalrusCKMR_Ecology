
library(fastmatch)
library(dplyr)
library(offarray)
library(mvbutils)

# USER INPUTS #
maxAge <- 44
years <- 2000:2028 # birth years of potential offspring
syears <- c(2013:2028) # years in which samples were collected
Amat <- 7
# END USER INPUTS #

load("./simulation/AllWalrusSamples.RData")

allsamples$Index <- 1:nrow(allsamples)
# find matches
mm <- fmatch(allsamples$Me, allsamples$Me)

recaps <- allsamples[mm != allsamples$Index,] |>
  select(Me, Sex, SampY, AgeLast) |>
  rename(Y = SampY, A = AgeLast)

Samps <- allsamples |>
  select(Me, Sex, SampY, AgeLast) |>
  rename(Y = SampY, A = AgeLast)

extract.named(Samps)

# find indices of matching pairs
for (i in unique(mm)){
  
  # pull out samples that match 
  ss <- allsamples[mm == i,]
  
  if(nrow(ss) == 1){next} else {
  # need to get all possible pairs
  combos <- t(combn(1:nrow(ss), 2))
  
  if(exists("selfrecaps")==FALSE){
  selfrecaps <- matrix(c(ss$Index[combos[,1]], ss$Index[combos[,2]]), ncol = 2)} else {
    selfrecaps <- rbind(selfrecaps, 
                        matrix(c(ss$Index[combos[,1]], ss$Index[combos[,2]]), ncol = 2))}
  
}}



m_self_YA <- offarray(table(Y=Y, A=A))

Rs_range <- Rju_range <- cq( Russia, USA)
As_range <- 0:maxAge 
Ys_range <- syears 

n_comp_SELF_RAY <- autoloop(
  rs1 = Rs_range, as1=As_range, ys1 = Ys_range, 
  rs2 = Rs_range, as2 = As_range, ys2 = Ys_range,{
    bs1 <- ys1 - as1
    bs2 <- ys2 - as2
    (rs1=='USA') * (rs2=='USA') *
      (bs1 == bs2) *
      m_self_YA[ys1, as1] * m_self_YA[ys2, as2] 
  })


Region <- rep( 'USA', nrow(selfrecaps))
# number of pairs
n_SELF_RAY <- offarray(table(
  rs1 = Region,
  as1 = A[selfrecaps[,1]],
  ys1 = Y[selfrecaps[,1]],
  rs2 = Region,
  as2 = A[selfrecaps[,2]],
  ys2 = Y[selfrecaps[,1]]
), template = n_comp_SELF_RAY)


self_out <- returnList( n_comp_SELF_RAY, n_SELF_RAY) # auto namer
save(self_out, file = "SELF.RData")
