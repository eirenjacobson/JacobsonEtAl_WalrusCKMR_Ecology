
# Script to take a simulated dataset with sampling from 2023-2028
# and subsample it according to various possible scenarios
# as detailed in the sampling tab of scenarios.xlsx

source("./scripts/subsample.R")

Ls <- c("L1", "L3")
Ds <- c("D1", "D2", "D3")

sdf <- expand.grid(Ds, Ls)

for (i in 1:nrow(sdf)){

  suffix <- paste0(sdf[i,1], "_", sdf[i,2], "_S0_Sd20241023")

  subsample(suffix, 1)
  subsample(suffix, 2)
  subsample(suffix, 3)
  subsample(suffix, 4)
  subsample(suffix, 5)
  subsample(suffix, 6)
  subsample(suffix, 7)
  subsample(suffix, 8)
  subsample(suffix, 9)
  
}
