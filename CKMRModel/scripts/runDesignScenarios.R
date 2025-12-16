
# script to run fit_walrus_ckmr and design_eg

Ds <- c("D1", "D2", "D3")
Ls <- c("L1", "L3")
Ss <- paste0("S", 1:9)

scenarios <- data.frame(expand.grid(Ds, Ls, Ss))

# load empty dataframes to hold design results
source("./scripts/createArrays.R")

for (i in 1:nrow(scenarios)){
 
  TEST_CASE <- paste0(scenarios[i,1], "_", scenarios[i,2], "_", scenarios[i,3])
   
  BASE_CASE <- paste0(scenarios[i,1], "_", scenarios[i,2], "_S0")
  
  d <- as.numeric(scenarios[i,1])
  l <- as.numeric(scenarios[i,2])
  s <- as.numeric(scenarios[i,3])
  
  source("./Scripts/fit_walrus_ckmr_07_01_scenarios.R")
  source("./Scripts/design_eg_ekj.R")
}

save(design_df, file = "./results/design_df_filled.RData")
