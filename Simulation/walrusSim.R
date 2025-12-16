# Walrus Simulation

library(fishSim)
library(dplyr)
library(tidyr)

seeds <- c(20241023, 20241121, 29571492, 76937593, 92759173, 
           41850183, 38519472, 35719375, 19276639, 71853917)

#seeds2 <- runif(40, 10000000, 99999999)
#seeds2 <- round(seeds2, digits = 0)
#save(seeds2, file = "./simulation/seeds2.RData")

#for (s in seeds2){

# clear allsamples (bc it get added to...)
rm(samples)
rm(allsamples)

s <- 20241023
set.seed(s)
# filename suffix
suffix <- paste0("D3_L3_S0_Sd", s)

source("./simulation/walrusFounders.R")
source("./simulation/walrusMate.R")
source("./simulation/walrusCapture.R")
source("./simulation/walrusPsicle.R")
 
# USER INPUTS
senescence <- TRUE # do females have a post-reproductive phase?
lethality <- TRUE # are lethal samples collected?
ageselect <- TRUE # whether or not to target specific juvenile ages (if FALSE, target devstage)

# Life-history parameters
maxAge <- 37 # nobody survives beyond this age
AFR_F <- 6 # age of first reproduction for females 
ALR_F <- 29 # age of last reproduction (if senescence == T)
AFR_M <- 15 # age of first reproduction for males
S0 <- 0.7 # total calf survival 
S1 <- 0.925 # juvenile survival 
if (senescence == FALSE) {S2 <- 0.9622} else { # fudged this
  S2 <- 0.99 # reproductive adult F survival 6+
  S3 <- 0.6} # post-reproductive adult F survival 30+
fec <- 0.3 # 0.13 is rate per repro adult F per year (not relevant with psicle)
pop <- 250000 # population size 
# nsamples <- 1500 # number of captures per year divided by 10 
years <- 1950:2030
syears <- c(2013:2017, 2023:2028) # years in which lethal and/or nonlethal samples were collected
if(lethality == TRUE){lsyears <- c(2012:2014, 2023:2028)} else {lsyears <- NA} # years in which lethal samples were collected
fsyears <- 2023:2028 # future sample years (where targets will matter)
# for future samples, targets are:
ntarget_byage <- c(310, 160, 130, 100, 100, 600) # target number of nonlethal samples for fsyears per age class
ntarget_bystage <- c(800, 600) # if nonselective within devstage
if(lethality==TRUE){nltarget <- c(50, 50, 50, 500, 500, 500, 500, 500, 500)} else {nltarget <- 0} # number of lethal samples for lsyears for adult F only
psi2 <- 0.1 # probability of breeding at 2-yr interval
psi3 <- 0.5 # probability of breeding at 3+ yr interval (+1st repro)
pUSA <- 0.5
pRUS <- 0.5

# END USER INPUTS

# ORGANIZE SAMPLE SIZES

# target number of samples in age classes 0, 1, 2, 3, 4-5, 6+
# load info on past samples:
load("./samplesizes/historicalsamplesummary.RData") 
# create a dataframe with info on historical and future samples
# targets for future samples
if(ageselect == TRUE){
newdf <- data.frame(Sex = "B", Year = fsyears, 
                    "0" = ntarget_byage[1], "1" = ntarget_byage[2], "2" = ntarget_byage[3], 
                    "3" = ntarget_byage[4], "4to5" = ntarget_byage[5], A = ntarget_byage[6])
names(newdf) <- names(df)[1:8]

#df[,3:7] <- df[,3:7]*10 # increase sample size by factor of 10

sampledf <- df %>% 
  select(-U) %>%
  rbind.data.frame(newdf) %>% 
  pivot_longer(cols = `0`:A, names_to = "AgeClass", values_to = "Number") 
} else {
  
  newdf <- data.frame(Sex = "B", Year = fsyears, 
                      J = ntarget_bystage[1], A = ntarget_bystage[2])
  
  stagedf <- df %>% 
    mutate("J" = (`0`+`1`+`2`+`3`+`4to5`)) %>%
    select(Sex, Year, J, A)
  
  sampledf <- rbind.data.frame(newdf, stagedf) %>%
    pivot_longer(cols = J:A, names_to = "AgeClass", values_to = "Number")
    
}
# END SAMPLE SIZES

# SET UP VARIOUS LIFE-HISTORY TABLES

# set up "movement" matrix for breeding psicle
moveMat <- matrix(data = c(0, 1, 0, 
                           psi2, 0, 1-psi2, 
                           psi3, 0, 1-psi3),
                  nrow = 3, byrow = TRUE)

# generate vectors of age and sex-specific rates
# fecundity (currently set to knife-edge)
# length = maxAge + 1 (one fecundity per age class)
# TODO add male curve so males must be 15
if(senescence == TRUE){
  # note this isn't really needed given psicle, but doesn't hurt either
  femaleCurve <- c(rep(0, AFR_F), rep(1, ALR_F - AFR_F), 
                   rep(0, maxAge - ALR_F + 1))} else
                   {femaleCurve <- c(rep(0, AFR_F), rep(1, maxAge - AFR_F), 0)}
maleCurve <- c(rep(0, AFR_M), rep(1, maxAge - AFR_M + 1))
# survival by age
# length = maxAge
if(senescence == TRUE){
  sCurve <- c(S0, rep(S1, AFR_F-1), rep(S2, ALR_F-AFR_F), rep(S3, maxAge-ALR_F))} else {
  sCurve <- c(S0, rep(S1, AFR_F-1), rep(S2, maxAge-AFR_F+1))
  }

# cumulative survival by age
csCurve <- rep(0, maxAge)
csCurve[1] <- S0 

for (i in 2:maxAge){
  csCurve[i] <- prod(sCurve[1:i])
}

sc <- csCurve/sum(csCurve)
# mortality by age
ageMort <- 1-sCurve

# get equil proportions give psi vals
Tbreedf <- rbind( 
  c( 0, psi2, psi3),
  c( 1, 0, 0),
  c( 0, 1-psi2, 1-psi3)
)

# No longer going the eigenroute; get stationary eivec by...
M <- Tbreedf - diag( 3)
M[3,] <- 1
Pbreedf <- solve( M, c( 0,0,1))

# END OF LIFE-HISTORY TABLE SET UP

# INITIALIZE SIMULATION

if(senescence == FALSE){
indiv <- walrusFounders(pop = pop,
                      stocks = Pbreedf,
                      maxAge = maxAge,
                      survCurv = sc,
                      year1 = years[1],
                      independenceAge = 4,
                      regions = c(pUSA, pRUS))} else {
          indiv <- walrusFounders(pop = pop,
                                  stocks = Pbreedf,
                                  maxAge = maxAge,
                                  survCurv = sc,
                                  year1 = years, 
                                  senescence = TRUE)                        
                      } # end else

# INITIALIZE TRACKERS
                        
ncaptures <- list()
nlcaptures <- list()
neligible <- list()
nalive <- list()
nadfemales <- list()

# RUN SIMULATION

for (i in 1:length(years)){
  
  indiv <- walrusMate(indiv,
                   year = years[i],
                   type = "ageSex",
                   femaleCurve = femaleCurve,
                   maleCurve = maleCurve,
                   batchSize = 1, # fecundity
                   fecundityDist = "binomial",
                   maxClutch = 1,
                   firstBreedF = AFR_F)
  
  nadfemales[[i]] <- length(which(is.na(indiv$DeathY) & indiv$Sex == "F" & indiv$AgeLast %in% AFR_F:maxAge)) # tracker
  nalive[[i]] <- length(which(is.na(indiv$DeathY))) # tracker
  
  # captures
  if (years[i] %in% syears){ 
    
    source("./simulation/captureCode.R")  # this is really ugly so moved to its own file
    
  } # end syears
  
  ncaptures[[i]] <- length(which(indiv$SampY == years[i] & indiv$Lethality == "NONLETHAL")) # tracker
  nlcaptures[[i]] <- length(which(indiv$SampY == years[i] & indiv$Lethality == "LETHAL"))

  indiv <- mort(indiv,
                year = years[i],
                maxAge = maxAge,
                type = "age",
                ageMort = c(ageMort, 1))
  
  # also kill off any calves whose mothers died (naturally or via lethal sampling)
  
  deadmothers <- which(indiv$DeathY == years[i] & indiv$Stock == 2)
  deadbabies <- which(indiv$Mum %in% indiv$Me[deadmothers] & indiv$BirthY == years[i])
  indiv[deadbabies,]$DeathY <- years[i]
  # OK here
  indiv <- psicle(indiv, moveMat)

  indiv <- birthdays(indiv)
  
  # put 4-year-old females into Stage 3, so they might get pregnant next year and have 1st calf at 6
  indiv[which(indiv$AgeLast == 4 & is.na(indiv$DeathY) & indiv$Sex == "F"),]$Stock <- 3
  
  # 4 yr olds must become independent
  indiv[which(indiv$AgeLast == 4 & is.na(indiv$DeathY)),]$Independent <- 1
  
  # if senescence == TRUE, put females into post-reproductive class when they hit 30
  if(senescence == TRUE){indiv[which(indiv$AgeLast == 30 & is.na(indiv$DeathY) & indiv$Sex == "F"),]$Stock <- 0}
  
   
#  if(years[i] == 2000){indiv2000 <- indiv; save(indiv2000, file = "./simulation/indiv2000.RData")}
#  if(years[i] == 2012){indiv2012 <- indiv; save(indiv2012, file = paste0("./simulation/results/indiv2012_", suffix, ".RData"))}

} # end for i

# END SIMULATION

# checks that numbers come out right
#unlist(nalive)
# note this may not be exactly equal to nsamples due to possible recaptures
#unlist(ncaptures)

#save(indiv, file = paste0("./simulation/results/WalrusSim_", suffix, ".RData")) # too big to commit!

samples <- indiv[which(!is.na(indiv$SampY)),]

save(samples, file = paste0("./simulation/results/WalrusSamples_", suffix, ".RData"))
save(allsamples, file = paste0("./simulation/results/AllWalrusSamples_", suffix, ".RData"))

future <- allsamples %>% 
  filter(SampY %in% fsyears) %>%
  mutate(Region = "USA") %>%
  rename(Year = SampY, Age = AgeLast) %>%
  group_by(Region, Age, Sex, Year) %>% 
  count() %>%
  rename(Number = n)

save(future, file = paste0("./simulation/results/futuresamplesizes_", suffix, ".RData"))

# check that sample numbers are correct
#View(allsamples %>% filter(Sex == "F") %>% group_by(SampY,AgeLast) %>% count())

#View(allsamples %>% mutate(DevStage = ifelse(AgeLast<6, "Juv", "Ad")) %>%
#       group_by(SampY, DevStage) %>% count())

Nfad <- unlist(nadfemales)
Nalive <- unlist(nalive)

# calculate growth rate
#mean(diff(log(Nfad[51:length(Nfad)])))

out <- c("Nfad_2000" = Nfad[51], 
         "RoI" = mean(diff(log(Nfad[51:length(Nfad)]))))

save(out, file = paste0("./simulation/results/Nfad_RoI_", suffix, ".RData"))

save(Nfad, file = paste0("./simulation/results/Nfad_", suffix, ".RData"))
#save(Nalive, file = paste0("./simulation/results/Nalive_", suffix, ".RData"))

#} # end seeds


# 
# metab <- table( allsamples$Me)
# selfptab <- names(which(metab > 1))
# 
# 
# selfrecaps <- allsamples %>% 
#   filter(Me %in% selfptab) %>%
#   select("Me", "SampY") 
# 
# selfcaptab <- data.frame()
# c <- 1
# for (i in unique(selfrecaps$Me)){
#   
#   df <- filter(selfrecaps, Me == i)
#   selfcaptab[c,1] <- i
#   selfcaptab[c,2] <- min(df$SampY)
#   selfcaptab[c,3] <- max(df$SampY)
#   selfcaptab[c,4] <- nrow(df)
#   c <- c+1
# }
# 
# filter(selfcaptab, V4 > 1)
