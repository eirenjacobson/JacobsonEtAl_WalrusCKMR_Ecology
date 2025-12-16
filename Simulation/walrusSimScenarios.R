# Walrus Simulation

library(fishSim)
library(dplyr)
library(tidyr)
library(readxl)

set.seed(20241023)

source("./simulation/walrusFounders.R")
source("./simulation/walrusMate.R")
source("./simulation/walrusCapture.R")
source("./simulation/walrusPsicle.R")

s <- 20241023
set.seed(s)

# USER INPUTS
# see scenarios.xls for life history parameter values for different demo scenarios
senescence <- TRUE # do females have a post-reproductive phase?
lethality <- TRUE # are lethal samples collected?
ageselect <- TRUE # whether or not to target specific juvenile ages (if FALSE, target devstage)
suffix <- paste0("D2_L3_S0_Sd", s) # filename suffix

# clear allsamples (bc it get added to...)
rm(samples)
rm(allsamples)

# Life-history parameters
maxAge <- 37 # nobody survives beyond this age
AFR_F <- 6 # age of first reproduction for females 
ALR_F <- 29 # age of last reproduction (if senescence == T)
AFR_M <- 15 # age of first reproduction for males
S0 <- 0.7 # total calf survival 
S1 <- 0.9 # juvenile survival 
senescence <- FALSE # do females have a post-reproductive phase?
if (senescence == FALSE) {S2 <- 0.9622} else { 
  S2 <- 0.99 # reproductive adult F survival 6+
  S3 <- 0.55} # post-reproductive adult F survival 30+
pop <- 250000 # population size 
years <- 1950:2030
hsyears <- 2013:2017 # previous years in which samples were collected
if(lethality == TRUE){lsyears <- 2023:2028} else {lsyears <- NA} # years in which lethal samples were collected
fsyears <- 2023:2028 # future sample years (where targets will matter)
psi2 <- 0.1 # probability of breeding at 2-yr interval
psi3 <- 0.5 # probability of breeding at 3+ yr interval (+1st repro)
# proportion of animals in USA/Russia
pUSA <- 1
pRUS <- 0

# END USER INPUTS

# SAMPLE SIZES

# historical samples sizes
historicalsamples <- read_excel("./samplesizes/Sample_Sizes_ForEiren_2025_03_05.xlsx", sheet = 1) %>% 
  rename("Number" = "SimUniqueIndividuals", "AgeClass" = "Age") %>%
  mutate(Type = "Live")

# future sample sizes -- no lethal

sheet <- ifelse(lethality == TRUE, 3, 2)

futuresamples <- read_excel("./samplesizes/Sample_Sizes_ForEiren_2025_03_05.xlsx", sheet = sheet) %>%
  pivot_longer(cols = 3:8, names_to = "AgeClass") %>%
  rename("Number" = "value") %>%
  expand_grid(Year = fsyears) %>%
  select(Year, AgeClass, Sex, Number, Type)

sampledf <- rbind.data.frame(historicalsamples, futuresamples)
syears <- c(hsyears, fsyears, lsyears)

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
                   maxClutch = 1)
  
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
  
} # end for i

# END SIMULATION

# checks that numbers come out right
#unlist(nalive)
# note this may not be exactly equal to nsamples due to possible recaptures
#unlist(ncaptures)

#save(indiv, file = "./simulation/WalrusSim.RData") # too big to commit!

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
#View(allsamples %>% group_by(SampY, Sex, AgeLast) %>% count())

# calculate growth rate

m <- glm(unlist(nalive) ~ years, family = "poisson")

Nfad <- unlist(nadfemales)
Nalive <- unlist(nalive)

out <- c("Nfad_2000" = Nfad[51], 
         "RoI" = mean(diff(log(Nfad[51:length(Nfad)]))))

save(out, file = paste0("./simulation/results/Nfad_RoI_", suffix, ".RData"))

save(Nfad, file = paste0("./simulation/results/Nfad_", suffix, ".RData"))
