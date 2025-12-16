
library(fishSim)
library(dplyr)
library(tidyr)

#source("./simulation/Scripts/walrusMate.R")
source("./simulation/Scripts/walrusCapture.R")

# Life-history parameters 
maxAge <- 40 # nobody survives beyond this age
AFR_F <- 8 # age of first reproduction for females
AFR_M <- 15 # age of first reproduction for males
S1 <- 0.85 # juvenile survival (ages 1-5)
S2 <- 0.95 # adult survival (ages 5+)
fec <- 0.33
pop <- 250000 # population size
nsamples <- 1500 # number of captures per year
nyears <- 50 # number of years to run the simulation for
syears <- 40:50 # years in which samples were collected

# generate vectors of age and sex-specific rates
# fecundity (currently set to knife-edge)
femaleCurve <- c(rep(0, AFR_F), rep(1, (maxAge) - AFR_F), 0)
maleCurve <- c(rep(0, AFR_M), rep(1, (maxAge) - AFR_M), 0)

# survival by age
sCurve <- c(rep(S1, 5), rep(S2, maxAge - 6), 0)
# cumulative survival by age
csCurve <- rep(0, maxAge)
csCurve[1] <- S1 

for (i in 2:maxAge){
  csCurve[i] <- prod(sCurve[1:i])
}

sc <- csCurve/sum(csCurve)
# mortality by age
ageMort <- 1-sCurve

# check that life-history parameters produce reasonable growth rate
check_growthrate(mateType = "ageSex", femaleCurve = femaleCurve, batchSize = fec,
                 mortType = "age", ageMort = ageMort)

indiv <- makeFounders(pop = pop,
                      stocks = 1,
                      maxAge = maxAge,
                      survCurv = sc)

ncaptures <- list()
neligible <- list()
nalive <- list()

for (i in 1:nyears){

  indiv <- altMate(indiv,
                   year = i,
                   type = "ageSex",
                   maleCurve = maleCurve,
                   femaleCurve = femaleCurve,
                   batchSize = 0.33, # fecundity
                   fecundityDist = "binomial",
                   maxClutch = 1)
  
  nalive[[i]] <- length(which(is.na(indiv$DeathY)))

  if (i %in% syears){
    indiv <- walrusCapture(indiv,
                     n = nsamples,
                     year = i,
                     fatal = FALSE,
                     sex = "F", 
                     age = 2:maxAge, 
                     replace = TRUE)
    neligible[[i]] <- length(which(is.na(indiv$DeathY) & indiv$Sex == "F" & indiv$AgeLast %in% 2:maxAge))
    ncaptures[[i]] <- length(which(indiv$SampY == i))
  }

  indiv <- mort(indiv,
                year = i,
                maxAge = maxAge,
                type = "age",
                ageMort = c(ageMort, 1))

  indiv <- birthdays(indiv)

} # end for i

unlist(nalive)
unlist(neligible)
unlist(ncaptures)

q <- quickin(indiv)

pairs <- findRelativesPar(indiv = indiv, sampled = TRUE, nCores = 4)

