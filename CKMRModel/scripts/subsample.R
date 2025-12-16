
# function to subsample from a scenario labelled DX_LX_S0
# where S0 indicates that sampling occurred in all possible years 2023-2028

library(tidyr)
library(dplyr)

subsample <- function(suffix, samplingscenario){
  
  load(paste0("./simulation/AllWalrusSamples_", suffix, ".RData"))
  
  if(samplingscenario == 1){ # leave out 2024, 2026, 2027, 2028
    
    allsamples <- allsamples[-which(allsamples$SampY %in% c(2024, 2026, 2027, 2028)),]
    
  } # end s1
  
  if(samplingscenario == 2){ # leave out 2024, 2027 and 2028
    
    allsamples <- allsamples[-which(allsamples$SampY %in% c(2024, 2027, 2028)),]
    
  } # end s2
  
  if(samplingscenario == 3){ # leave out 2024, 2028 
    
    allsamples <- allsamples[-which(allsamples$SampY %in% c(2024, 2028)),]
    
  } # end s3
  
  if(samplingscenario == 4){ # no effort in 2024
    
    allsamples <- allsamples[-which(allsamples$SampY %in% c(2024)),]
    
  } # end s4
  
  if(samplingscenario == 5){ # leave out 2024, 2026, 2027, 2028 and 0.75 effort in 2025
    s2025 <- allsamples[which(allsamples$SampY == 2025 & allsamples$Lethality == "NONLETHAL"),]
    c25 <- round(nrow(s2025)*.75)

    allsamples <- rbind.data.frame(
      allsamples[-which(allsamples$SampY %in% c(2024:2028)),], 
      s2025[1:c25,])
  } # end s5
  
  if(samplingscenario == 6){ # leave out 2024, 2027, 2028 and 0.75 effort in 2025 and 2026
    s2025 <- allsamples[which(allsamples$SampY == 2025),]
    s2026 <- allsamples[which(allsamples$SampY == 2026),]
    c25 <- round(nrow(s2025)*.75)
    c26 <- round(nrow(s2026)*.75)
    
    allsamples <- rbind.data.frame(
      allsamples[-which(allsamples$SampY %in% c(2024:2028)),], 
      s2025[1:c25,], s2026[1:c26,])
  } # end s6
  
  if(samplingscenario == 7){ # no effort in 2024 or 2028, 75% effort in 2025, 2026, 2027
    
    s2025 <- allsamples[which(allsamples$SampY == 2025),]
    s2026 <- allsamples[which(allsamples$SampY == 2026),]
    s2027 <- allsamples[which(allsamples$SampY == 2027),]
    c25 <- round(nrow(s2025)*.75)
    c26 <- round(nrow(s2026)*.75)
    c27 <- round(nrow(s2027)*.75)
    
    allsamples <- rbind.data.frame(allsamples[-which(allsamples$SampY %in% 2024:2028),],
                                   s2025[1:c25,], s2026[1:c26,], s2027[1:c27,])
  } # end s7
  
  if(samplingscenario == 8){ # no effort in 2024, 75% effort in 2025:2028
    
    s2025 <- allsamples[which(allsamples$SampY == 2025),]
    s2026 <- allsamples[which(allsamples$SampY == 2026),]
    s2027 <- allsamples[which(allsamples$SampY == 2027),]
    s2028 <- allsamples[which(allsamples$SampY == 2028),]
    c25 <- round(nrow(s2025)*.75)
    c26 <- round(nrow(s2026)*.75)
    c27 <- round(nrow(s2027)*.75)
    c28 <- round(nrow(s2028)*.75)
    
    allsamples <- rbind.data.frame(allsamples[-which(allsamples$SampY %in% 2024:2028),],
                                   s2025[1:c25,], s2026[1:c26,], s2027[1:c27,], s2028[1:c28,])
  } # end s8
  
  if(samplingscenario == 9){ # no effort 2026:2028
    
    allsamples <- allsamples[-which(allsamples$SampY %in% 2026:2028),]
    
  } # end s9
  
  # now edit allsamples so each individual only appears once in samples
  # create samples from allsamples
  metab <- table( allsamples$Me) 
  multicaps <- names( metab[ metab>1]) # find IDs that were recaptured
  # divide into singles and multiples
  allsamples_nok <- allsamples[which(allsamples$Me %in% multicaps),]
  allsamples_ok <- allsamples[-which(allsamples$Me %in% multicaps),] 
  for (i in unique(multicaps)){
    lastcap <- max(which(allsamples_nok$Me == i))
    allsamples_ok <- rbind.data.frame(allsamples_ok, allsamples_nok[i,])
  } # end for multicaps
    
  samples <- arrange(allsamples_ok, BirthY)
    
  # save outputs for samples and allsamples
  save(samples, file = paste0("./simulation/WalrusSamples_", 
                              substr(suffix, 1,6), "S", samplingscenario,
                              substr(suffix, 9, 19), ".RData"))
  save(allsamples, file = paste0("./simulation/AllWalrusSamples_", 
                              substr(suffix, 1,6), "S", samplingscenario, 
                              substr(suffix, 9, 19), ".RData"))
  
} # end function