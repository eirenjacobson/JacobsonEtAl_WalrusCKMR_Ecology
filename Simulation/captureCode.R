# note possible bug: need to add na.omit in case any of the sexes (F, M, B) don't exist in the sdf

sdf <- subset(sampledf, Year == years[i])

  fdf <- subset(sdf, Sex == "F")
  captured0_f <- walrusCapture(indiv, n = max(fdf$Number[which(fdf$AgeClass == "0")],0), year = years[i], fatal = FALSE, femaleages = 0, maleages = NULL, replace = TRUE, self = TRUE)
  captured1_f <- walrusCapture(captured0_f$indiv, n = max(fdf$Number[which(fdf$AgeClass == "1" & fdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = 1, maleages = NULL, replace = TRUE, self = TRUE)
  captured2_f <- walrusCapture(captured1_f$indiv, n = max(fdf$Number[which(fdf$AgeClass == "2"& fdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = 2, maleages = NULL, replace = TRUE, self = TRUE)
  captured3_f <- walrusCapture(captured2_f$indiv, n = max(fdf$Number[which(fdf$AgeClass == "3" & fdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = 3, maleages = NULL, replace = TRUE, self = TRUE)
  captured45_f <- walrusCapture(captured3_f$indiv, n = max(fdf$Number[which(fdf$AgeClass == "4to5" & fdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = 4:5, maleages = NULL, replace = TRUE, self = TRUE)
  captured6_f <- walrusCapture(captured45_f$indiv, n = max(fdf$Number[which(fdf$AgeClass == "6" & fdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = AFR_F:maxAge, maleages = NULL, replace = TRUE, self = TRUE)
  
  mdf <- subset(sdf, Sex == "M")
  captured0_m <- walrusCapture(captured6_f$indiv, n = max(mdf$Number[which(mdf$AgeClass == "0" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 0, replace = TRUE, self = TRUE)
  captured1_m <- walrusCapture(captured0_m$indiv, n = max(mdf$Number[which(mdf$AgeClass == "1" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 1, replace = TRUE, self = TRUE)
  captured2_m <- walrusCapture(captured1_m$indiv, n = max(mdf$Number[which(mdf$AgeClass == "2" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 2, replace = TRUE, self = TRUE)
  captured3_m <- walrusCapture(captured2_m$indiv, n = max(mdf$Number[which(mdf$AgeClass == "3" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 3, replace = TRUE, self = TRUE)
  captured45_m <- walrusCapture(captured3_m$indiv, n = max(mdf$Number[which(mdf$AgeClass == "4to5" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 4:5, replace = TRUE, self = TRUE)
  captured6_m <- walrusCapture(captured45_m$indiv, n = max(mdf$Number[which(mdf$AgeClass == "6" & mdf$Type == "Live")],0), year = years[i], fatal = FALSE, femaleages = NULL, maleages = 6:44, replace = TRUE, self = TRUE)
  
  indiv <- captured6_m$indiv
  
  if (!exists("allsamples")){allsamples <- rbind.data.frame(captured0_f$samples, captured0_m$samples,
                                                            captured1_f$samples, captured1_m$samples,
                                                            captured2_f$samples, captured2_m$samples,
                                                            captured3_f$samples, captured3_m$samples, 
                                                            captured45_f$samples, captured45_m$samples,
                                                            captured6_f$samples, captured6_m$samples)} else{
                                                              allsamples <- rbind.data.frame(allsamples, captured0_f$samples, captured0_m$samples,
                                                                                             captured1_f$samples, captured1_m$samples,
                                                                                             captured2_f$samples, captured2_m$samples, 
                                                                                             captured3_f$samples, captured3_m$samples, 
                                                                                             captured45_f$samples, captured45_m$samples, 
                                                                                             captured6_f$samples, captured6_m$samples)}

# lethal captures

if (years[i] %in% lsyears){
  
  lethalcaptured6_f <- walrusCapture(indiv, n = max(fdf$Number[which(fdf$AgeClass == "6" & fdf$Type == "Dead")],0), year = years[i], fatal = TRUE, femaleages = AFR_F:maxAge, maleages = NULL, replace = FALSE, self = TRUE)
  lethalcaptured6_m <- walrusCapture(lethalcaptured6_f$indiv, n = max(mdf$Number[which(mdf$AgeClass == "6" & mdf$Type == "Dead")],0), year = years[i], fatal = TRUE, femaleages = NULL, maleages = 6:37, replace = FALSE, self = TRUE)
  
  indiv <- lethalcaptured6_m$indiv
  if(exists("allsamples")){allsamples <- rbind.data.frame(allsamples, lethalcaptured6_f$samples, lethalcaptured6_m$samples)} else (allsamples <- lethalcaptured6_m$samples)
} # end lethal sampling
