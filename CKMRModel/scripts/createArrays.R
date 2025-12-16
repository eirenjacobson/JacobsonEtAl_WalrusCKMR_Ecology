# create arrays to fill

library(tidyr)

compdesign <- offarray(rep(NA, 3*2*9*2*2*7*3),
  dimseq = list(D = 1:3, L = c("1", "3"), S = 1:9, 
                CKMR = c("Yes", "No"), SelfP = c("Yes", "No"),
                Type = c("Nfad_2015", "Nfad_2020", "Nfad_2025", 
                         "fadsurv", "fjusurv", "ppn_breedy",
                         "rel_ad_15_25"),
                Value = c("True", "Est", "SE")))


IDs <- rep(NA, 9*3*2)
for (i in 1:length(IDs)){
  IDs[i] <- paste0(scenarios[i,1], "_", scenarios[i,2], "_", scenarios[i,3])
}

design_df <- expand_grid(cbind(factor(IDs), scenarios), 
                         CKMR = factor(c("Yes", "No")), Value = factor(c("TRUE", "Est", "SE", "CV")),
                         Nfad_2015 = 0, Nfad_2020 = 0, Nfad_2025 = 0, 
                         fadsurv = 0, fjusurv = 0, ppn_breedy = 0, 
                         rel_ad_15_25 = 0)

names(design_df)[1:4] <- c("ID", "D", "L", "S")
