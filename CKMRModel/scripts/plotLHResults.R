library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(PNWColors)

load("./results/design_df_filled.RData")

# take out scenario 9
design_df <- design_df[-which(design_df$S == "S9"),]

# filter to keep L1 and L3
design_df <- design_df %>% filter(L == "L1" | L == "L3")

lhdf <- filter(design_df, Value == "Est") %>%
  select(ID, CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "Est") %>%
  filter(Par == "fadsurv") %>%
  group_by(CKMR, L) %>%
  summarize(Mean = round(mean(Est, na.rm=TRUE), digits = 2)) 

cvdf <- na.omit(filter(design_df, Value == "CV")) %>%
  select(ID, CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "CV") %>%
  mutate(CKMR = as.factor(CKMR)) %>%
  mutate(L = factor(L, levels=c("L1", "L3"), labels=c("No", "Yes")))


sedf <- na.omit(filter(design_df, Value == "SE")) %>%
  select(ID, CKMR, Value, D, L, S, fadsurv, fjusurv, ppn_breedy) %>%
  pivot_longer(cols = c(fadsurv, fjusurv, ppn_breedy), names_to = "Par", 
               values_to = "SE") %>%
  mutate(CKMR = as.factor(CKMR)) %>%
  mutate(L = factor(L, levels=c("L1", "L3"), labels=c("No", "Yes")))


### REsULTS TEXT
# mean decrease in SE on fadsurv with CKMR
sedf %>% filter(Par == "fadsurv") %>% 
  pivot_wider(names_from = CKMR, values_from = SE) %>%
  mutate(Gain = Yes - No) %>%
  summarize(round(mean(Gain), digits = 2))

# mean expected CVs under different demo scenarios
sedf %>% filter(Par == "fadsurv") %>% 
  group_by(D) %>% 
  summarize(Mean = round(mean(SE), digits = 2), 
            Min = round(min(SE), digits = 2), 
            Max = round(max(SE), digits = 2))

# mean decrease in fjusurv with CKMR
sedf %>% filter(Par == "fjusurv") %>% 
  pivot_wider(names_from = CKMR, values_from = SE) %>%
  mutate(Gain = Yes - No) %>%
  summarize(round(mean(Gain), digits = 2))

# mean expected SEs under different demo scenarios
sedf %>% filter(Par == "fjusurv") %>% 
  group_by(D) %>% 
  summarize(Mean = round(mean(SE), digits = 2), 
            Min = round(min(SE), digits = 2), 
            Max = round(max(SE), digits = 2))

# mean decrease in ppnbreedy with CKMR
sedf %>% filter(Par == "ppn_breedy") %>% 
  pivot_wider(names_from = CKMR, values_from = SE) %>%
  mutate(Gain = Yes - No) %>%
  summarize(round(mean(Gain), digits = 2))

# mean expected CVs under different demo scenarios
sedf %>% filter(Par == "ppn_breedy") %>% 
#  group_by(D) %>% 
  summarize(Mean = round(mean(SE), digits = 2), 
            Min = round(min(SE), digits = 2), 
            Max = round(max(SE), digits = 2))


# Appendices

lhappendix <- na.omit(filter(design_df, Value == "SE")) %>%
  select(D, L, S, CKMR, fadsurv, fjusurv, ppn_breedy) %>%
  mutate(L = factor(L, levels=c("L3", "L1"), labels=c("Yes", "No"))) %>%
  rename("Lethal Samples" = L, "Sampling Scenario" = S,
         "Demographic Scenario" = D) %>%
  mutate(fadsurv = round(fadsurv, digits = 2),
         fjusurv = round(fjusurv, digits = 2),
         ppn_breedy = round(ppn_breedy, digits = 2)) %>%
  rename("Adult Female Survival" = fadsurv, 
         "Juvenile Female Survival" = fjusurv,
         "P. Adult Female in State 2" = ppn_breedy)

write.csv(lhappendix, file = "./results/lhappendix.csv")

Nappendix <- na.omit(filter(design_df, Value == "CV")) %>%
  select(D, L, S, CKMR, Nfad_2015, Nfad_2020, Nfad_2025) %>%
  mutate(L = factor(L, levels=c("L3", "L1"), labels=c("Yes", "No"))) %>%
  rename("Lethal Samples" = L, "Sampling Scenario" = S,
         "Demographic Scenario" = D) %>%
  mutate(Nfad_2015 = round(Nfad_2015, digits = 2),
         Nfad_2020 = round(Nfad_2020, digits = 2),
         Nfad_2025 = round(Nfad_2025, digits = 2)) %>%
  rename("2015 Adult Females" = Nfad_2015, 
         "2020 Adult Females" = Nfad_2020,
         "2025 Adult Females" = Nfad_2025)

write.csv(Nappendix, file = "./results/Nappendix.csv")
