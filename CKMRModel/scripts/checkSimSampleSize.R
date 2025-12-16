# script to check that sample sizes match targets set in 
# Sample_Sizes_ForEiren_2025_03_05.xlsx

library(readxl)

# historical samples sizes
historicalsamples <- read_excel("./simulation/Sample_Sizes_ForEiren_2025_03_05.xlsx", sheet = 1) %>% 
  rename("Number" = "SimUniqueIndividuals", "AgeClass" = "Age") %>%
  mutate(Type = "Live")


# future sample sizes -- no lethal

sheet <- ifelse(lethality == TRUE, 3, 2)

futuresamples <- read_excel("./simulation/Sample_Sizes_ForEiren_2025_03_05.xlsx", sheet = sheet) %>%
  pivot_longer(cols = 3:8, names_to = "AgeClass") %>%
  rename("Number" = "value") %>%
  expand_grid(Year = fsyears) %>%
  select(Year, AgeClass, Sex, Number, Type)

sampledf <- rbind.data.frame(historicalsamples, futuresamples)
syears <- c(hsyears, fsyears, lsyears)