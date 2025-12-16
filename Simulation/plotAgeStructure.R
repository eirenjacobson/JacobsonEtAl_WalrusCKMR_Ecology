load("./simulation/WalrusSim_RealisticNoLethality.RData")

library(tidyr)
library(dplyr)
library(ggplot2)

# age structure at end of time period

fagestructure <- indiv %>% 
  filter(Sex == "F") %>%
  filter(is.na(DeathY)) %>%
  select(AgeLast) %>%
  table() %>%
  as.data.frame() %>%
  mutate("Perc" = Freq/sum(Freq))

ggplot(data = agestructure) +
  geom_col(aes(x=AgeLast, y=Perc))

ageclasses <- indiv %>%
  filter(Sex == "F") %>%
  filter(is.na(DeathY)) %>%
  select(AgeLast) %>%
  mutate("DevStage" = ifelse(AgeLast %in% 1:6, "Juv", "Adult")) %>%
  count(DevStage) 
