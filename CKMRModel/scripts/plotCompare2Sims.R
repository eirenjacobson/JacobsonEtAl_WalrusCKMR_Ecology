# code for plotting stuff in compare2sims.R

library(dplyr)
library(tidyr)
library(ggplot2)

#

as.data.frame(Dlglk0) %>% 
  mutate(Parameter = row.names(as.data.frame(Dlglk0))) %>%
  pivot_longer(cols = 1:50, names_to = "SimID", values_to = "Dlglk0") %>%
  group_by(Parameter) %>%
  summarize(Mean = mean(Dlglk0), Var = var(Dlglk0))
  


# plot observed v expected numbers of kin pairs

Obs <- as.data.frame(compcheck[,SLICE = "Obs",]) 

Obs.long <- Obs %>% 
  mutate("Type" = row.names(Obs)) %>%
  pivot_longer(cols = 1:50, names_to = "SimID", values_to = "Obs")

Exp <- as.data.frame(compcheck[,SLICE="Exp",])

Exp.long <- Exp %>%
  mutate("Type" = row.names(Exp)) %>%
  pivot_longer(cols = 1:50, names_to = "SimID", values_to = "Exp") %>%
  group_by(Type) %>%
  summarize(Mean = mean(Exp)) %>%
  mutate(LExp = Mean - 2*sqrt(Mean),
         UExp = Mean + 2*sqrt(Mean))

ggplot() +
  geom_histogram(data = Obs.long, aes(x=Obs), fill = "gray") +
  geom_vline(data=Exp.long, aes(xintercept=Mean), color = "red")+
  geom_vline(data=Exp.long, aes(xintercept = LExp), color = "red", lty = 2) +
  geom_vline(data=Exp.long, aes(xintercept = UExp), color = "red", lty = 2)+
  facet_wrap(~Type, scales = "free", nrow = 3) +
  ylab("Number of simulated datasets")+
  xlab("Number of kin pairs")+
  theme_bw()

ggsave(plot = last_plot(), file = "./figures/ObsVExpKinPairs_wBounds.png", 
       width = 8, height = 6, units = "in")

# plot observed v expected maternal age

df <- data.frame("Age" = 1:37,
                 "Observed" = as.numeric(fcOnl),
                 "Expected" = as.numeric(fcEnl)) %>%
  pivot_longer(cols = 2:3, names_to = "Type", values_to = "Value")
ggplot(df) +
  geom_point(aes(x=Age, y= Value, color = Type))+
  ylab("Number of mother-offspring pairs")+
  theme_bw()

