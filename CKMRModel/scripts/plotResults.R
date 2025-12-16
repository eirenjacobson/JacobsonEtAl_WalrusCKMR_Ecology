
library(ggplot2)
library(dplyr)
library(tidyr)

# check expected versus observed comparison types


for (i in 1:length(suffixes)){
  
  c <- data.frame(compcheck[,,suffixes[i]])
  c$ID <- suffixes[i]
  c$Type <- c("MOP", "XmHSP", "SelfP")
  if(i == 1){allcomparisons <- c} else {
  
  allcomparisons <- rbind.data.frame(allcomparisons, c)
  
  }}

p <- ggplot(allcomparisons) +
  geom_histogram(aes(P), binwidth=0.1) +
  facet_wrap(~Type, nrow = 3) +
  xlab("P-Value")+
  ylab("Number of Simulated Datasets")+
  theme_bw()

ggsave(plot = p, file = "./figures/comphist_D0_L1_S0.png", 
       width = 6, height = 4, units = "in")


#############

ggplot(allcomparisons) +
  geom_boxplot(aes(x=Type, y = P)) +
  ylim(c(0,1))+
  xlab(NULL) +
  ylab("P-value")+
  theme_bw()

ggsave(plot = last_plot(), file = "./figures/compcheck_D0_L1_S0.png", 
       width = 6, height = 4, units = "in")

all_Dnonprob_lglk <- data.frame(matrix(rep(NA, 6*10), nrow = 10, ncol = 6))
names(all_Dnonprob_lglk) <- names(ptru)

for (i in 1:length(suffixes)){
  
  load(paste0("./results/Dnonprob_lglk_", suffixes[i], ".RData"))
  
    all_Dnonprob_lglk[i,] <- Dnonprob_lglk
    
}

all_Dnonprob_lglk %>%
  pivot_longer(names_to = "Type", values_to = "Value", cols = 1:6) %>%
  filter(Type != "RoI") %>%
  ggplot() +
  geom_boxplot(aes(x=Type, y = Value)) +
  theme_bw()


