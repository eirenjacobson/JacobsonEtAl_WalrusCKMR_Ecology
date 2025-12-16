library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

data <- read_excel("./results/Summary_D1_L1.xlsx")

longdata <- data %>% pivot_longer(cols = 2:7, names_to = "Parameter", values_to = "SE")

ggplot(longdata) +
  geom_line(aes(x=Years, y= SE, color = CKMR)) +
  facet_wrap(~Parameter, scales = "free", nrow = 6) +
  theme_bw()
