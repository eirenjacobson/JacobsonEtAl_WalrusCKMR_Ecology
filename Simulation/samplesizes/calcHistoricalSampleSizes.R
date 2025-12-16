
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

load("./simulation/WalrusSim.RData")

df <- read_excel("./samplesizes/Sample_Summaries_byAgeSex.xlsx") %>%
  pivot_longer(cols = c("A.2013", "A.2014", "A.2015", "A.2016", "A.2017"), 
               names_to = "CaptureYear", values_to = "Number") %>%
  mutate(CR = as.numeric(str_split_i(CaptureYear, "A.", i = 2))) %>%
  filter(Year == CR) %>%
  select(Age, Sex, Year, Number) %>%
  pivot_wider(names_from = Age, values_from = Number)

save(df, file = "./samplesizes/historicalsamplesummary.RData")

calves <- data.frame("Region" = "USA", "Age" = 0, "Sex" = df$Sex, "Year" = df$Year, "Number" = df$`0`)
oners <- data.frame("Region" = "USA", "Age" = 1, "Sex" = df$Sex, "Year" = df$Year, "Number" = df$`1`)
twoers <- data.frame("Region" = "USA", "Age" = 2, "Sex" = df$Sex, "Year" = df$Year, "Number" = df$`2`)
threers <- data.frame("Region" = "USA", "Age" = 3, "Sex" = df$Sex, "Year" = df$Year, "Number" = df$`3`)

# find proportions of individuals of different ages sampled in each year

# four and five year olds

ff <- subset(samples, AgeLast == 4 | AgeLast == 5 )
ff.sampbyyear <- expand.grid("Year" = 2013:2017, "Sex" = c("M", "F"), "Age" = c(4, 5))
ff.sampbyyear$Number <- NA

for (i in 1:nrow(ff.sampbyyear)){
  
  sim.s <- filter(indiv, AgeLast == 4 | AgeLast == 5 , is.na(DeathY)) |> 
#   filter(SampY == ff.sampbyyear$Year[i]) %>%
    count(AgeLast)
  
  hist.s <- df$`4to5`[df$Year == ff.sampbyyear$Year[i] & df$Sex == ff.sampbyyear$Sex[i]]
  
  p <- sim.s$n[sim.s$AgeLast == ff.sampbyyear$Age[i]]/sum(sim.s$n)
  
  ff.sampbyyear$Number[i] <- round(p*hist.s)

}

fourfive <- data.frame("Region" = "USA", 
                       "Age" = ff.sampbyyear$Age, 
                       "Sex" = ff.sampbyyear$Sex,
                       "Year" = ff.sampbyyear$Year, 
                       "Number" = ff.sampbyyear$Number)

# six plus

sp <- subset(samples, AgeLast %in% 6:44)
sp.sampbyyear <- expand.grid("Year" = 2013:2017, "Sex" = c("M", "F"), "Age" = 6:44)
sp.sampbyyear$Number <- NA

for (i in 1:nrow(sp.sampbyyear)){
  
  sim.s <- filter(indiv, AgeLast %in% 6:44, is.na(DeathY)) |> 
 #   filter(SampY == sp.sampbyyear$Year[i]) %>%
    count(AgeLast)
  
  hist.s <- df$A[df$Year == sp.sampbyyear$Year[i] & df$Sex == sp.sampbyyear$Sex[i]]
  
  p <- sim.s$n[sim.s$AgeLast == sp.sampbyyear$Age[i]]/sum(sim.s$n)
  
  sp.sampbyyear$Number[i] <- round(p*hist.s)
  
}

sixplus <- data.frame("Region" = "USA", 
                      "Age" = sp.sampbyyear$Age, 
                      "Sex" = sp.sampbyyear$Sex,
                      "Year" = sp.sampbyyear$Year, 
                      "Number" = sp.sampbyyear$Number)

historical <- rbind.data.frame(calves, oners, twoers, threers, fourfive, sixplus)

save(historical, file = "historicalsamplesizes.RData")


