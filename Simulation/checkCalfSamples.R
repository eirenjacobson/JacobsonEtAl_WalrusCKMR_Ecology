# check whether calves were sampled with mother

calfsamples <- filter(allsamples, AgeLast == 0)
calfsamples$MumToo <- NA

for (i in 1:nrow(calfsamples)){

  m <- calfsamples$Mum[i]
  calfsamples$MumToo[i] <- any(allsamples$Me == m & allsamples$SampY == calfsamples$SampY[i])
  
}

sum(calfsamples$MumToo)/nrow(calfsamples)
