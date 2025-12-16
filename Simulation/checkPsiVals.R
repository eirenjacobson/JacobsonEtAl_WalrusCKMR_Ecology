

psi2 <- 0.1 # probability of breeding at 2-yr interval
psi3 <- 0.5 # probability of breeding at 3+ yr interval (+1st repro)
moveMat <- matrix(data = c(0, 1, 0, psi2, 0, 1-psi2, psi3, 0, 1-psi3),
                  nrow = 3, byrow = TRUE)

state <- rep(0, 100)
state[1] <- 3

for (i in 2:100){
  s <- state[i-1]
  state[i] <- sample(x = 1:3, size = 1, prob = moveMat[s,])
}


length(which(state == 2))



# Breeding psicle
# col= FROM, row= TO (Mark's version) so that T %*% oldprob = newprob
# Birth happens when Mum is in state 2
Tbreedf <- rbind( 
  c( 0, psi2, psi3),
  c( 1, 0, 0),
  c( 0, 1-psi2, 1-psi3)
)

# No longer going the eigenroute; get stationary eivec by...
M <- Tbreedf - diag( 3)
M[3,] <- 1
Pbreedf <- solve( M, c( 0,0,1))
recip_ppn_breedy <- 1 / Pbreedf[2] # of the popln, at any one time
