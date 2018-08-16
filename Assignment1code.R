nA <- 15
nB <- 10
nAB <- 3
nO <- 1 
N <- nA +nB + nAB + nO
realPhenFreqA <- nA / N 
realPhenFreqB <- nB / N
realPhenFreqAB <- nAB / N
realPhenFreqO <- nO / N

pA <- 0/5
pB <- 0/5
pO <- 5/5

newPA <- 0
newPB <- 0
newPO <- 0

counter <- 1

emBlood <- function(pA, pB, pO, nA, nB, nAB, nO) {
  N <- nA +nB + nAB + nO
  while(TRUE){
    # if(counter %% 2 == 1){
    #   nAA <- pA**2 * N
    #   nAO <- nA - nAA
    #   nBB <- pB**2 * N 
    #   nBO <- nB - nBB
    # } else {
    nAO <- 2* pA * pO * N
    nAA <- nA - nAO
    nBO <- 2* pB * pO * N
    nBB <- nB - nBO
    # }
    
    newPA <- (2*nAA + nAO + nAB) / (2*N)
    newPB <- (2*nBB + nBO + nAB) / (2*N)
    newPO <- (2*nO + nBO + nAO) / (2*N)
    
    # print(paste("pA:", newPA, "pB:", newPB, "pO:", newPO))
    
    if (newPA != pA || newPB != pB || newPO != pO) {
      pA <- newPA
      pB <- newPB
      pO <- newPO
    } else {
      break
    }
  }
  alleleFreq <- c(pA, pB, pO)
  names(alleleFreq) <- c("pA", "pB", "pO")
  return(alleleFreq)
}

for(i in 1:50){
  ps <- c(0,0,0)
  place <- sample(1:3, 3)
  ps[place[1]] <- runif(1, 0, 1)
  ps[place[2]] <- runif(1, 0, 1-ps[place[1]])
  ps[place[3]] <- 1-ps[place[1]]-ps[place[2]]
  print(paste("Initial probs: pA", ps[1], "pB:", ps[2], "pO:", ps[3]))
  print(emBlood(ps[1], ps[2], ps[3], nA, nB, nAB, nO))
}

nAEst <- (2* pA * pO + pA**2) * N
nBEst <- (2* pB * pO + pB**2) * N
nABEst <- 2 * pA * pB * N
nOEst <- pO**2 * N


