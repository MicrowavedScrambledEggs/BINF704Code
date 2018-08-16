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
  
  print(paste("pA:", newPA, "pB:", newPB, "pO:", newPO))
  
  if (newPA != pA || newPB != pB || newPO != pO) {
    pA <- newPA
    pB <- newPB
    pO <- newPO
    counter <- counter + 1
  } else {
    break
  }
}

nAEst <- (2* pA * pO + pA**2) * N
nBEst <- (2* pB * pO + pB**2) * N
nABEst <- 2 * pA * pB * N
nOEst <- pO**2 * N


