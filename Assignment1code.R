nA <- 15
nB <- 10
nAB <- 3
nO <- 1 
N <- nA +nB + nAB + nO
realPhenFreqA <- nA / N 
realPhenFreqB <- nB / N
realPhenFreqAB <- nAB / N
realPhenFreqO <- nO / N

# Expectation Maximization for blood type allele frequencies
# Takes a vector of initial allele frequency guesses (pInit) and observed
# phenotype counts (nA, nB, nAB, nO)
emBlood <- function(pInit, nA, nB, nAB, nO) {
  pA <- pInit[1]
  pB <- pInit[2]
  pO <- pInit[3]
  N <- nA +nB + nAB + nO
  while(TRUE){
    nAO <- 2* pA * pO * N
    nAA <- nA - nAO
    nBO <- 2* pB * pO * N
    nBB <- nB - nBO
    
    newPA <- (2*nAA + nAO + nAB) / (2*N)
    newPB <- (2*nBB + nBO + nAB) / (2*N)
    newPO <- (2*nO + nBO + nAO) / (2*N)
    
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

# Create 50 random sets of initial allele frequencies
# Use to test if different inputs result in finding different local maxima
startfreq <- matrix(nrow = 50, ncol = 3)

for(i in 1:50){
  ps <- c(0,0,0)
  place <- sample(1:3, 3)
  bar <- runif(2, 0, 1)
  ps[place[1]] <- min(bar)
  ps[place[2]] <- max(bar) - min(bar)
  ps[place[3]] <- 1-max(bar)
  startfreq[i,] <- ps
}

hist(startfreq[,2])

maxLAlleleFreq <- apply(startfreq, 1, emBlood, nA = nA, nB = nB, nO = nO, nAB = nAB)
maxLAlleleFreq <- t(maxLAlleleFreq)
freqCounts <- list(table(maxLAlleleFreq[,1]), table(maxLAlleleFreq[,2]), table(maxLAlleleFreq[,3]))
names(freqCounts) <- colnames(maxLAlleleFreq)

pA <- maxLAlleleFreq[1,1]
pB <- maxLAlleleFreq[1,2]
pO <- maxLAlleleFreq[1,3]
nAEst <- (2* pA * pO + pA**2) * N
nBEst <- (2* pB * pO + pB**2) * N
nABEst <- 2 * pA * pB * N
nOEst <- pO**2 * N




