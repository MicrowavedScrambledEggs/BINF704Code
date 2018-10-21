# Data
genotype <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 
              0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2)
phenotype <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
genPhen <- cbind(genotype, phenotype)

# Logistic regression model Likelihood
logistReg <- function(beta0, beta1, genPhen)
{
  # number representing genotype
  G <- c(0,1,2)
  # p for each genotype given the beta values
  p <- exp(beta0 + beta1 * G) / (1 + exp(beta0 + beta1 * G))
  names(p) <- G
  # number of cases and controls for each genotype
  nCase <- table(genPhen[genPhen[,2] == 1, 1])
  nCon <- table(genPhen[genPhen[,2] == 0, 1])
  # combined log likelihoods for cases of each genotype
  logLCase <- sapply(
    as.character(G), function(g) nCase[g]*log(p[g]))
  # combined log likelihoods for controls of each genotype
  logLCon <- sapply(
    as.character(G), function(g) nCon[g]*log(1 - p[g])) 
  # Remove NA values. Can happen if there are no cases or controls for a 
  # particular genotype
  logLCase <- logLCase[!is.na(logLCase)]
  logLCon <- logLCon[!is.na(logLCon)]
  
  # Sum all the log likelihoods
  sumLogL <- sum(logLCase) + sum(logLCon)
  return(sumLogL)
}

# log likelihood of beta0 = 0.1 and beta1 = 0.2 given all the data
logistReg(0.1, 0.2, genPhen)

# ML estimate of beta0 and beta1 using all the data
# starting with beta0 = 0.1 and beta1 = 0.2 just as arbitrary start points
# control = c(fnscale = -1) needed to make optim do maximization
mlAllData <- optim(c(0.1, 0.2), function(b) logistReg(b[1], b[2], genPhen), 
                   control = c(fnscale = -1))
mlAllData

# Try a couple different starting beta values to see if that effects 
# ML estimate
optim(c(0.5, 0.3), function(b) logistReg(b[1], b[2], genPhen), 
      control = c(fnscale = -1))
optim(c(0.8, 0.8), function(b) logistReg(b[1], b[2], genPhen), 
      control = c(fnscale = -1))
# Log likelihood values all idendical 
# parameter estimates all the same when rounding to 3 dp
# Looks good

# Generate bootstrap samples
# 100 samples seems like enough to measure variance???
B = 100

set.seed(123)
bootSamps <- list()
for(i in 1:B){
  bootSamps[[i]] <- genPhen[sample.int(nrow(genPhen), replace = TRUE),]
}

b0s <- c()
b1s <- c()
# Get the beta estimates for each sample
for(bootSamp in bootSamps){
  mlResult <- optim(
    c(0.1, 0.2), function(b) logistReg(b[1], b[2], bootSamp), 
    control = c(fnscale = -1))
  b0s <- c(b0s, mlResult$par[1])
  b1s <- c(b1s, mlResult$par[2])
}
# beta estimate means
avb0 <- mean(b0s)
avb1 <- mean(b1s)
# beta estimate variances
varsb0 <- (b0s - avb0)**2
varsb1 <- (b1s - avb1)**2
# estimated variance of the betas
estVarb0 <- sum(varsb0) / (length(b0s) - 1)
estVarb1 <- sum(varsb1) / (length(b1s) - 1)

# take a look at the distributions of the estimates
hist(b0s, xlab = "Beta0 ML estimates", 
     main = paste("Beta0 ML estimate distribution for", B, 
                  "bootstap samples"))
# Mark the beta0 estimate for all the data
abline(v = mlAllData$par[1], col = "red")
legend("topleft", legend = c("Beta0 ML estimate\nfor full dataset"), 
       col = c("red"), lty = c(1))

hist(b1s, xlab = "Beta1 ML estimates",
     main = paste("Beta1 ML estimate distribution for", B, 
                  "bootstap samples"))
# Mark the beta0 estimate for all the data
abline(v = mlAllData$par[2], col = "red")
legend("topleft", legend = c("Beta1 ML estimate\nfor full dataset"), 
       col = c("red"), lty = c(1))

# test for association
testAssoc <- avb1**2 / estVarb1
pAssoc <- pchisq(testAssoc, 1, lower.tail = FALSE)
# or
logBeta1.as.0 <- optim(
  0.1, function(b) logistReg(b, 0, bootSamp), method = "Brent",
  lower = -100, upper = 100,
  control = c(fnscale = -1))
testAssoc2 <- 2*abs(logBeta1.as.0$value - mlAllData$value)
pAssoc2 <- pchisq(testAssoc2, 1, lower.tail = FALSE)
