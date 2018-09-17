# Returns the log-likelihood for the logistic model given genotype counts
logisticlk.genotypes <- function(betas, nCaseAA = 6, nCaseAa = 3, nCaseaa= 10, nCtrlAA = 1, nCtrlAa = 5, nCtrlaa = 6)
{
  # Parameters of the logistic regression model
  b0 = betas[1];
  b1 = betas[2];
  
  
  # Calculate the probabilities of the six possible outcomes 
  # note that exp(z) / (1 + exp(z)) is equal to 1 / (exp(-z) + 1) ( by multiplying top and bottom by exp(-z) )
  pCaseAA = 1 / ( 1 + exp( -b0 - 0 * b1 ) );
  pCtrlAA = 1 - pCaseAA;
  pCaseAa = 1 / ( 1 + exp( -b0 - 1 * b1) );
  pCtrlAa = 1 - pCaseAa;
  pCaseaa = 1/(1 + exp( -b0 - 2 * b1));
  pCtrlaa = 1 - pCaseaa;
  
  
  # Calculate the log-likelihood
  loglkval =
    nCaseAA * log(pCaseAA) +
    nCaseAa * log(pCaseAa) +
    nCaseaa * log(pCaseaa) +
    nCtrlAA * log(pCtrlAA) +
    nCtrlAa * log(pCtrlAa) +
    nCtrlaa * log(pCtrlaa) ;
  
  return(-loglkval);
}

## use optim function to find the maximum.  See ? optim

