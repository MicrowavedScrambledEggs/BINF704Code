calculateLociScore <- function(loci, impute=FALSE){
  # Deciding how to deal with missing genotypes
  # Either impute using the provided effect allele frequency and HWE
  # Or treat as 0
  alFrq <- loci[["EAF_HapmapCEU"]]
  missGen <- if(impute) 2*alFrq**2 + 2*alFrq*(1-alFrq) else 0 
  # Assign genotype numbers to genotypes depending on if the Ref or the Alt 
  # allele is the effect allele
  refEffect <- c("0/0"=2, "0/1"=1, "1/1"=0, "./."=missGen, "1/2"=0)
  altEffect <- c("0/0"=0, "0/1"=1, "1/1"=2, "./."=missGen, "1/2"=1)
  g <- if(!is.na(loci[["REF"]]) & loci[["REF"]] == loci[["A1"]]) 
    refEffect else altEffect
  g.altai <- if(!is.na(loci[["REF.Altai"]]) & 
                loci[["REF.Altai"]] == loci[["A1"]]) 
    refEffect else altEffect
  g.denisovan <- if(!is.na(loci[["REF.Denisovan"]]) & 
                     loci[["REF.Denisovan"]] == loci[["A1"]]) 
    refEffect else altEffect
  g.vindija <- if(!is.na(loci[["REF.Vindija"]]) & 
                  loci[["REF.Vindija"]] == loci[["A1"]]) 
    refEffect else altEffect
  # Treat NA genotypes as missing genotypes
  gens <- loci[c(14:39,42,45,48)]
  gens[which(is.na(gens))] <- "./."
  loci[c(14:39,42,45,48)] <- gens
  # Calculate loci's contribution to the PRS for all individuals
  b <- as.numeric(loci[["b"]])
  score.hs <- sapply(loci[14:39], function(x) g[x]*b)
  score.a <- g.altai[loci[["Altai"]]]*b
  score.d <- g.denisovan[loci[["Denisovan"]]]*b
  score.v <- g.vindija[loci[["Vindija"]]]*b
  scores <- c(score.hs, score.a, score.d, score.v)
  names(scores) <- c(names(loci[14:39]), "Altai", "Denisovan", "Vindija")
  return(scores)
}

calculate.m <- function(indiv){
  counts <- table(indiv)
  return(sum(counts[!names(counts) %in% c("./.")]))
}

clump.pvalue <- function(superTable, clumpTable, pvalue){
  filterClump <- clumpTable[clumpTable$P < pvalue,]
  filterSuper <- superTable[superTable$rsid %in% filterClump$SNP, ]
  
  clump.score.mat <- apply(filterSuper, 1, calculateLociScore)
  clump.score <- apply(clump.score.mat,1,sum)
  m.clump <- apply(filterSuper[,c(14:39, 42, 45,48)], 2, calculate.m)
  norm.clump.score <- clump.score / m.clump
  
  toReturn <- list(clump.score, m.clump, norm.clump.score)
  names(toReturn) <- paste0(c("PRS.", "#.SNP.", "Norm.PRS."), as.character(pvalue))
  return(toReturn)
}
