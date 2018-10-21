# Badi James (bjam575) BIOINF704 Assignment 4 Code

# Read in effect SNP table with summary statistic
effects <- read.table("Data/bmi-sig5e-3.txt", sep = "", header = TRUE, 
                      stringsAsFactors = F)

# Read in the genotype data for all the individuals
genFileNams <- paste0("Data/bmi/bmi.chr", 1:22, ".human.genotypes.txt")
chromGen <- sapply(genFileNams, function(x) 
  list(x=read.table(x, sep = "", header = TRUE, stringsAsFactors = F)))
chromGenMat <- chromGen[[1]]
for(i in 2:22)chromGenMat <- rbind(chromGenMat, chromGen[[i]])
# Read in the genotype data for the other homonids 
altai <- read.table("Data/bmi/bmi.Altai.genotypes.txt", sep = "", 
                    header = TRUE, stringsAsFactors = F)
denisovan <- read.table("Data/bmi/bmi.Denisovan.genotypes.txt", sep = "", 
                        header = TRUE, stringsAsFactors = F)
vindija <- read.table("Data/bmi/bmi.Vindija33.19.genotypes.txt", sep = "", 
                      header = TRUE, stringsAsFactors = F)

# Merge the tables together 
colnames(altai)[4:6] <- c("REF.Altai", "ALT.Altai", "Altai")
colnames(denisovan)[4:6] <- c("REF.Denisovan", "ALT.Denisovan", "Denisovan")
colnames(vindija)[4:6] <- c("REF.Vindija", "ALT.Vindija", "Vindija")
colnames(effects)[2] <- "CHROM"
colnames(chromGenMat)[3] <- "rsid"
chromGenMat <- merge(chromGenMat, altai, all = TRUE)
chromGenMat <- merge(chromGenMat, denisovan, all = TRUE)
chromGenMat <- merge(chromGenMat, vindija, all = TRUE)

superTable <- merge(effects, chromGenMat)

# Calculates the contribution of a loci (row in the table) to every 
# individual's PRS score
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

# Calculate score using data for all SNPs
score.mat <- apply(superTable, 1, calculateLociScore)
scores <- apply(score.mat,1,sum)
scores

# Normalisation

# Count the number of SNPs where the genotype data is NOT missing for an
# individual
calculate.m <- function(indiv){
  counts <- table(indiv)
  return(sum(counts[!names(counts) %in% c("./.")]))
}
m <- apply(superTable[,c(14:39, 42, 45,48)], 2, calculate.m)
# Normalise scores by number of SNPs for each individual
norm.scores <- scores / m
norm.scores

# Calculate score after LD clumping
clumpFileNams <- paste0("Data/bmi.clump/bmi.chr", 1:22, ".clumped")
clumpGen <- sapply(clumpFileNams, function(x) 
  list(x=read.table(x, sep = "", header = TRUE, stringsAsFactors = F)))
clumpGenMat <- clumpGen[[1]]
for(i in 2:22)clumpGenMat <- rbind(clumpGenMat, clumpGen[[i]])

clumpTab.5e3 <- superTable[superTable$rsid %in% clumpGenMat$SNP,]
clump.score.mat.5e3 <- apply(clumpTab.5e3, 1, calculateLociScore)
clump.score.5e3 <- apply(clump.score.mat.5e3,1,sum)
m.clump.5e3 <- apply(clumpTab.5e3[,c(14:39, 42, 45,48)], 2, calculate.m)
norm.clump.score.5e3 <- clump.score.5e3 / m.clump.5e3 

par(mar=c(5,4,4,2))
barplot(norm.scores, las=2, ylim = c(-1.5e-3, 1e-3), cex.names = 0.8,
        main = "Normalised PRS all SNP (p<5x10^-3)")
barplot(norm.clump.score.5e3, las=2, col = 2, cex.names = 0.8, ylim = c(0, 2.5e-3),
        main = "Normalised PRS LD clumped SNP (p<5x10^-3)")

# Filter by LD clumped SNPs by p-value then caculate PRS score
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
pvals <- c(1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7)
filter.scores <- sapply(pvals, function(p) {
  nam=as.character(p)
  l=list(clump.pvalue(superTable, clumpGenMat, p))
  names(l)=nam
  return(l)
})

# Output score table and plot distributions
q3TableNorm <- data.frame(
  names(norm.scores),norm.scores, norm.clump.score.5e3, 
  filter.scores$`0.001`$Norm.PRS.0.001, filter.scores$`5e-04`$`Norm.PRS.5e-04`, 
  filter.scores$`1e-04`$`Norm.PRS.1e-04`, filter.scores$`5e-05`$`Norm.PRS.5e-05`, 
  filter.scores$`1e-05`$`Norm.PRS.1e-05`, filter.scores$`5e-06`$`Norm.PRS.5e-06`, 
  filter.scores$`1e-06`$`Norm.PRS.1e-06`, filter.scores$`5e-07`$`Norm.PRS.5e-07`, 
  filter.scores$`1e-07`$`Norm.PRS.1e-07`
)
colnames(q3TableNorm) <- c(
  "ID","All SNP p<5e-3", "LD Clump p<5e-3", "LD Clump p<1e-3", "LD Clump p<5e-4",
  "LD Clump p<1e-4", "LD Clump p<5e-5", "LD Clump p<1e-5", "LD Clump p<5e-6",
  "LD Clump p<1e-6", "LD Clump p<5e-7", "LD Clump p<1e-7")
write.table(q3TableNorm, file = "BMI.txt", col.names = T, row.names = F, sep = "\t")

colnames(q3TableNorm) <- c("ID","All SNP p<5e-3", "p<5e-3", "p<1e-3", "p<5e-4",
                           "p<1e-4", "p<5e-5", "p<1e-5", "p<5e-6",
                           "p<1e-6", "p<5e-7", "p<1e-7")
boxplot(q3TableNorm[1:26,3:12], las =2, cex = 0.8)
stripchart(q3TableNorm[27,3:12], add = T, col = 2, vertical = T, pch=2)
stripchart(q3TableNorm[28,3:12], add = T, col = 6, vertical = T, pch=3)
stripchart(q3TableNorm[29,3:12], add = T, col = 4, vertical = T, pch =4)
legend("topleft", names(scores)[27:29], pch = c(2,3,4), col = c(2,6,4))
