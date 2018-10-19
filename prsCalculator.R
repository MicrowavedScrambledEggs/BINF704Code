effects <- read.table("Data/bmi-sig5e-3.txt", sep = "", header = TRUE)
genFileNams <- paste0("Data/bmi/bmi.chr", 1:22, ".human.genotypes.txt")

chromGen <- sapply(genFileNams, function(x) list(x=read.table(x, sep = "", header = TRUE)))
chromGenMat <- chromGen[[1]]
for(i in 2:22)chromGenMat <- rbind(chromGenMat, chromGen[[i]])
colnames(effects)[1] <- "ID"
colnames(chromGenMat)[1] <- "CHR"
superTable <- merge(effects, chromGenMat)

calculateLociScore <- function(loci, impute=FALSE){
  alFrq <- loci$EAF_HapmapCEU
  missGen <- if(impute) 2*alFrq**2 + 2*alFrq*(1-alFrq) else 0 
  refEffect <- c("0/0"=2, "0/1"=1, "1/1"=0, "./."=missGen)
  altEffect <- c("0/0"=0, "0/1"=1, "1/1"=2, "./."=missGen)
  g <- if(loci$REF == loci$A1) refEffect else altEffect
  
}
