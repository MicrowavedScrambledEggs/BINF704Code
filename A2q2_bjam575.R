library(genetics)

fms.Dat <- read.table("FMS_data.txt")

# a) Testing for hardy weiberg equilibrium

snp.to.test <- genotype(fms.Dat$akt1_t10726c_t12868c, sep = 1)
hweTest <- HWE.test(snp.to.test)
hweTest

# Stratifying by Race
races <- levels(fms.Dat$Race)
raceTests <- vector("list", length = length(races))
names(raceTests) <- races
for(race in races){
  raceGen <- fms.Dat$akt1_t10726c_t12868c[which(fms.Dat$Race == race)]
  if(length(raceGen) >= 30){
    snp.to.test.race <- genotype(raceGen, sep = 1)
    hweTest.race <- HWE.test(snp.to.test.race)
    raceTests[[race]] <- hweTest.race
  }
}
raceTests

# # looking for Wahlund effect
# mi <- c()
# pAi <- c()
# pai <- c()
# races <- levels(fms.Dat$Race)
# races <- c(races, NA)
# for(race in races){
#   if(is.na(race))
#     raceGen <- fms.Dat$akt1_t10726c_t12868c[which(is.na(fms.Dat$Race))]
#   else
#     raceGen <- fms.Dat$akt1_t10726c_t12868c[which(fms.Dat$Race == race)]
#   if(length(raceGen) >= 30){
#     raceGen <- raceGen[which(!is.na(raceGen))]
#     mi <- c(mi, length(raceGen) / nrow(fms.Dat))
#     pAi <- c(pAi, (table(raceGen)[1]*2 + table(raceGen)[2]) / (2*length(raceGen)))
#     pai <- c(pai, (table(raceGen)[3]*2 + table(raceGen)[2]) / (2*length(raceGen)))
#   }
# }
# EpA <- sum(mi*pAi)
# pAA <- EpA**2 + var(pAi)
# Epa <- sum(mi*pai)
# paa <- Epa**2 + var(pAi)
# pAa <- 1 - pAA - paa
# # Expected gen counts
# length(which(!is.na(fms.Dat$akt1_t10726c_t12868c)))*c(pAA, pAa, paa)
# # compare to actual
# table(fms.Dat$akt1_t10726c_t12868c)

# b) testing association between genotypes and phenotype

individuals <- which(!is.na(fms.Dat$esr1_rs1042717) & !is.na(fms.Dat$pre.BMI))
loci.to.test <- fms.Dat$esr1_rs1042717[individuals]
phen.to.test <- rep(0, length(individuals))
phen.to.test[fms.Dat$pre.BMI[individuals] > 25] <- 1
phen.to.test <- factor(phen.to.test)

association.Test <- chisq.test(loci.to.test, phen.to.test)
observed <- association.Test$observed
armitage.test <- prop.trend.test(observed[,2], apply(observed, 1, sum),
                                 score = c(2,1,0))

# stratifying by race
race.Assoc.Tests <- vector("list", length = length(races))
names(race.Assoc.Tests) <- races
race.Trend.Tests <- vector("list", length = length(races))
names(race.Trend.Tests) <- races
for(race in races){
  raceRows <- which(fms.Dat$Race[individuals] == race)
  raceGen <- fms.Dat$esr1_rs1042717[raceRows]
  racePhen <- rep(0, length(raceGen))
  racePhen[fms.Dat$pre.BMI[raceRows] > 25] <- 1
  racePhen <- factor(racePhen)
  if(length(raceGen) >= 30){
    observed <- table(raceGen, racePhen)
    if(TRUE %in% (5 > observed))
      assoc.test <- fisher.test(observed)
    else
      assoc.test <- chisq.test(observed)
    trend.test <- prop.trend.test(observed[,2], apply(observed, 1, sum),
                                  score = c(2,1,0))
    race.Assoc.Tests[[race]] <- assoc.test
    race.Assoc.Tests[[race]]$observed <- observed
    race.Trend.Tests[[race]] <- trend.test
  }
}
race.Assoc.Tests

