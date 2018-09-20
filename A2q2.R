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

# b) testing association between genotypes and phenotype

loci.to.test <- fms.Dat$esr1_rs1042717
phen.to.test <- rep(0, nrow(fms.Dat))
phen.to.test[fms.Dat$pre.BMI > 25] <- 1
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
  raceRows <- which(fms.Dat$Race == race)
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
    race.Trend.Tests[[race]] <- trend.test
  }
}
race.Assoc.Tests

