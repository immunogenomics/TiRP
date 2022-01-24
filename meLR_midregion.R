library(plyr)
library(dplyr)
library(lme4)
library(broom.mixed)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
i = as.numeric(as.character(args[1]))
data = readRDS("data/seay2016/seay_training.rds")
data$cell_state = factor(data$cell_state, levels=c("Tconv", "Treg"))

data$donor = as.character(data$donor)
drs = unique(data$donor)
jobs = c("seaytraining", "seayCASE", "seayCONTROL", drs)
job = jobs[i]

if (job=="seayCASE"){
  data = data[data$pheno=="T1D",]
} else if (job=="seayCONTROL"){
  data[data$pheno=="CONTROL",]
} else if (job!="seaytraining"){
  data = data[data$donor == job,]
}

covs = c(paste("perc_mid", c("V", "C", "D", "E", "F", "H", "I", "L", "N", "Q", "R", "S", "T", "W", "Y"), sep="_"), "length")
if (length(unique(data$site))>1){
  covs = c(covs, "(1|site)")
}
if (length(unique(data$donor))>1){
  covs = c(covs, "(1|donor)", "(1|donor:site)")
}

data$length = as.character(data$length)
data$length = factor(data$length, levels=c("15", "12", "13", "14", "16", "17"))
ncovs = covs[!(covs %in% c("(1|donor:site)", "length", "(1|site)", "(1|donor)"))]
if (length(ncovs)>0){
  for (j in 1:length(ncovs)){
    ind = which(colnames(data)==ncovs[j])
    data[,ind] = scale(data[,ind])[,1]
  }
}

f1 = as.formula(paste("cell_state~", paste0(covs, collapse = " + ")))

if (length(unique(data$site))>1){
fit1 = glmer(f1, data=data, family=binomial)
} else {
fit1 = glm(f1, data=data, family=binomial)
}
stats = tidy(fit1)
setwd("results/TCRfeat_effectsizes/midregion_percentAA")
saveRDS(stats, paste(job, "15prcAAs.rds", sep="_"))

