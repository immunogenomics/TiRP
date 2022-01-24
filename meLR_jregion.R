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
jobs = c("all", "CASE", "CONTROL", drs)
job = jobs[i]
tag = paste(job, i, sep="_")
if (job=="CASE"){
  data = data[data$pheno=="T1D",]
} else if (job=="CONTROL"){
  data = [data$pheno=="CONTROL",]
} else if (job!="all"){
  data = data[data$donor == job,]
}

covs = c("Jmotif", "p113")
if (length(unique(data$site))>1){
  covs = c(covs, "(1|site)")
}

obs = unique(as.character(data$Jmotif))
data$Jmotif = factor(data$Jmotif, levels=c("TEAFF", obs[obs!="TEAFF"]))
obs = unique(as.character(data$p113))
data$p113 = factor(data$p113, levels=c("N", obs[obs!="N"]))

f1 = as.formula(paste("cell_state~", paste0(covs, collapse = " + ")))

if (length(unique(data$site))>1){
fit1 = glmer(f1, data=data, family=binomial)
} else {
fit1 = glm(f1, data=data, family=binomial)
}
stats = tidy(fit1)
setwd("/PHShome/kl162/TiRP_github/results/TCRfeat_effectsizes/Jregion/seay2016")
saveRDS(stats, paste(tag, "jregion.rds", sep="_"))
