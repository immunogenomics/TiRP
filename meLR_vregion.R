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

covs = c("vgene", "p107")
if (length(unique(data$site))>1){
  covs = c(covs, "(1|site)")
}
if (length(unique(data$donor))>1){
  vgsr = readRDS("~/TiRP_github/results/VJgene_selectionrates/seay2016_VGSC.rds")
  colnames(vgsr)[2:3] = c("vgene", "VGSC")
  data = left_join(data, vgsr, by=c("donor", "vgene"))
  covs = c(covs, "VGSC", "(1|donor)", "(1|donor:site)")
}

obs = unique(as.character(data$vgene))
data$vgene = factor(data$vgene, levels=c("TCRBV05-01", obs[obs!="TCRBV05-01"]))
obs = unique(as.character(data$p107))
data$p107 = factor(data$p107, levels=c("S", obs[obs!="S"]))

f1 = as.formula(paste("cell_state~", paste0(covs, collapse = " + ")))

if (length(unique(data$site))>1){
fit1 = glmer(f1, data=data, family=binomial)
} else {
fit1 = glm(f1, data=data, family=binomial)
}
stats = tidy(fit1)
setwd("/PHShome/kl162/TiRP_github/results/TCRfeat_effectsizes/Vregion/seay2016")
saveRDS(stats, paste(tag, "vregion.rds", sep="_"))
