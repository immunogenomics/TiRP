add_loc <- function(var){
  ind = which(colnames(data)==var)
  data$loc = as.character(data[,ind])
  obs = unique(data$loc)
  if (var!="p108"){
    data$loc = factor(data$loc, levels=c("G", obs[obs!="G"]))
  } else {
    data$loc = factor(data$loc, levels=c("L", obs[obs!="L"]))
  }
  return(data)
}

library(dplyr)
library(lme4)
library(broom.mixed)

source("utils.R")
data = readRDS("data/seay_training.rds")
lens = c(12, 13, 14, 15, 16, 17)
poss = vector(mode="character")
for (l in 1:length(lens)){
  poss = c(poss, paste(lens[l], get_positions(lens[l], mid=TRUE), sep=".p"))
}
feats = paste(poss, "aa", sep=".")

args = commandArgs(trailingOnly=TRUE)

drs = c("seay70")
job_grid = expand.grid(feats, drs)
jobs = paste(job_grid$Var1, job_grid$Var2)

i = as.numeric(as.character(args[1]))
job = jobs[i]
jinfo = strsplit(job, " ")[[1]]
feat = jinfo[1]
dr = jinfo[2]

if (dr != "seay70"){
  data = data[data$donor == dr,]
}

info = strsplit(feat, "\\.")[[1]]
var = info[2]
len = info[1]
ft = info[3]

if (len!="all"){
  data = data[data$length==len,]
}

data <- data
covs = c("loc", "(1|donor)","(1|site)", "(1|donor:site)")
data = add_loc(var)

covs = c(covs, "perc_mid_C", "perc_mid_H", "perc_mid_E", "perc_mid_D", "perc_mid_R", "perc_mid_W", "perc_mid_Y", "perc_mid_F", "perc_mid_L", "perc_mid_I", "perc_mid_T", "perc_mid_N", "perc_mid_S", "perc_mid_Q")

tag = paste0(c(dr, len, var, ft), collapse="_")
data$length = as.character(data$length)
data$length = factor(data$length, levels=c("15", "12", "13", "14", "16", "17"))
ncovs = covs[!(covs %in% c("(1|donor:site)", "length", "loc", "(1|site)", "(1|donor)"))]
if (length(ncovs)>0){
for (j in 1:length(ncovs)){
  ind = which(colnames(data)==ncovs[j])
  print(ncovs[j])
  print(ind)
  data[,ind] = scale(data[,ind])[,1]
}
}
data$cell_state = factor(data$cell_state, levels=c("Tconv", "Treg"))

setwd("results/repr_seay70_posbetas")
if (dr!="seay70"){
  covs = covs[!(grepl("(1|donor", covs))]
} 
f1 = as.formula(paste("cell_state~", paste0(covs, collapse = " + ")))

if (length(unique(data$site))>1){
  fit1 = glmer(f1, data=data, family=binomial)
} else {
  covs = covs[covs!="(1|site)"]
  fit1 = glm(f1, data=data, family=binomial)
}
print(f1)
print(getwd())
stats = tidy(fit1)
covs = covs[covs!="loc"]

f0 = as.formula(paste("cell_state~", paste0(covs, collapse = " + ")))
if (length(unique(data$site))>1){
  fit0 = glmer(f0, data=data, family=binomial)
} else {
  covs = covs[covs!="(1|site)"]
  fit0 = glm(f0, data=data, family=binomial)
}
print(f0)
print(getwd())
an = anova(fit1, fit0, test="Chisq")
t = table(data$loc)
save(stats, an, t, file=paste(tag, "res.RData", sep=""))
