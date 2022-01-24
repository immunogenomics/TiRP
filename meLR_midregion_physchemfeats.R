source("utils.R")

args = commandArgs(trailingOnly=TRUE)
i = as.numeric(as.character(args[1]))
jobs = c("12_p108", "12_p109", "13_p108", "13_p109", "13_p110", "14_p108", "14_p109", "14_p110", "14_p112", "15_p108", "15_p109", "15_p110", "15_p111", "15_p112", "16_p108", "16_p109", "16_p110", "16_p111", "16_p112.1", "16_p112", "17_p108", "17_p109", "17_p110", "17_p111", "17_p111.1", "17_p112.1", "17_p112")

job = jobs[i]
len = strsplit(job, "_")[[1]][1] 
pos = strsplit(job, "_")[[1]][2]

rf = readRDS("data/aminoacid_physiochemical_feats.rds")
aa_features <- colnames(rf)[2:ncol(rf)]
aa_maps <- list()
for (i in 1:3){
  aa_maps[[i]] <- create_hashmap(as.character(rf$AA), as.numeric(as.character(rf[,(i+1)])))
}

data = readRDS("data/seay2016/seay_training.rds")
data = data[data$length==len,]

ind = which(colnames(data)==pos)
data$loc = as.character(data[,ind])
data$hydrophob = scale(sapply(data$loc, function(x) get_feat_score(x, aa_maps[[2]])))[,1]
data$pI = scale(sapply(data$loc, function(x) get_feat_score(x, aa_maps[[1]])))[,1]
data$volume = scale(sapply(data$loc, function(x) get_feat_score(x, aa_maps[[3]])))[,1]
data$cell_state = factor(data$cell_state, levels=c("Tconv", "Treg"))
fit = glmer(cell_state ~ hydrophob + pI + volume + (1|donor) + (1|site) + (1|donor:site), data=data, family="binomial")
stats = tidy(fit)

setwd("results")
saveRDS(stats, paste(job, "res.rds", sep="_"))
