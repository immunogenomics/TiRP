

source("utils.R")

rf = readRDS("data/aminoacid_physiochemical_feats.rds")
aa_features <- colnames(rf)[2:ncol(rf)]
aa_maps <- list()
for (i in 1:3){
  aa_maps[[i]] <- create_hashmap(as.character(rf$AA), as.numeric(as.character(rf[,(i+1)])))
}

lengths = c("12", "13", "14", "15", "16", "17")
args = commandArgs(trailingOnly=TRUE)
i = as.numeric(as.character(args[1]))
len = lengths[i]

data = readRDS("data/seay2016/seay_training.rds")
data = data[!(is.na(data$cdr1)),]
data = data[data$length==len,]

poss = get_positions(seq_length=as.numeric(as.character(len)), CDR12 = TRUE)
poss = paste("p", poss, sep="")
feats = c("pI", "hydrophob", "volume")
x <- matrix(nrow=nrow(data), ncol=0)
for (i in 1:length(poss)){
  ind = which(colnames(data)==poss[i])
  for (j in 1:3){
    vals = scale(sapply(data[,ind], function(x) get_feat_score(x, aa_maps[[j]])))[,1]
    x = cbind(x, vals)
    colnames(x)[ncol(x)] = paste(poss[i], feats[j], sep="_")
  }
}

data$batch = paste(data$donor, data$site)
batches = sort(unique(as.character(data$batch)))
for (i in 2:length(batches)){
  vals = sapply(data$batch, function(x) ifelse(x==batches[i], 1, 0))
  x = cbind(x, vals)
  colnames(x)[ncol(x)] = paste("batch_is", batches[i])
}

y = ifelse(data$cell_state=="Treg", 1, 0)
fit = glmnet(x, y, data=data, alpha=0, lambda=0.01, family="binomial")

stats = tidy(fit)

setwd("results/TCRfeat_effectsizes/physicochemical_features/ridge_regularized/seay2016")
saveRDS(stats, paste(len, "glmnet.rds", sep="_"))
