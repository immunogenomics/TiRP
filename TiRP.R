
path_to_data = "/PHShome/kl162/TiRP_data/data"
path_to_results = "/PHShome/kl162/TiRP_github/TiRP/results"

reformat_vgene <- function(vg){
  vg = as.character(vg)
  vgene = gsub("TRBV", "TCRBV", vg)
  info = strsplit(vgene, "-")[[1]]
  family = info[1]
  member = ifelse(length(info)==2, info[2], "01")
  family = ifelse(nchar(family)==6, paste("TCRBV0", substr(family, 6, 6), sep=""), family)
  member = ifelse(substr(member, 1, 1)=="0", member, paste("0", member, sep=""))
  return(paste(family, member, sep="-"))
}

get_perc_score <- function(percents, touse, perc_term){
  df = data.frame(perc_term, percents)
  df = df[df$perc_term %in% touse$feat,]
  df = left_join(df, touse, by=c("perc_term"= "feat"))
  df$prod = df$percents * df$metabeta
  return(sum(df$prod))
}

get_pos_score <- function(x, touse, pos_terms){
  terms = pos_terms[[nchar(x)-1]]
  terms = sapply(1:length(terms), function(y) paste(terms[y], substr(x, y, y), sep="_"))
  return(sum(touse$metabeta[touse$feat %in% terms]))
}

##data.frame with first column corresponding to v_gene, second column corresponding to CDR3 sequence
TiRP <- function(data){
  weights = readRDS(paste(path_to_results, "TiRP_scoring/TiRP_weights.rds", sep="/"))

  suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
  })

  ######################## defining TCR features ###########################

  ##standardizing TRBV gene names
  ex = data[,1][!(grepl("nresolved", data[,1]))][1]
  #if (substr(data[1,1], 1, 4)=="TRBV"){
  if (!(grepl("-0", ex))){
     data$vgene = sapply(data[,1], function(x) reformat_vgene(x))
  } else {
    data$vgene = as.character(data[,1])
  }

  ##CDR3 sequence
  data$cdr3 = as.character(data[,2])
  data = data[!(is.na(data$cdr3)),]

  ##CDR3 length
  data$length = sapply(data$cdr3, function(x) nchar(x))
  data = data[data$length >=12 & data$length<=17,]

  ##Jmotif
  jmotifs = gsub("Jmotif", "", weights$feat[grepl("Jmotif", weights$feat)])
  data$Jmotif = sapply(data$cdr3, function(x) substr(x, nchar(x)-4, nchar(x)))
  data$Jmotif[!(data$Jmotif %in% jmotifs)] = "other"

  ##amino acid composition in the CDR3 middle region
  data$cdr3MR = sapply(data$cdr3, function(x) substr(x, 5, nchar(x)-6))
  aminos = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  ref_stats = readRDS(paste(path_to_results, "TCR_structure/heldout_means_sds.rds", sep="/"))
  for (i in 1:length(aminos)){
    new = sapply(data$cdr3MR, function(x) str_count(x, aminos[i])/nchar(x))
    mn = ref_stats$mean[ref_stats$amino==aminos[i]]
    sd = ref_stats$sd[ref_stats$amino==aminos[i]]
    new = (new - mn)/sd
    data = cbind(data, new)
    colnames(data)[ncol(data)] = paste("perc_mid", aminos[i], sep="_")
  }


  ######################## scoring TCR features ###########################

  ##including only terms that were Bonferroni-significant in meta-analysis:
  touse = weights[weights$metap<(0.05/nrow(weights)),]

  perc_terms <- touse$feat[grepl("perc_mid", touse$feat)]
  data$perc_score = sapply(1:nrow(data), function(x) sum(data[x, perc_terms] * touse[perc_terms, c("metabeta")]))

  pos_terms <- list()
  pos_terms[[1]] <- c("12_5", "12_6")
  pos_terms[[2]] <- c("13_5", "13_6", "13_7")
  pos_terms[[3]] <- c("14_5", "14_6", "14_7", "14_-7")
  pos_terms[[4]] <- c("15_5", "15_6", "15_7", "15_8", "15_-7")
  pos_terms[[5]] <- c("16_5", "16_6", "16_7", "16_8", "16_-8", "16_-7")
  pos_terms[[6]] <- c("17_5", "17_6", "17_7", "17_8", "17_9", "17_-8", "17_-7")
  data$pos_score = sapply(data$cdr3MR, function(x) get_pos_score(x, touse, pos_terms))

  data$feat = paste("vgene", data$vgene, sep="")
  data = left_join(data, touse[,c("feat", "metabeta")], by="feat")
  colnames(data)[ncol(data)] = "vgene_score"
  data$vgene_score = sapply(data$vgene_score, function(x) ifelse(is.na(x), 0, x))
  if (nrow(data[data$vgene_score!=0,])==0){
    stop("unidentifiable TRBV genes. match TRBV gene formatting in TiRP_demo.ipynb or TiRP_weights.rds")
  }

  data$feat = paste("p107", substr(data$cdr3, 4, 4), sep="")
  data = left_join(data, touse[,c("feat", "metabeta")], by="feat")
  colnames(data)[ncol(data)] = "p107_score"
  data$p107_score = sapply(data$p107_score, function(x) ifelse(is.na(x), 0, x))

  data$feat = paste("Jmotif", data$Jmotif, sep="")
  data = left_join(data, touse[,c("feat", "metabeta")], by="feat")
  colnames(data)[ncol(data)] = "Jmotif_score"
  data$Jmotif_score = sapply(data$Jmotif_score, function(x) ifelse(is.na(x), 0, x))

  data$feat = sapply(data$cdr3, function(x) substr(x, nchar(x)-5, nchar(x)-5))
  data$feat = paste("p113", data$feat, sep="")
  data = left_join(data, touse[,c("feat", "metabeta")], by="feat")
  colnames(data)[ncol(data)] = "p113_score"
  data$'p113_score' = sapply(data$'p113_score', function(x) ifelse(is.na(x), 0, x))

  data$feat = paste("length", data$length, sep="")
  data = left_join(data, touse[,c("feat", "metabeta")], by="feat")
  colnames(data)[ncol(data)] = "length_score"
  data$length_score = sapply(data$length_score, function(x) ifelse(is.na(x), 0, x))
  data$feat <- NULL

  ######################## summation ###########################
  data$total_score = data$vgene_score + data$Jmotif_score + data$p107_score + data$p113_score + data$pos_score + data$perc_score + data$length_score
  ##scaling by the mean and standard deviation of originally held-out data to standardize the TiRP scale
  data$vTiRP = (data$vgene_score + data$p107_score + 0.1459054)/0.2364
  data$mTiRP = (data$perc_score + data$pos_score + data$length_score - 0.03846178)/0.2364
  data$jTiRP = (data$Jmotif_score + data$p113_score - 0.07454362)/0.2364
  data$TiRP = (data$total_score + 0.0329)/0.2364
  data = data[,!(grepl("perc_mid", colnames(data)))]
  return(data)
}
