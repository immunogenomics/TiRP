path_to_data = "/PHShome/kl162/TiRP_data/data"
path_to_results = "/PHShome/kl162/TiRP_github/TiRP/results"

suppressPackageStartupMessages({
library(stringr)
library(DescTools)
library(pheatmap)
library(lme4)
library(broom.mixed)
library(ggseqlogo)
library(VIF)
library(glmnet)
library(harmony)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(dplyr)
library(infotheo)
library(pals)
library(hash)
library(viridis)
library(Morpho)
library(class)
})

source("TiRP.R")
source("https://raw.githubusercontent.com/dorianps/corrplot2/master/corrplot2.R")
options(warn=-1)

reformat_vgene <- function(vg){
  vgene = gsub("TRBV", "TCRBV", vg)
  info = strsplit(vgene, "-")[[1]]
  family = info[1]
  member = ifelse(length(info)==2, info[2], "01")
  family = ifelse(nchar(family)==6, paste("TCRBV0", substr(x, 6, 6), sep=""), family)
  member = ifelse(substr(member, 1, 1)=="0", member, paste("0", member, sep=""))
  return(paste(family, member, sep="-"))
}

aminos <- function(){
  return(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
}

amino_color_scheme <- function(){
  return(make_col_scheme(chars=aminos(), groups=c('neutral', 'hydrophobic', 'acidic', 'acidic', 'very hydrophobic', 'neutral', 'basic', 'hydrophobic', 'basic', 'very hydrophobic', 'hydrophobic', 'neutral', 'neutral', 'neutral', 'basic', 'neutral', 'neutral', 'neutral', 'very hydrophobic', 'very hydrophobic'),

                  cols=c('black', 'darkgoldenrod', 'navyblue', 'navyblue', 'darkgoldenrod2', 'black', 'coral2', 'darkgoldenrod', 'coral2', 'darkgoldenrod2', 'darkgoldenrod', 'black', 'black', 'black', 'coral2', 'black', 'black', 'black', 'darkgoldenrod2', 'darkgoldenrod2')))

}

get_positions <- function(seq_length, mid=FALSE, CDR12 = FALSE){
  poss = c("104", "105", "106", "107", "108", "109", "110", "111", "111.1", "112.1", "112", "113", "114", "115", "116", "117", "118")
  if (mid){
    poss = poss[!(poss %in% c("104", "105", "106", "107", "113", "114", "115", "116", "117", "118"))]
  }
  if (seq_length==12){
    poss = poss[!(poss %in% c("110", "111", "111.1", "112.1", "112"))]
  }
  if (seq_length==13){
    poss = poss[!(poss %in% c("111", "111.1", "112.1","112"))]
  }
  if (seq_length==14){
    poss = poss[!(poss %in% c("111.1", "112.1", "111"))]
  }
  if (seq_length==15){
    poss = poss[!(poss %in% c("111.1", "112.1"))]
  }
  if (seq_length==16){
    poss = poss[!(poss %in% c("111.1"))]
  }
  if (CDR12){
    poss = c("27", "28", "29", "37", "38", "56", "57", "58", "63", "64", "65", poss)
  }
  return(poss)
}

mutinf_heatmap <- function(seq_length, add_conditional_MI=FALSE, recomputeMI = FALSE){
  if (seq_length!="all"){
    poss = paste("p", get_positions(seq_length), sep="")
    sub = data[data$length==seq_length,]
  } else {
    poss = c("vgene", "p27", "p28", "p29", "p37", "p38", "p56", "p57", "p58", "p63", "p64", "p65", "Vmotif", "p107", "p113", "Jmotif", "jgene")
    sub = data[!(is.na(data$cdr1)),]
  }
  if (add_conditional_MI){
    mat = get_condMImatrix(poss, seq_length, recomputeMI)
  } else {
    mat = matrix(nrow=length(poss), ncol=length(poss))
  }
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
        if (!(add_conditional_MI) | j<=i){
          indi = which(colnames(sub)==poss[i])
          indj = which(colnames(sub)==poss[j])
          mat[i,j] = MutInf(sub[,indi], sub[,indj])/mean(c(Entropy(table(sub[,indi])), Entropy(table(sub[,indj]))))
        }
      }
    }
  rownames(mat) = poss
  colnames(mat) = poss
  return(pheatmap(mat, cluster_rows = FALSE, cluster_cols=FALSE))
}

add_positionalAA <- function(exc_VJ=FALSE, inc_CDR12 = FALSE){
  if ("cdr3" %in% colnames(data)){
    data$sequence = as.character(data$cdr3)
  }
  data$length <- sapply(data$sequence, function(x) nchar(x))
  data <- data[data$length>=12 & data$length<=17,]
  data$Vmotif = sapply(data$sequence, function(x) substr(x, 1, 3))
  poss = get_positions(17)
  for (i in 1:6){
    if (!(exc_VJ & i<4)){
      data$new <- sapply(data$sequence, function(x) substr(x, i, i))
      colnames(data)[ncol(data)] <- paste("p", poss[i], sep="")
    }
  }
  data$p110 <- sapply(data$sequence, function(x) ifelse(nchar(x)>12, substr(x, 7, 7), "*"))
  data$p111 <- sapply(data$sequence, function(x) ifelse(nchar(x)>14, substr(x, 8, 8), "*"))
  data$'p111.1' <- sapply(data$sequence, function(x) ifelse(nchar(x)>16, substr(x, 9, 9), "*"))
  data$'p112.1' <- sapply(data$sequence, function(x) ifelse(nchar(x)>15, substr(x, nchar(x)-7, nchar(x)-7), "*"))
  data$'p112' <- sapply(data$sequence, function(x) ifelse(nchar(x)>13, substr(x, nchar(x)-6, nchar(x)-6), "*"))
  for (i in 5:0){
    if (!(exc_VJ & i<5)){
      data$new <- sapply(data$sequence, function(x) substr(x, nchar(x)-i, nchar(x)-i))
      colnames(data)[ncol(data)] <- paste("p", poss[length(poss)-i], sep="")
    }
  }
  data$Jmotif = sapply(data$sequence, function(x) substr(x, nchar(x)-4, nchar(x)))
  jmotifs = readRDS(paste(path_to_results, "TCR_structure/Jmotifs.rds", sep="/"))
  data$Jmotif[!(data$Jmotif %in% jmotifs)] = "other"
  if (inc_CDR12){
    cdr1cdr2_ref = readRDS(paste(path_to_data, "CDR1bCDR2b_reference.rds", sep="/"))
    cdr1cdr2_ref$vgene = sapply(cdr1cdr2_ref$vgene, function(x) reformat_vgene(x))
    data = left_join(data, cdr1cdr2_ref[,1:14], by="vgene")
  }
  return(data)
}

add_percentAA <- function(){
  aminos = aminos()
  data$midseq <- sapply(data$sequence, function(x) substr(x, 5, nchar(x)-6))
  for (i in 1:length(aminos)){
    data$new <- sapply(data$midseq, function(x) 100*str_count(x, aminos[i])/nchar(x))
    colnames(data)[ncol(data)] <- paste("perc_mid", aminos[i], sep="_")
  }
  return(data)
}

get_condMImatrix <- function(poss, length, recompute=FALSE){
  if (!recompute){
   mat = readRDS(paste(paste(path_to_results, "conditional_mutual_information/seaydisc_condMI", sep="/"), paste(length, "rds", sep="."), sep="_"))
  } else {
  if (length!="all"){
    sub = data[data$length==length,]
  }
  mat = matrix(nrow=length(poss), ncol=length(poss))
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      max = 0
      for (k in 1:length(poss)){
        if (k==i | k==j){ next }
        indi = which(colnames(sub)==poss[i])
        indj = which(colnames(sub)==poss[j])
        indk = which(colnames(sub)==poss[k])
        cmi = condinformation(X=sub[,indi], Y=sub[,indj], S=sub[,indk])
        entx = condentropy(X=sub[,indi], Y=sub[,indk])
        enty = condentropy(X=sub[,indj], Y=sub[,indk])
        ment = (entx+enty)/2
        if (ment==0) { x = 1 } else { x = cmi/ment }
        if (x >= max){ max = x }
      }
      mat[i,j] = max
    }
  }
  }
  return(mat)
}

one_hot_encode <- function(vars, region="midregion", dataset="seay_disc", recompute=FALSE){
  if (!recompute & region=="Vregion" & dataset=="seay_disc"){
    return(readRDS(paste(path_to_data, "seay2016/seaydisc_Vregion_OHE.rds", sep="/")))
  }
  if (!recompute & region=="Jregion"  & dataset=="seay_disc"){
    return(readRDS(paste(path_to_data, "seay2016/seaydisc_Jregion_OHE.rds", sep="/")))
  }
  if (!recompute & region=="midregion" & dataset=="seay_disc"){
    return(readRDS(paste(path_to_data, "seay2016/seaydisc_midregion_OHE.rds", sep="/")))
  }
  if ("p27" %in% vars){
    data = data[!(is.na(data$p27)),]
  }
  dat = data.frame(matrix(nrow=nrow(data), ncol=0))
  for (i in 1:length(vars)){
    if (grepl("perc_mid", vars[i])){
      dat$new = as.numeric(as.character(data[,which(colnames(data)==vars[i])]))
      colnames(dat)[ncol(dat)] = vars[i]
    } else if (vars[i]=="vgene"){
      vs = unique(data$vgene[data$vgene!="TCRBV05-01"])
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data$vgene==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste("vgene_is", vs[j], sep="_")
      }
    } else if (vars[i]=="jgene"){
      vs = unique(data$jgene[data$jgene!="TCRBJ02-07"])
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data$jgene==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste("jgene_is", vs[j], sep="_")
      }
    } else if (vars[i]=="Vmotif"){
      vmotifs = readRDS(paste(path_to_results, "Vmotifs.rds", sep="/"))
      vs = vmotifs[vmotifs!="CAS"]
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data$Vmotif==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste("Vmotif_is", vs[j], sep="_")
      }
    } else if (vars[i]=="Jmotif"){
      vmotifs = readRDS(paste(path_to_results, "Jmotifs.rds", sep="/"))
      vs = jmotifs[jmotifs!="TEAFF"]
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data$Jmotif==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste("Jmotif_is", vs[j], sep="_")
      }
    } else if (vars[i] == "length"){
      vs = c("12", "13", "14", "16", "17")
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data$length==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste("length_is", vs[j], sep="_")
      }
    } else {
      ref = "G"
      if (vars[i] == "p104"){ ref = "C" }
      if (vars[i] == "p105"){ ref = "A" }
      if (vars[i] == "p106"){ ref = "S" }
      if (vars[i] == "p107"){ ref = "S" }
      if (vars[i] == "p108"){ ref = "L" }
      if (vars[i] == "p113"){ ref = "N" }
      if (vars[i] == "p27"){ ref = "S" }
      if (vars[i] == "p29"){ ref = "H" }
      if (vars[i] == "p37"){ ref = "N" }
      if (vars[i] == "p38"){ ref = "S" }
      if (vars[i] == "p56"){ ref = "S" }
      if (vars[i] == "p57"){ ref = "Q" }
      if (vars[i] == "p58"){ ref = "N" }
      if (vars[i] == "p63"){ ref = "E" }
      if (vars[i] == "p64"){ ref = "E" }
      if (vars[i] == "p65"){ ref = "Q" }
      ind = which(colnames(data)==vars[i])
      obs = unique(data[,ind])
      vs = obs[obs!=ref]
      for (j in 1:length(vs)){
        dat$new = 0
        dat$new[data[,ind]==vs[j]] = 1
        colnames(dat)[ncol(dat)] = paste(vars[i], paste("is", vs[j], sep="_"), sep="_")
      }
    }
  }
  return(dat)
}

get_cormat <- function(vars, region="midregion", dataset="seay_disc", recompute=FALSE){
  if (!recompute & region=="midregion" & dataset=="seay_disc"){
    return(readRDS(paste(path_to_results, "TCRfeat_correlation/seaydisc_midregion_cormat.rds", sep="/")))
  }
  if (!recompute & region=="Vregion" & dataset=="seay_disc"){
    return(readRDS(paste(path_to_results, "TCRfeat_correlation/seaydisc_Vregion_cormat.rds", sep="/")))
  }
  if (!recompute & region=="Jregion" & dataset=="seay_disc"){
    return(readRDS(paste(path_to_results, "TCRfeat_correlation/seaydisc_Jregion_cormat.rds", sep="/")))
  }
  dat = one_hot_encode(vars, region, dataset, recompute)
  return(cor(dat))
}

create_hashmap <- function(key_vector, val_vector){
  h = hash()
  for (i in 1:length(key_vector)){
    h[[key_vector[i]]] = val_vector[i]
  }
  return(h)
}

get_feat_score <- function(x, amap){
  sum = 0
  for (i in 1:nchar(x)){
    sum = sum + amap[[substr(x, i, i)]]
  }
  return(sum/nchar(x))
}

plot_varex <- function(dat, region){
  if (region=="V") { col = "#628d55" }
  if (region=="mid") { col = "#e7862b" }
  if (region=="J") { col = "#cb7392" }
  sub = dat[dat$region==region,]
  srs = unique(sub$series)
  hl= max(sub$ve_norm)
  bonus = 0.3
  if (region=="mid"){ bonus = 5}
  hj = max(sub$step)/10
  g = ggplot()
  sub$label[sub$label=="length" & sub$type=="primary"] = ""
  g = g + geom_hline(yintercept=hl, color=col, size=2)  + geom_line(aes(x=sub$cdf[sub$type=="alternative"], y=sub$ve_norm[sub$type=="alternative"], group=sub$series[sub$type=="alternative"]), linetype="dashed", size=0.5) +  geom_line(aes(x=sub$cdf[sub$type=="primary"], y=sub$ve_norm[sub$type=="primary"], group=sub$series[sub$type=="primary"]), size=0.5) + geom_point(aes(x=sub$cdf, y=sub$ve_norm), color=col) +  geom_point(aes(x=sub$cdf, y=sub$ve_norm), color="black", pch=1) + ylim(c(0,100))  + theme_classic(base_size=12) + xlab("") + ylab("") + xlim(c(-0.3,max(sub$cdf)+bonus))
  return(g + xlab("model degrees of freedom") + ylab("percent of explained variance"))
}

collect_midreg_effectsizes <- function(dir, label){
    setwd(dir)
    files = list.files()
    estimate = vector(mode="numeric")
    term = vector(mode="character")
    std.error = vector(mode="numeric")
    p.value = vector(mode="numeric")
    statistic = vector(mode="numeric")
    c = 1
    for (i in 1:length(files)){
      info = strsplit(gsub(".RData", "", files[i]), "_")[[1]]
      load(files[i])
      ref = "G"
      if (info[4]=="5") { ref = "L" }
      t = t[names(t)!=ref]
      for (j in 1:length(t)){
        if (t[j]/sum(t) > 0.005){
          ind = which(stats$term == paste("loc", names(t)[j], sep=""))
          estimate[c] = stats$estimate[ind]
          std.error[c] = stats$std.error[ind]
          p.value[c] = stats$p.value[ind]
          statistic[c] = stats$statistic[ind]
          term[c] = paste(paste(info[2], info[4], sep="_"), names(t)[j], sep="_")
          c = c + 1
        }
      }
    }
    spos = data.frame(term, estimate, std.error, statistic, p.value)
    colnames(spos)[2:ncol(spos)] = paste(label, colnames(spos[2:ncol(spos)]), sep="_")
    return(spos)
}

convert_to_IMGTpos <- function(pos_betas){
    pos_betas$type = "mid"
    pos_betas$term = as.character(pos_betas$term)
    pos_betas$pos = sapply(pos_betas$term, function(x) strsplit(x, "_")[[1]][2])
    pos_betas$pos[pos_betas$pos=="5"] = "p108"
    pos_betas$pos[pos_betas$pos=="6"] = "p109"
    pos_betas$pos[pos_betas$pos=="7"] = "p110"
    pos_betas$pos[pos_betas$pos=="8"] = "p111"
    pos_betas$pos[pos_betas$pos=="9"] = "p111.0"
    pos_betas$pos[pos_betas$pos=="-8"] = "p111.1"
    pos_betas$pos[pos_betas$pos=="-7"] = "p112"
    return(pos_betas)
}

collect_vj_effectsizes <- function(file, label, region){
    load(file)
    sv = stats[,3:7]
    if (region=="V"){
        rare = names(t4)[t4/sum(t4) < 0.005]
        sv = sv[!(sv$term %in% paste("pos4", rare, sep="")),]
        rare = names(tvg)[tvg/sum(tvg) < 0.005]
        sv = sv[!(sv$term %in% paste("vgene", rare, sep="")),]
    } else if (region=="J" & label=="GT"){
        rare = names(t)[t/sum(t) < 0.005]
        sv = sv[!(sv$term %in% paste("pos.6", rare, sep="")),]
    }
    colnames(sv)[2:ncol(sv)] = paste(label, colnames(sv)[2:ncol(sv)], sep="_")
    sv = sv[!(grepl("vgeneF", sv$term)),]
    sv = sv[!(grepl("length", sv$term)),]
    sv = sv[!(grepl("Intercept", sv$term)),]
    return(sv)
}

read_heldout_TiRP <- function(){
    seay = readRDS(paste(path_to_results, "TiRP_scoring/seay2016_heldout_TiRP.rds", sep="/"))
    GT = readRDS(paste(path_to_results, "TiRP_scoring/gt2017_TiRP.rds", sep="/"))
    GT = GT[GT$donor %in% c("HD1", "HD2", "T1D6"),]
    GT$sequence = as.character(GT$CDR3)
    GT$site = "PBMC"
    GT$cell_state = "Tconv"
    GT$cell_state[GT$celltype=="Treg"] = "Treg"
    seay$dataset = "seay"
    GT$dataset = "GT"
    seay$total_score = seay$length_score + seay$vgene_score + seay$pos4_score + seay$perc_score + seay$pos_score + seay$Jseq_score + seay$'pos.6_score'
    cn = c("sequence", "total_score", "donor", "dataset", "cell_state", "site")
    bulkscored = rbind(seay[,cn], GT[,cn])
    bulkscored = bulkscored[!(duplicated(bulkscored)),]
    return(bulkscored)
}

get_percentile <- function(x){
  return(as.numeric(as.character(gsub("%", "", names(ptls[min(which(ptls>=x))])))))
}


plot_adtlevels <- function(feature, pct = 0.95) {
  ind = which(rownames(prot)==feature)
  vals = as.numeric(as.character(prot[ind,]))
  max.cutoff = quantile(vals, pct)
  min.cutoff = quantile(vals, 1-pct)
  tmp <- sapply(X = vals, FUN = function(x) {
    return(ifelse(test = x > max.cutoff, yes = max.cutoff,
                  no = x))
  })
  tmp <- sapply(X = tmp, FUN = function(x) {
    return(ifelse(test = x < min.cutoff, yes = min.cutoff,
                  no = x))
  })
  umap_res_plot <- cbind(rmd, tmp)
  return(ggplot(data = as.data.frame(umap_res_plot)[sample(nrow(umap_res_plot)),] , aes(x = UMAP_1, y = UMAP_2)) +
           geom_point_rast(mapping = aes(color = tmp), shape = ".", show.legend=FALSE) +
           scale_color_viridis(end =0.9) +
           theme_classic() + ggtitle(feature))
}
