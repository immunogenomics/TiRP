source("utils.R")

data <- readRDS("data/seay2016/seay_discovery.rds")

args = commandArgs(trailingOnly=TRUE)
i = as.numeric(as.character(args[1]))

if (i==1){
  vars = c("p27", "p28", "p29", "p37", "p38", "p56", "p57", "p58", "p63", "p64", "p65", "vgene", "Vmotif", "p105", "p106", "p107")
  tag = "Vregion"
}

if (i==2){
  vars = c("Jmotif", "jgene", "p113", "p114", "p115", "p116", "p117", "p118")
  tag =  "Jregion"
}

if (i==3){
  vars = c(colnames(data)[grepl("perc_mid", colnames(data))], "length", "p108", "p109", "p110", "p111", "p111.1", "p112.1", "p112")
  tag = "midregion"
}

dat = one_hot_encode(vars, region=tag)
setwd("results")
saveRDS(cor(dat), paste(tag, "cormat.rds", sep="_"))
