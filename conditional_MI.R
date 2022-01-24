source("utils.R")

data <- readRDS("data/seay_discovery.rds")

args = commandArgs(trailingOnly=TRUE)
i = as.numeric(as.character(args[1]))

lens = c("all", "12", "13", "14", "15", "16", "17")

seq_length = lens[i]

if (seq_length!="all"){
    poss = paste("p", get_positions(seq_length), sep="")
    sub = data[data$length==seq_length,]
  } else {
    poss = c("vgene", "p27", "p28", "p29", "p37", "p38", "p56", "p57", "p58", "p63", "p64", "p65", "Vmotif", "p107", "p113", "Jmotif", "jgene")
    sub = data[!(is.na(data$cdr1)),]
  }

mat = get_condMImatrix(poss, seq_length, recompute=TRUE)

saveRDS(mat, paste(paste("results/condMI", seq_length, sep="_"), "rds", sep="."))
