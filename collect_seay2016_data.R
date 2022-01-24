
setwd("~/TiRP_github/downloads/immunoSEQ_seay2016")
files = list.files()

data = NULL
for (i in 1:length(files)){
  dat = read.table(files[i], sep="\t", header=TRUE)
  info = strsplit(files[i], "_")[[1]]
  if (info[1]=="6323"){
    dat$donor = "6323"
    dat$pheno = "T1D"
    dat$site = info[2]
    dat$cell_state = gsub(".tsv", "", info[3])
  } else {
    dat$donor = info[2]
    dat$pheno = info[1]
    dat$site = info[3]
    if (length(info)==3) { 
      dat$cell_state = "unsorted"
    } else {
      dat$cell_state = gsub(".tsv", "", info[4])
    }
  }
  data = rbind(dat, data)
}

data$locus = "TCRB"
data$frame_type = as.character(data$sequenceStatus)
data$sequence = as.character(data$aminoAcid)
data$vgene = as.character(data$vGeneName)
data$jgene = as.character(data$jGeneName)
data = data[data$cell_state!="unsorted",]
save(data[,c("cell_state", "donor", "pheno", "site", "locus", "sequence", "frame_type", "vgene", "jgene")], "~/TiRP_github/data/seay2016/seay_data.rds")

