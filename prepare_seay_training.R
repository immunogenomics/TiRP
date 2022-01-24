source("../utils.R")

data <- readRDS(paste(path_to_data, "seay2016/seay_discovery.rds", sep="/"))

##to focus on unambiguous Treg and Tconv TCRs, we filter from this dataset the <5% of TCRs that are observed in both the Treg and Tconv 

data$ident <- paste(data$sequence, data$donor)
ex = data %>% group_by(ident) %>% summarise(nTconvs = length(cell_state[cell_state=="Tconv"]), nTregs = length(cell_state[cell_state=="Treg"]))
ex$clonetype = "mixed"
ex$clonetype[ex$nTconvs==0] = "Treg"
ex$clonetype[ex$nTregs==0] = "Tconv"

nrow(data[(data$ident %in% ex$ident[ex$clonetype=="mixed"]),])/nrow(data)

data <- data[!(data$ident %in% ex$ident[ex$clonetype=="mixed"]),]

data <- add_percentAA()

saveRDS(data, paste(path_to_data, "seay2016/seay_training.rds", sep="/"))


