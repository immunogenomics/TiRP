
source("utils.R")

heldout = read_heldout_TiRP()
ptls <- quantile(heldout$total_score, probs=seq(0.01, 1, by=0.01))

heldout$Treg_percentile = sapply(heldout$total_score, function(x) get_percentile(x))
grouped = heldout %>% group_by(Treg_percentile) %>% dplyr::summarise(Treg_odds = length(cell_state[cell_state=="Treg"])/length(cell_state[cell_state=="Tconv"]))
grouped$Treg_percentile = grouped$Treg_percentile - 1

saveRDS(grouped, paste(path_to_results, "TiRP_scoring/heldout_TiRP_decilegroups.rds", sep="/"))


