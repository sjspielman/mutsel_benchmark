require(cowplot)
require(dplyr)
dat <- read.csv("aa_frequency_comparison_1B4T_A_simulated.csv")
# keep<-c(93,114,126)  # here, define sites to plot. Three would be nice (low, med, high dnds)
dat %>% filter(site %in% keep) -> subdat
subdat$truednds <- round(subdat$truednds, 4)

subdat %>% ggplot(aes(x = aminoacid, y = freq, group=method, fill=method)) + geom_bar(stat="identity",position="dodge") + facet_grid(~truednds) + theme(axis.text = element_text(size = 8))
