require(dplyr)
require(tidyr)
require(cowplot)
require(readr)
require(reshape2)
require(grid)
require(lme4)
require(lmerTest)

result_directory <- "../../results/"

dnds.results <- read_csv(paste0(result_directory,"dnds_results.csv"))
sim.results <- dnds.results %>% filter(type == "simulation") %>% select(-type)
emp.results <- dnds.results %>% filter(type == "empirical")  %>% select(-type)
jsd.results <- read_csv(paste0(result_directory, "jsd_results.csv"))


#### Boxplots of Jensen-Shannon distance among methods for simulated datasets ####
sim.results <- left_join(sim.results, jsd.results) 
sim.results %>% filter(method != "fel1") %>% select(-dnds) %>% na.omit() %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.data
jsd.data$method <- factor(jsd.data$method, levels = c("nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbmutsel"), labels = c("No penalty", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbMutSel"))
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 13),  axis.title = element_text(size = 16)))
jsd.boxplots <- ggplot(jsd.data, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4))
save_plot("plots/jsd_boxplot_raw.pdf", jsd.boxplots, base_width = 8, base_height = 4.5 )



##### Does level of constraint correlate with error? #####
# Marginally. Error is a bit higher when constraint is lower, but the effect is minimal and probably not worth mentioning.
#sim.results %>% filter(method == "true") %>% spread(method, dnds) %>% select(-jsd) -> sim.true
#sim.results %>% filter(method == "nopenal") %>% spread(method, dnds) -> nopenal.jsd
#nopenal.jsd <- left_join(nopenal.jsd, sim.true)
#summary(lmer(jsd ~ true + 1(|dataset), data = nopenal.jsd))
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept) 5.325e-02  1.679e-03 8.000e+01  31.709  < 2e-16 ***
#     true    2.460e-02  4.423e-03 3.682e+03   5.562 2.86e-08 ***
#ggplot(nopenal.jsd, aes(x = true, y = jsd)) + geom_point() 



##### Scatterplots for a representative dataset #####
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15)))
representative_sim <- "1IBS_A_simulated" # Note that this dataset was chosen because it has the most sites
plotme <- sim.results %>% select(-jsd) %>% filter(dataset == representative_sim) %>% spread(method, dnds) 
p1 <- ggplot(plotme, aes(x = true, y = fel1)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("FEL dN/dS") + scale_y_continuous(limits=c(0,1))
p2 <- ggplot(plotme, aes(x = true, y = nopenal)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1))
p3 <- ggplot(plotme, aes(x = true, y = d1.0)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("swMutSel d1.0 dN/dS") + scale_y_continuous(limits=c(0,1)) + theme(axis.title.y = element_text(size = 13))
p4 <- ggplot(plotme, aes(x = true, y = pbmutsel)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1))
sim.scatterplots <- plot_grid(p1, p2, p3, p4, nrow=1, labels=c("A", "B", "C", "D"), label_size = 15)
save_plot("plots/scatterplots_simulated.pdf", sim.scatterplots, base_width = 12, base_height = 3 )



##### Frequency comparison barplots for representative dataset #####
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 10),  axis.title = element_text(size = 16), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=14, face="bold")))
aafreq <- read.csv("aa_frequency_comparison_1IBS_A_simulated.csv") # Same representative dataset as in scatterplots
keep <- c(283,176,56)  # here, define sites to plot.
aafreq %>% filter(site %in% keep) -> sub_aafreq
sub_aafreq$truednds <- round(sub_aafreq$truednds, 4)
sub_aafreq <- mutate(sub_aafreq, true_facet = paste0("dN/dS = ",truednds))
sub_aafreq$method <- factor(sub_aafreq$method, levels=c("true", "mvn10", "phylobayes"), labels = c("True", "swMutSel", "pbMutSel"))
sub_aafreq$freq <- sub_aafreq$freq + 0.001 # For clarity in plot

sub_aafreq %>% ggplot(aes(x = aminoacid, y = freq, group=method, fill=method)) + geom_bar(stat="identity",position="dodge",width=0.75) + facet_grid(~true_facet) + ylab("Equilibrium frequency") + xlab("Amino Acid") + scale_fill_manual(name = "", values = c("black", "firebrick1", "mediumblue")) -> freq_barplot
save_plot("plots/aafreq_barplot_1IBS_A.pdf", freq_barplot, base_width = 11.5, base_height = 3.75)


##### Scatterplots for empirical datasets #####

# These lines of code give some summary statistics to decide what to use in scatterplot. We want the least and most biased.
# emp.results %>% spread(method, dnds) %>% group_by(dataset) %>% do(bias.raw = glm(nopenal ~ offset(fel1), dat = .), cor.raw = cor(.$fel1, .$nopenal)) %>% mutate(estbias = summary(bias.raw)$coeff[1], r = cor.raw[1]) %>% select(-bias.raw, -cor.raw) -> emp.r.bias
# emp.stats <- read_csv("empirical_data_statistics.csv")
# emp.stats <- left_join(emp.r.bias, emp.stats)

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)))

representative_emp <- "PF00593" # least estimator bias = -0.034, r = 0.84
p1 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = nopenal)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))

p2 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


representative_emp <- "PF04055" # highest estimator bias 0.26, r = 0.73
p3 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = nopenal)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
p4 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


representative_emp <- "pb2" # everybody has analyzed this, so we should too
p5 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = nopenal)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
p6 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


p<-plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, labels = c("A", "B", "C", "D", "E", "F"), label_size=17, scale=0.95)
save_plot("plots/scatterplots_empirical_raw.pdf", p, base_width = 8, base_height = 10)







######## OLD PLOTS ARE DOWN HERE!! ##########


##### Correlation boxplots for simulated data #####
# simcorrs <- NULL
# simnames <- unique(sim.results$dataset)
# first <- TRUE
# i <- 1
# for (dataname in simnames){
#   newdname <- strsplit(dataname, "_")[[1]][1]
#   temp <- sim.results %>% spread(method, dnds) %>% filter(dataset == dataname, fel1 != 1) %>% select(-dataset, -site) %>% na.omit()
#   corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 %in% c("true", "fel1"), value != 1) %>% arrange(Var1)
#   corrs$dataset <- newdname
#   colnames(corrs)[3] <- "r"
#   if (first){
#     simcorrs <- corrs 
#   }
#   else{
#     simcorrs <- rbind(simcorrs, corrs)
#   }
#   first <- FALSE
#   i <- i+1
# }
# # correlations between true and inferred
# #      Var2 median(r)
# # 1   d0.01 0.9854866
# # 2    d0.1 0.9452034
# # 3    d1.0 0.8427745
# # 5   mvn10 0.9870753  # best but not far off
# # 6  mvn100 0.9865595
# # 7 mvn1000 0.9865439
# # 8 nopenal 0.9865319
# simcorrs$Var2 <-  factor(simcorrs$Var2, levels = c("fel1", "nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbmutsel"), labels = c("FEL", "No penalty", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbMutSel"))
# theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 12),  axis.title = element_text(size = 16)))
# simcorrs %>% filter(Var1 == "true") %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("\nInference method") + ylab("Pearson Correlation") -> sim.boxplots
# save_plot("plots/correlation_boxplots_simulated_raw.pdf", sim.boxplots, base_width = 8, base_height = 4.5 )


# ##### Correlation boxplots for empirical data #####
# empcorrs <- NULL
# empnames <- unique(emp.results$dataset)
# first <- TRUE
# i <- 1
# for (dataname in empnames){
#   temp <- emp.results %>% filter(dataset == dataname) %>% spread(method, dnds) %>% select(-dataset, -site) %>% na.omit()
#   corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 == "fel1", value != 1) %>% arrange(Var1)
#   corrs$dataset <- dataname
#   colnames(corrs)[3] <- "r"
#   if (first){
#     empcorrs <- corrs 
#   }
#   else{
#     empcorrs <- rbind(empcorrs, corrs)
#   }
#   first <- FALSE
#   i <- i+1
# }
# 
# empcorrs$Var2 <-  factor(empcorrs$Var2, levels = c("fel1", "nopenal", "pbmutsel"), labels = c("FEL", "swMutSel", "pbMutSel"))
# theme_set(theme_cowplot() + theme(axis.text = element_text(size = 11),  axis.title = element_text(size = 13)))
# empcorrs %>% filter(Var1 == "fel1") %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("Inference method") + ylab("Spearman Correlation") + scale_y_continuous(limits=c(0,1), labels=c(0.1, 0.3, 0.5, 0.7, 0.9)) -> emp.boxplots
# save_plot("plots/correlation_boxplots_empirical_raw.pdf", emp.boxplots, base_width = 4, base_height = 3)

#emp.stats <- read_csv("empirical_data_statistics.csv")
#empcorrs <- left_join(empcorrs, emp.stats)
#empcorrs %>% filter(Var2 == "swMutSel", dataset !="pb2")-> sanspb2
#summary(lm(r~meanpairwise, data=sanspb2))



