require(dplyr)
require(tidyr)
require(cowplot)
require(readr)
#require(reshape2)
require(grid)

result_directory <- "../../results/"

results <- read_csv(paste0(result_directory,"results.csv"))
sim.results <- results %>% filter(type == "simulation") %>% select(-type)
emp.results <- results %>% filter(type == "empirical")  %>% select(-type)
#sim.results$method <- factor(levels = c("true", "fel1", "nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbmutsel"), labels = c("True", "FEL", "No penalty", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbMutSel"))
#emp.results$method <- factor(levels = c("fel1", "mvn10", "pbmutsel"), labels = c("FEL", "swMutSel", "pbMutSel"))



##### Correlation boxplots for simulated data #####
simcorrs <- NULL
simnames <- unique(sim.results$dataset)
first <- TRUE
i <- 1
for (dataname in simnames){
    newdname <- strsplit(dataname, "_")[[1]][1]
    temp <- sim.results %>% spread(method, dnds) %>% filter(dataset == dataname, fel1 != 1) %>% select(-dataset, -site) %>% na.omit()
    corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 %in% c("true", "fel1"), value != 1) %>% arrange(Var1)
    corrs$dataset <- newdname
    colnames(corrs)[3] <- "r"
    if (first){
        simcorrs <- corrs 
    }
    else{
        simcorrs <- rbind(simcorrs, corrs)
    }
    first <- FALSE
    i <- i+1
}
# correlations between true and inferred
#      Var2 median(r)
# 1   d0.01 0.9854866
# 2    d0.1 0.9452034
# 3    d1.0 0.8427745
# 5   mvn10 0.9870753  # best but not far off
# 6  mvn100 0.9865595
# 7 mvn1000 0.9865439
# 8 nopenal 0.9865319
simcorrs$Var2 <-  factor(simcorrs$Var2, levels = c("fel1", "nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbmutsel"), labels = c("FEL", "No penalty", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbMutSel"))


theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 12),  axis.title = element_text(size = 16))
simcorrs %>% filter(Var1 == "true") %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("\nInference method") + ylab("Pearson Correlation") -> sim.boxplots
save_plot("correlation_boxplots_simulated_raw.pdf", sim.boxplots, base_width = 8, base_height = 4.5 )





##### Scatterplots for a representative dataset #####
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15)))
representative_sim <- "1IBS_A_simulated" # Note that this dataset was chosen because it has the most sites
plotme <- sim.results %>% filter(dataset == representative_sim) %>% spread(method, dnds)
p1 <- ggplot(plotme, aes(x = true, y = fel1)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True") + ylab("FEL") + scale_y_continuous(limits=c(0,1))
p2 <- ggplot(plotme, aes(x = true, y = mvn10)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True") + ylab("swMutSel mvn10") + scale_y_continuous(limits=c(0,1))
p3 <- ggplot(plotme, aes(x = true, y = d1.0)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True") + ylab("swMutSel d1.0") + scale_y_continuous(limits=c(0,1))
p4 <- ggplot(plotme, aes(x = true, y = pbmutsel)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True") + ylab("pbMutSel") + scale_y_continuous(limits=c(0,1))
sim.scatterplots <- plot_grid(p1, p2, p3, p4, nrow=1, labels=c("A", "B", "C", "D"), label_size = 16)
save_plot("scatterplots_simulated_raw.pdf", sim.scatterplots, base_width = 12, base_height = 3 )




##### Correlation boxplots for empirical data #####
empcorrs <- NULL
empnames <- unique(emp.results$dataset)
first <- TRUE
i <- 1
for (dataname in empnames){
    temp <- emp.results %>% filter(dataset == dataname) %>% spread(method, dnds) %>% select(-dataset, -site) %>% na.omit()
    corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 == "fel1", value != 1) %>% arrange(Var1)
    corrs$dataset <- dataname
    colnames(corrs)[3] <- "r"
    if (first){
        empcorrs <- corrs 
    }
    else{
        empcorrs <- rbind(empcorrs, corrs)
    }
    first <- FALSE
    i <- i+1
}

empcorrs$Var2 <-  factor(empcorrs$Var2, levels = c("fel1", "mvn10", "pbmutsel"), labels = c("FEL", "swMutSel", "pbMutSel"))
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 11),  axis.title = element_text(size = 13)))
empcorrs %>% filter(Var1 == "fel1") %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("Inference method") + ylab("Pearson Correlation") + scale_y_continuous(limits=c(0,1), labels=c(0.1, 0.3, 0.5, 0.7, 0.9)) -> emp.boxplots
save_plot("correlation_boxplots_empirical_raw.pdf", emp.boxplots, base_width = 4, base_height = 3)


emp.stats <- read_csv("empirical_data_statistics.csv")
empcorrs <- left_join(empcorrs, emp.stats)
#empcorrs %>% filter(Var2 == "swMutSel", dataset !="pb2")-> sanspb2
#summary(lm(r~meanpairwise, data=sanspb2))


##### Scatterplots for empirical datasets #####
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 17)))

representative_emp <- "PF00593"
p1 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = mvn10)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("swMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + annotate("text", x = 0.8, y = )

p2 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("pbMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


representative_emp <- "PF01926"
p3 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = mvn10)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("swMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
p4 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("pbMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


representative_emp <- "pb2"
p5 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = mvn10)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("swMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
p6 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = fel1, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL") + ylab("pbMutSel") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))


p<-plot_grid(p1, p2, p3, p4, p5, p6, nrow=3, labels = c("A", "B", "C", "D", "E", "F"), label_size=16)







