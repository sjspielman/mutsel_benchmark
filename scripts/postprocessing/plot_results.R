require(dplyr)
require(tidyr)
require(cowplot)
require(readr)
require(reshape2)

results <- read_csv("../../results/results.csv")
results %>% spread(method, dnds) -> results.spread
results %>% filter(type == "simulation") -> sim.results
results %>% filter(type == "empirical") -> emp.results


#### 






















##### Correlation boxplots for simulated data #####
simcorrs <- NULL
simnames <- unique(sim.results$dataset)
first <- TRUE
i <- 1
for (dataname in simnames){
    newdname <- strsplit(dataname, "_")[[1]][1]
    print(newdname)
    temp <- results.spread %>% filter(dataset == dataname, fel1 != 1) %>% select(-dataset, -type, -site) %>% na.omit()
    corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 %in% c("empirical", "true", "fel1"), value != 1) %>% arrange(Var1)
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
simcorrs$Var2 <-  factor(simcorrs$Var2, levels = c("true", "fel1", "empirical", "nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0"))
simcorrs %>% filter(Var1 == "true", Var2 %in% c("fel1", "empirical", "nopenal", "mvn10", "d1.0") %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("Inference method") + ylab("Pearson Correlation")


##### Scatterplots for a representative dataset, 2BCG_Y_simulated #####
representative_sim <- "2BCG_Y_simulated"
results.spread %>% filter(dataset == representative_sim, fel1 != 1)  -> plotme # fel1 are not informative.
r.true_fel <- round(cor(plotme$true, plotme$fel1),2)
r.true_empirical <- round(cor(plotme$true, plotme$empirical),2)
r.fel_empirical <- round(cor(plotme$fel1, plotme$empirical),2)
r.true_nopenal <- round(cor(plotme$true, plotme$nopenal),2)
r.true_mvn10 <- round(cor(plotme$true, plotme$mvn10),2)
r.true_d1.0 <- round(cor(plotme$true, plotme$d1.0),2)
#r.true_pb <- round(cor(plotme$true, plotme$pb),2)
p1 <- ggplot(plotme, aes(x = true, y = fel1)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("FEL1") + ggtitle(paste(dataname,r))
p2 <- ggplot(plotme, aes(x = true, y = empirical)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("Empirical dN/dS")
p3 <- ggplot(plotme, aes(x = true, y = nopenal)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("swMutSel nopenal dN/dS")
plot.new()
p<-plot_grid(p1, p2, p3, nrow=1)
print(p)














##### Correlation boxplots for empirical data #####
empcorrs <- NULL
empnames <- unique(emp.results$dataset)
first <- TRUE
i <- 1
for (dataname in empnames){
    temp <- results.spread %>% filter(dataset == dataname) %>% select(-dataset, -type, -site, -true) %>% na.omit()
    corrs <- as.matrix(cor(temp)) %>% melt() %>% filter(Var1 %in% c("empirical", "fel1"), value != 1) %>% arrange(Var1)
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
empcorrs$Var2 <-  factor(empcorrs$Var2, levels = c("fel1", "empirical", "nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0"))
empcorrs %>% filter(Var1 == "fel1", Var2 %in% c("empirical", "nopenal", "mvn10", "d1.0")) %>% ggplot(aes(x = Var2, y = r)) + geom_boxplot() + xlab("Inference method") + ylab("Pearson Correlation") # outliers are either amine or pb2





##### Scatterplots for a representative dataset #####
d <- "PF07715"
results.spread %>% filter(dataset == d)  -> plotme
p1 <- ggplot(plotme, aes(x = fel1, y = empirical)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL1 dN/dS") + ylab("Empirical dN/dS") + ggtitle(d) + scale_y_log10() + scale_x_log10()
p2 <- ggplot(plotme, aes(x = fel1, y = nopenal)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("FEL1 dN/dS") + ylab("swMutSel nopenal dN/dS")
p3 <- ggplot(plotme, aes(x = empirical, y = nopenal)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("Empirical dN/dS") + ylab("swMutSel nopenal dN/dS")
plot_grid(p1, p2, p3, nrow=1)








