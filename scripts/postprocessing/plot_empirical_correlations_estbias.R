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
emp.results <- dnds.results %>% filter(type == "empirical")  %>% select(-type) %>% spread(method,dnds) %>% filter(slac <= 1) %>% gather(method,dnds,d0.01:slac) 

emp.r.bias <- data.frame("dataset" = as.character(),
                         "method"  = as.character(),
                         "corr"    = as.numeric(),
                         "estbias" = as.numeric())
datasets <- unique(emp.results$dataset)
for (d in datasets){
   emp.results %>% filter(dataset == d) -> temp
   print(d)
   for (m in c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000", "pbmutsel")){
       
       r <- cor(temp$dnds[temp$method == "slac"], temp$dnds[temp$method == m]) 
       braw <- lm(temp$dnds[temp$method == m] ~ offset(temp$dnds[temp$method == "slac"]))
       b <- summary(braw)$coeff[1]
       
       temp.r.bias <- data.frame("dataset" = d,
                                 "method"  = m,
                                 "corr"    = r,
                                 "estbias" = b)
       emp.r.bias <- rbind(emp.r.bias, temp.r.bias)
       
}}

# multivariate penalties are basically the same as no penalty
emp.r.bias %>% filter(method %in% c("nopenal", "d0.01", "d0.1", "d1.0", "pbmutsel")) -> sub.emp.r.bias

sub.emp.r.bias$method <- factor(sub.emp.r.bias$method)
colors <- c("skyblue", "dodgerblue", "blue", "navy", "red")

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size=7,angle=30))) 

correlations <- ggplot(sub.emp.r.bias, aes(x = dataset, y = corr, color = method)) + geom_point(size=3, alpha=0.7) + xlab("dataset") + ylab("correlation with slac1") + scale_color_manual(values = colors)

estbiases <- ggplot(sub.emp.r.bias, aes(x = dataset, y = estbias, color = method)) + geom_hline(yintercept=0) + geom_point(size=3, alpha=0.7) + xlab("dataset") + ylab("estimator bias") + scale_color_manual(values = colors)

both <- plot_grid(correlations, estbiases, nrow=1)
save_plot("empirical_correlations_estimatorbias.pdf", both, base_width=8, base_height=3)
both