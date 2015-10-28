require(dplyr)
require(tidyr)
require(cowplot)
require(readr)
require(reshape2)
require(grid)
require(lme4)
require(lmerTest)

# datasets <- unique(emp.results$dataset)
# for (d in datasets){
#     emp.results %>% filter(dataset == d) -> temp
#     print(d)
#     for (m in c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000", "pbmutsel")){
#         
#         r <- cor(temp$dnds[temp$method == "slac"], temp$dnds[temp$method == m]) 
#         braw <- lm(temp$dnds[temp$method == m] ~ offset(temp$dnds[temp$method == "slac"]))
#         b <- summary(braw)$coeff[1]
#         print(paste(m, r, b))
# }}


result_directory <- "../../results/"

dnds.results <- read_csv(paste0(result_directory,"dnds_results.csv"))
sim.results <- dnds.results %>% filter(type == "simulation") %>% select(-type)
emp.results <- dnds.results %>% filter(type == "empirical")  %>% select(-type)
jsd.results <- read_csv(paste0(result_directory, "jsd_results.csv"))


#### Boxplots of Jensen-Shannon distance among methods for simulated datasets ####
left_join(sim.results, jsd.results)  %>% filter(method != "slac") %>% select(-dnds) %>% na.omit() %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.data
jsd.data$method <- factor(jsd.data$method, levels = c("nopenal", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbmutsel"), labels = c("No penalty", "mvn10", "mvn100", "mvn1000", "d0.01", "d0.1", "d1.0", "pbMutSel"))
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 13),  axis.title = element_text(size = 16)))
jsd.boxplots <- ggplot(jsd.data, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2))
save_plot("plots/jsd_boxplot_raw.pdf", jsd.boxplots, base_width = 8, base_height = 4.5 )



##### Scatterplots for a representative dataset #####
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15)))
representative_sim <- "1ZNN_A_simulated" # Chosen since its r, bias are most similar to average for pbMutSel
plotme <- sim.results %>% filter(dataset == representative_sim) %>% spread(method, dnds) 
p1 <- ggplot(plotme, aes(x = true, y = slac)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("SLAC dN/dS") + scale_y_continuous(limits=c(0,1))
p2 <- ggplot(plotme, aes(x = true, y = nopenal)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1))
p3 <- ggplot(plotme, aes(x = true, y = d1.0)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("swMutSel d1.0 dN/dS") + scale_y_continuous(limits=c(0,1)) + theme(axis.title.y = element_text(size = 13))
p4 <- ggplot(plotme, aes(x = true, y = pbmutsel)) + geom_point(alpha=0.6) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1))
sim.scatterplots <- plot_grid(p1, p2, p3, p4, nrow=1, labels=c("A", "B", "C", "D"), label_size = 15)
save_plot("plots/scatterplots_simulated_1ZNN.pdf", sim.scatterplots, base_width = 12, base_height = 3 )




##### Violin plots of correlation and estimator bias across simulated datasets #####
sim.results %>% spread(method, dnds) %>% group_by(dataset) %>% do(bias.raw = glm(pbmutsel ~ offset(true), dat = .), cor.raw = cor(.$true, .$pbmutsel)) %>% mutate(bias = summary(bias.raw)$coeff[1], r2 = (cor.raw[1])^2) %>% select(-bias.raw, -cor.raw) %>% ungroup() %>% mutate(method = "pbMutSel", meanbias = mean(bias), meanr2 = mean(r2)) -> sim.r2.bias.pb
sim.results %>% spread(method, dnds) %>% group_by(dataset) %>% do(bias.raw = glm(nopenal ~ offset(true), dat = .), cor.raw = cor(.$true, .$nopenal)) %>% mutate(bias = summary(bias.raw)$coeff[1], r2 = (cor.raw[1])^2) %>% select(-bias.raw, -cor.raw) %>% ungroup() %>% mutate(method = "swMutSel", meanbias = mean(bias), meanr2 = mean(r2)) -> sim.r2.bias.sw
sim.results %>% spread(method, dnds) %>% group_by(dataset) %>% do(bias.raw = glm(slac ~ offset(true), dat = .), cor.raw = cor(.$true, .$slac)) %>% mutate(bias = summary(bias.raw)$coeff[1], r2 = (cor.raw[1])^2) %>% select(-bias.raw, -cor.raw) %>% ungroup() %>% mutate(method = "SLAC", meanbias = mean(bias), meanr2 = mean(r2)) -> sim.r2.bias.slac
sim.r2.bias <- rbind(sim.r2.bias.pb, sim.r2.bias.sw, sim.r2.bias.slac)
sim.r2.bias$method <- factor(sim.r2.bias$method, levels=c("SLAC", "swMutSel", "pbMutSel"))

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 13), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 17)))
r2.violin <- ggplot(sim.r2.bias, aes(x = method, y= r2)) + geom_violin(scale="width") + geom_point(aes(x = method, y = meanr2), size=3) + xlab("Inference Method") + ylab("Variance Explained") + scale_y_continuous(limits = c(0.6,1.0), breaks=c(0.6,0.7,0.8,0.9,1.0))
bias.violin <- ggplot(sim.r2.bias, aes(x = method, y= bias)) + geom_violin() + geom_point(aes(x = method, y = meanbias),size=3) + xlab("Inference Method") + ylab("Estimator Bias") + geom_hline(yintercept=0, color="grey40", size=0.75) + scale_y_continuous(breaks = c(-0.05, -0.025, 0, 0.025, 0.05, 0.075))
violins.r2.bias <- plot_grid(r2.violin, bias.violin, nrow=1, labels = c("A","B"), label_size=18, scale=0.95)
save_plot("plots/violins_r2_bias.pdf", violins.r2.bias, base_width=9, base_height=4)


# summary(lm(abs(bias) ~ method, data = sim.r2.bias))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.017635   0.002057   8.574 7.71e-12 ***
#   methodswMutSel -0.010074   0.002909  -3.463  0.00102 ** 
#   methodpbMutSel  0.055891   0.002909  19.215  < 2e-16 ***
# 
# summary(lm(r2 ~ method, data = sim.r2.bias))
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.811105   0.012083  67.125  < 2e-16 ***
#   methodswMutSel  0.161260   0.017089   9.437 3.01e-13 ***
#   methodpbMutSel -0.002122   0.017089  -0.124    0.902  





##### Scatterplots for empirical datasets #####
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15)))

plot_empirical_scatter_row <- function(representative_emp)
{
  p1 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = slac, y = nopenal)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("SLAC dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
  
  p2 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = slac, y = pbmutsel)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("SLAC dN/dS") + ylab("pbMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
  
  p3 <- emp.results %>% filter(dataset == representative_emp) %>% spread(method, dnds) %>% ggplot(aes(x = pbmutsel, y = nopenal)) + geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("pbMutSel dN/dS") + ylab("swMutSel dN/dS") + scale_y_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0)) + scale_x_continuous(limits=c(0,1.1), breaks=c(0, 0.25, 0.50, 0.75, 1.0))
  
  list(p1,p2,p3)
}

row1 <- plot_empirical_scatter_row("PF00593") # swMutSel better than pbMutSel
row2 <- plot_empirical_scatter_row("PF04055") # pbMutSel better than swMutSel
row3 <- plot_empirical_scatter_row("pb2")     # commonly analyzed
p <- plot_grid(row1[[1]], row1[[2]], row1[[3]], row2[[1]], row2[[2]], row2[[3]], row3[[1]], row3[[2]], row3[[3]], nrow=3, ncol=3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), label_size=19, scale=0.95)
save_plot("plots/scatterplots_empirical_raw.pdf", p, base_width = 13, base_height = 10)
