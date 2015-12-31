require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(grid)


repr_sim <- "1B4T_A"

result_directory <- "../../results/summarized_results/"

jsd.results <- read_csv(paste0(result_directory, "simulation_jsd.csv"))

theme_set(theme_cowplot() + theme(panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12)))

#### Boxplots of Jensen-Shannon distance among methods for simulated datasets ####
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 13),  axis.title = element_text(size = 16)))

jsd.results %>% filter(del == "strong") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.strong.summary
jsd.strong.summary$method <- factor(jsd.strong.summary$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))
jsd.results %>% filter(del == "weak") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.weak.summary
jsd.weak.summary$method <- factor(jsd.weak.summary$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))
jsd.results %>% filter(del == "strong", dataset == repr_sim) -> jsd.repr.sim
jsd.repr.sim$method <- factor(jsd.repr.sim$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))



## Figure 1A 
repr.jsd.boxplots <- ggplot(jsd.repr.sim, aes(x = method, y = jsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Site JSD") + scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) + ggtitle("JSD for representative dataset")
## Figure 1B
mean.jsd.strong.boxplots <- ggplot(jsd.strong.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average site JSD") + scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) + ggtitle("Mean JSD across datasets")
## Full Figure 1
fig1 <- plot_grid(repr.jsd.boxplots, mean.jsd.strong.boxplots, nrow=2, labels=c("A", "B"), scale=0.95)



#### SI weak jsd boxplots
mean.jsd.weak.boxplots <- ggplot(jsd.weak.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))



##### Figures 2,3: dN/dS scatterplots and boxplots for simulated datasets #####
sim.dat <- read.csv(paste0(result_directory, "simulation_derived_dnds.csv"))
sim.dat %>% spread(method,dnds) %>% gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, site, true, dnds, method)-> dnds

dnds.repr.strong <- dnds %>% filter(dataset == repr_sim, del == "strong")
dnds.repr.strong$method <- factor(dnds.repr.strong$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))
dnds.repr.weak <- dnds %>% filter(dataset == repr_sim, del == "weak")
dnds.repr.weak$method <- factor(dnds.repr.weak$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))
dnds %>% group_by(dataset, del, method) %>% do(rraw = cor(.$true, .$dnds), braw = glm(dnds ~ offset(true), dat=.)) %>% mutate(r = rraw[1], r2 = r^2, b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) %>% ungroup() -> sim.dat.r.b
sim.dat.r.b$method <- factor(sim.dat.r.b$method, levels = c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes"), labels = c("No penalty", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel"))

## Figure 2
fig2 <- ggplot(dnds.repr.strong, aes(x = true, y = dnds)) + geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("Predicted dN/dS") + scale_y_continuous(limits=c(0,0.8)) + scale_x_continuous(limits=c(0,0.8)) + facet_grid(~method) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))


## Figure 3A
sim.dat.r.b %>% filter(del == "strong") %>% ggplot(aes(x = method, y = r2)) + geom_boxplot() + xlab("Inference Method") + ylab("Variance Explained") +  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 16))-> boxplot.sim.r2
## Figure 3B
sim.dat.r.b %>% filter(del == "strong") %>% ggplot(aes(x = method, y = b)) + geom_boxplot() + xlab("Inference Method") + ylab("Estimator Bias") + geom_hline(yintercept=0 ) +  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 16)) -> boxplot.sim.b
## Full Figure 3
fig3 <- plot_grid(boxplot.sim.r2, boxplot.sim.b, nrow=2, labels=c("A", "B"), scale = 0.95)

# #code if we want to combine figs2,3
# ggdraw() +
#   draw_plot(scatter.dnds.repr.strong, 0, 0.53, 1, 0.48) +
#   draw_plot(boxplot.sim.r2, 0, 0, 0.48, 0.48) +
#   draw_plot(boxplot.sim.b, 0.52, 0, 0.48, 0.48) +
#   draw_plot_label(c("A", "B", "C"), c(0, 0, 0.52), c(1, 0.5, 0.5), size = 17) ->p
# 










































####################### NOOOOOOOOOOOO!!!!!!!!!!!!!! ########################






sc.strong.raw <- read_csv(paste0("../../results/summarized_results/1B4T_A_delstrong_selection_coefficients.csv"))
sc.strong.raw %>% filter(dataset == "1B4T_A", del == "strong") %>% select(-realcoeff) %>% spread(method, binnedcoeff) %>% gather(method, binnedcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, binnedcoeff) -> sc.strong.binned
sc.strong.raw %>% filter(dataset == "1B4T_A", del == "strong") %>% select(-binnedcoeff) %>% spread(method, realcoeff) %>% gather(method, realcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, realcoeff) -> sc.strong.real

sc.weak.raw <- read_csv(paste0("../../results/summarized_results/1B4T_A_delweak_selection_coefficients.csv"))
sc.weak.raw %>% filter(dataset == "1B4T_A", del == "weak") %>% select(-realcoeff) %>% spread(method, binnedcoeff) %>% gather(method, binnedcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, binnedcoeff) -> sc.weak.binned
sc.weak.raw %>% filter(dataset == "1B4T_A", del == "weak") %>% select(-binnedcoeff) %>% spread(method, realcoeff) %>% gather(method, realcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, realcoeff) -> sc.weak.real

sc.strong.binned$method <- factor(sc.strong.binned$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
sc.weak.binned$method <- factor(sc.weak.binned$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
sc.strong.real$method <- factor(sc.strong.real$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
sc.weak.real$method <- factor(sc.weak.real$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))


plot.dnds.strong <- ggplot(dnds.strong, aes(x = true, y = dnds)) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(~method)
plot.dnds.weak <- ggplot(dnds.weak, aes(x = true, y = dnds)) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(~method)

plot.sc.strong.binned <- ggplot(sc.strong.binned, aes(x = binnedcoeff)) + geom_density(alpha = 0.8,fill="blue", color="blue",aes(x=true)) + geom_density(alpha = 0.7, fill = "orange", color="orange") + facet_grid(~method)
plot.sc.weak.binned <- ggplot(sc.weak.binned, aes(x = binnedcoeff)) + geom_density(fill="blue", color="blue", aes(x=true)) + geom_density(alpha = 0.7, fill = "orange", color="orange") + facet_grid(~method)


strong.grid <- plot_grid(plot.dnds.strong, plot.sc.strong.binned, nrow=2)
weak.grid   <- plot_grid(plot.dnds.weak, plot.sc.weak.binned, nrow=2)





# 
# methods <- c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01")
# for (name in sim.datasets){
#     for (del in c("weak", "strong")){
#         fullname <- paste0(name, "_del",del)
#         sc <- read_csv(paste0("../../results/summarized_results/", name, "_selection_coefficients.csv"))
#         sc$method <- factor(sc$method, levels=c("true", "nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# 
#         a <- dnds %>% filter(dataset == name) %>% ggplot(aes(x = true, y = dnds)) + facet_grid(~method) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0)
#         b <- sc %>% ggplot(aes(x=binnedcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
#         p <- plot_grid(a,b,nrow=2)
#         print(p)
#         readline()
#     }
# }    
# 
# 
# 
# 
# dnds_raw <- read.csv("../../results/summarized_results/empirical_derived_dnds.csv")
# dnds_raw %>% spread(method,dnds) -> dnds
# 
# for (name in emp.datasets){
#     print(name)
#     sc <- read_csv(paste0("../../results/", name, "_selection_coefficients.csv"))
#     sc$method <- factor(sc$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# 
#     a <- dnds %>% filter(dataset == name) %>% ggplot(aes(x = nopenal, y = phylobayes)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0)
#     b <- sc %>% filter(method %in% c("nopenal", "phylobayes")) %>% ggplot(aes(x=binnedcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
#     p <- plot_grid(a,b,nrow=2)
#     print(p)
#     readline()
# }    
# 
#     
# #     a <- dnds %>% filter(dataset == dataset) %>% ggplot(aes(x = nopenal, y = d0.01)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d0.01")
# #     b <- dnds %>% filter(dataset == dataset) %>%ggplot(aes(x = nopenal, y = d0.1)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d0.1")
# #     cc <- dnds %>% filter(dataset == dataset) %>%ggplot(aes(x = nopenal, y = d1.0)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d1.0")
# #     
# # 
# #     d <- sc %>% filter(method %in% c("nopenal", "d0.01")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
# #     e <- sc %>% filter(method %in% c("nopenal", "d0.1")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
# #     f <- sc %>% filter(method %in% c("nopenal", "d1.0")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
#     
#     plotgrid <- plot_grid(a,b,cc,d,e,f,nrow=2)
#     print(plotgrid)
#     readline()
# }
# 
# 
# 
# 
# # Some possibly useful code, REQUIRES PHYLOBAYES TO BE DONE! # 
# extract_corr_bias <- function(df, type)
# {
#     df %>% spread(method,dnds) %>% 
#            group_by(dataset) %>% 
#            do(rraw = cor(.$nopenal, .$phylobayes), braw = lm(phylobayes ~ offset(nopenal), data=.)) %>% 
#            mutate(r = rraw[1], b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) %>% arrange(b,r) -> newdf
#     newdf$type <- type
#     newdf
# }
# 
# 
# dnds.sim <- read.csv("../../results/summarized_results/simulation_derived_dnds.csv")
# dnds.sim %>% filter(del == "weak") %>% extract_corr_bias("weak") -> sim.weak.r.b
# dnds.sim %>% filter(del == "strong") %>% extract_corr_bias("strong") -> sim.strong.r.b
# dnds.emp <- read.csv("../../results/summarized_results/empirical_derived_dnds.csv")
# emp.r.b <- extract_corr_bias(dnds.emp, "empirical")
# full.r.b <- rbind(sim.weak.r.b, sim.strong.r.b, emp.r.b)
# full.r.b$type <- factor(full.r.b$type, levels=c("strong", "weak", "empirical"), labels=c("Simulated, strong", "Simulated, weak", "Empirical"))
# rplot_swpb <- ggplot(full.r.b, aes(x = type, y = r, group = type, fill=type)) + geom_boxplot() + scale_fill_manual(values=c("orangered2", "royalblue3")) + xlab("Data type") + ylab("swMutSel-pbMutSel Correlation")
# bplot_swpb <- ggplot(full.r.b, aes(x = type, y = b, group = type, fill=type)) + geom_boxplot() + scale_fill_manual(values=c("orangered2", "royalblue3")) + xlab("Data type") + ylab("Average swMutSel > pbMutSel")
# rb <- plot_grid(rplot, bplot, nrow=1)
# 
