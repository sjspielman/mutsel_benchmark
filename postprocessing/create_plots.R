require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(grid)


repr_sim <- "1IBS_A" # Selected because this dataset has the most codon columns (291)

result_directory <- "../results/summarized_results/"
true_directory   <- "../scripts/simulation/true_simulation_parameters/"
plot_directory   <- "plots/"

jsd.results <- read.csv(paste0(result_directory, "simulation_jsd.csv"))
sim.dnds <- read.csv(paste0(result_directory, "simulation_derived_dnds.csv"))
sim.dnds %>% spread(method,dnds) %>% gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, site, true, dnds, method)-> spread.dnds
methods_levels <- c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")
methods_labels <- c("Unpenalized", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel")

#emp.dnds <- read.csv(paste0(result_directory, "empirical_derived_dnds.csv"))

theme_set(theme_cowplot() + theme(panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12)))




##### selection coefficients across methods ####
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=10), axis.title = element_text(size=13), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 11), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))

true.sc <- read_csv(paste0(true_directory, repr_sim, "_delstrong_true_selcoeffs.csv"))
strong.sc.full %>% filter(method != "true") -> strong.sc.full
strong.sc.full$true.binned <- true.sc$binnedcoeff
strong.sc.full$true.real <- true.sc$realcoeff
strong.sc.full$method <- factor(strong.sc.full$method, levels = methods_levels, labels = methods_labels)

figsc <- ggplot(strong.sc.full, aes(x = binnedcoeffF)) + geom_density(aes(x = true.binned), fill = "grey40") + geom_density(fill = "yellow", alpha = 0.4) + facet_grid(~method) + xlab("Scaled Selection Coefficients") + ylab("Density")
ggsave(paste0(plot_directory, "sc_across_methods_raw.pdf"), figsc, width=12, height=2.5)



#### Boxplots of Jensen-Shannon distance among methods for simulated datasets ####
jsd.results %>% filter(del == "strong") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.strong.summary
jsd.strong.summary$method <- factor(jsd.strong.summary$method, levels = methods_levels, labels = methods_labels)
jsd.results %>% filter(del == "weak") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.weak.summary
jsd.weak.summary$method <- factor(jsd.weak.summary$method, levels = methods_levels, labels = methods_labels)
jsd.results %>% filter(del == "strong", dataset == repr_sim) -> jsd.repr.sim
jsd.repr.sim$method <- factor(jsd.repr.sim$method, levels = methods_levels, labels = methods_labels)

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)))
## Figure 2A 
repr.jsd.boxplots <- ggplot(jsd.repr.sim, aes(x = method, y = jsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Site JSD") + scale_y_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) #+ ggtitle("JSD for representative dataset")
## Figure 2B
mean.jsd.strong.boxplots <- ggplot(jsd.strong.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) #+ ggtitle("Mean JSD across datasets")
## Full Figure 2
fig2 <- plot_grid(repr.jsd.boxplots, mean.jsd.strong.boxplots, nrow=2, labels=c("A", "B"), label_size = 17, scale=0.925)

ggsave(paste0(plot_directory, "jsd_maintext.pdf"), fig2, width = 7, height = 6) 

#### SI weak jsd boxplots
mean.jsd.weak.boxplots <- ggplot(jsd.weak.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
ggsave(paste0(plot_directory, "jsd_weak_SI.pdf"), mean.jsd.weak.boxplots, width = 6.5, height = 4) 




##### Figures 3,4: dN/dS scatterplots and boxplots for simulated datasets #####
dnds.repr.strong <- spread.dnds %>% filter(dataset == repr_sim, del == "strong")
dnds.repr.strong$method <- factor(dnds.repr.strong$method, levels = methods_levels, labels = methods_labels)
dnds.repr.weak <- spread.dnds %>% filter(dataset == repr_sim, del == "weak")
dnds.repr.weak$method <- factor(dnds.repr.weak$method, levels = methods_levels, labels = methods_labels)

spread.dnds %>% filter(del == "strong") %>% group_by(dataset, method) %>% do(rraw = cor(.$true, .$dnds), braw = glm(dnds ~ offset(true), dat=.)) %>% mutate(r = rraw[1], r2 = r^2, b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) %>% ungroup() -> sim.dnds.r.b
sim.dnds.r.b$method <- factor(sim.dnds.r.b$method, levels = methods_levels, labels = methods_labels)

## Figure 3
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10), axis.title = element_text(size = 12), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=12)))
fig3 <- ggplot(dnds.repr.strong, aes(x = true, y = dnds)) + geom_point(size=1) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("Predicted dN/dS") + scale_y_continuous(limits=c(0,0.85)) + scale_x_continuous(limits=c(0,0.85)) + facet_grid(~method) 
ggsave(paste0(plot_directory, "repr_dnds_scatter.pdf"), fig3, width=10, height=2)


## Figure 4A
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14)))
sim.dnds.r.b %>% ggplot(aes(x = method, y = r2)) + geom_boxplot() + xlab("Inference Method") + ylab("Variance Explained") -> boxplot.sim.r2
## Figure 4B
sim.dnds.r.b %>%  ggplot(aes(x = method, y = b)) + geom_boxplot() + xlab("Inference Method") + ylab("Estimator Bias") + geom_hline(yintercept=0 ) -> boxplot.sim.b
## Full Figure 4
fig4 <- plot_grid(boxplot.sim.r2, boxplot.sim.b, nrow=2, labels=c("A", "B"), label_size=17, scale=0.925)
ggsave(paste0(plot_directory, "r2_bias.pdf"), fig4, width=7, height=6)

# #code if we want to combine figs2,3
# ggdraw() +
#   draw_plot(scatter.dnds.repr.strong, 0, 0.53, 1, 0.48) +
#   draw_plot(boxplot.sim.r2, 0, 0, 0.48, 0.48) +
#   draw_plot(boxplot.sim.b, 0.52, 0, 0.48, 0.48) +
#   draw_plot_label(c("A", "B", "C"), c(0, 0, 0.52), c(1, 0.5, 0.5), size = 17) ->p
# 




## Figure 5

# function to return pvalue from an lm object
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


## Figure 5A
true.jsd <- sim.dnds %>% filter(method=="true") %>% mutate(truednds = dnds) %>% select(-method,-dnds) %>% left_join(jsd.results)
strong <- true.jsd %>% filter(del == "strong", dataset == repr_sim)
strong$method <- factor(strong$method, levels = methods_levels, labels = methods_labels)
fig5a <- ggplot(strong, aes(x = truednds, y = jsd)) + geom_point(size=1) + facet_grid(~method) + geom_smooth(method="lm", color="red") + scale_y_continuous(limits=c(0, 0.42)) + xlab("True dN/dS") + ylab("Site JSD") + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10), axis.title = element_text(size = 13), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=10))


## Figure 5B
alpha <- 0.01/length(methods_levels) #  Bonferroni correction.
true.jsd %>% filter(del == "strong") %>% group_by(dataset,method) %>% do(fit=lm(jsd~truednds, data=.)) %>% mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < alpha) %>% select(-fit) -> jsd.true.slope.p
jsd.true.slope.p$method <- factor(jsd.true.slope.p$method, levels = methods_levels, labels = methods_labels)
fig5b <- ggplot(jsd.true.slope.p, aes(x = method, y = slope, shape=sig)) + geom_point(position = position_jitter(w = 0.28)) + scale_shape_manual(values=c(1,19)) + geom_hline(yintercept=0) + theme(legend.position="none") + xlab("Inference Method") + ylab("Slope")+ theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10), axis.title = element_text(size = 13))

## Full figure 5
fig5 <- plot_grid(fig5a, fig5b, nrow=2, labels=c("A", "B"), label_size=17)
ggsave(paste0(plot_directory, "jsd_dnds.pdf"), fig5, width=9, height=4)



## Figure 5: weak and strong dnds, selection coefficients
theme_set(theme_cowplot() + theme(axis.text = element_text(size=11), axis.title = element_text(size=12), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 10), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))

strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
weak.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delweak_selection_coefficients.csv"))
strong.sc <- filter(strong.sc.full, method %in% c("true", "nopenal", "phylobayes"))
weak.sc <- filter(weak.sc.full, method %in% c("true", "nopenal", "phylobayes"))
sc <- rbind(strong.sc, weak.sc)
sc$method <- factor(sc$method, levels=c("true", "nopenal", "phylobayes"))
spread.dnds$del <- factor(spread.dnds$del, levels=c("strong", "weak"), labels = c("Highly deleterious", "Weakly deleterious"))


fig6a <- spread.dnds %>% filter(dataset == repr_sim, method == "nopenal") %>% ggplot(aes(x = true, y = dnds)) + geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + facet_grid(~del) + xlab("True dN/dS") + ylab("swMutSel dN/dS") + scale_x_continuous(limits = c(0,0.8)) + scale_y_continuous(limits=c(0, 0.8))
fig6b <- sc %>% filter(method %in% c("true","nopenal")) %>% ggplot(aes(x = binnedcoeff, fill = method)) + geom_density() + facet_grid(~del) + scale_fill_manual(values =c("grey40", rgb(1, 1, 0, 0.4))) + ylab("Density") + xlab("Scaled Selection Coefficients") + theme(strip.text = element_blank())
fig6c <- spread.dnds %>% filter(dataset == repr_sim, method == "phylobayes") %>% ggplot(aes(x = true, y = dnds)) + geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + facet_grid(~del) + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + scale_x_continuous(limits = c(0,0.8)) + scale_y_continuous(limits=c(0, 0.8))
fig6d <- sc %>% filter(method %in% c("true","phylobayes")) %>% ggplot(aes(x = binnedcoeff, fill = method)) + geom_density() + facet_grid(~del) + scale_fill_manual(values = c("grey40", rgb(1, 1, 0, 0.4))) + ylab("Density") + xlab("Scaled Selection Coefficients") + theme(strip.text = element_blank())

## Figure6
fig6 <- plot_grid(fig6a, fig6b, fig6c, fig6d, nrow=1, labels=c("A", "B", "C", "D"), label_size = 15, scale=0.95)
ggsave(paste0(plot_directory, "dnds_sc_weakstrong_raw.pdf"), fig6, width=10, height=5)














############################################################################################
############################################################################################
############################################################################################
############################################################################################
################################# NOOOOOOOOOOOO!!!!!!!!!!!!!! ##############################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################



# 
# ## Figure 6
# # Something to do with empirical datasets.
# # We can plot pbmutsel against nopenalty for a simulated and for an empirical. Then show violin plots for r2 and bias between pbmutsel and nopenalty, for both simulation and empirical
# sim.dnds %>% filter(del == "strong") %>% select(-del) %>% spread(method,dnds) %>% select(dataset,site,nopenal,phylobayes) %>% mutate(type = "Simulation")  -> sim
# emp.dnds %>% unique() %>% filter(method %in% c("nopenal", "phylobayes")) %>% spread(method,dnds) %>% select(dataset, site, nopenal, phylobayes) %>% mutate(type = "Empirical") %>% rbind(sim)-> sim.emp.dnds
# sim.emp.dnds %>% group_by(dataset, type) %>% do(rraw = cor(.$nopenal, .$phylobayes), braw = glm(nopenal ~ offset(phylobayes), dat=.)) %>% mutate(r = rraw[1], r2 = r^2, b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) -> hmm
# ggplot(hmm, aes(x = type, y = r)) + geom_boxplot() + xlab("Dataset") + ylab("Pearson Correlation")  -> r #NS
# ggplot(hmm, aes(x = type, y = b)) + geom_boxplot()  + xlab("Dataset") + ylab("swMutSel > pbMutsel") -> b     # P=2.167e-05,   mean of x : -0.08040115 mean of y : -0.10315696
#  
# rb <- plot_grid(r, b, nrow=1)
# ggsave(paste0(plot_directory, "figure6.pdf"), rb, width=7, height=4)
# 
# 
# 




# 
# 
# sc.strong.raw <- read_csv(paste0("../../results/summarized_results/1B4T_A_delstrong_selection_coefficients.csv"))
# sc.strong.raw %>% filter(dataset == "1B4T_A", del == "strong") %>% select(-realcoeff) %>% spread(method, binnedcoeff) %>% gather(method, binnedcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, binnedcoeff) -> sc.strong.binned
# sc.strong.raw %>% filter(dataset == "1B4T_A", del == "strong") %>% select(-binnedcoeff) %>% spread(method, realcoeff) %>% gather(method, realcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, realcoeff) -> sc.strong.real
# 
# sc.weak.raw <- read_csv(paste0("../../results/summarized_results/1B4T_A_delweak_selection_coefficients.csv"))
# sc.weak.raw %>% filter(dataset == "1B4T_A", del == "weak") %>% select(-realcoeff) %>% spread(method, binnedcoeff) %>% gather(method, binnedcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, binnedcoeff) -> sc.weak.binned
# sc.weak.raw %>% filter(dataset == "1B4T_A", del == "weak") %>% select(-binnedcoeff) %>% spread(method, realcoeff) %>% gather(method, realcoeff, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, dummy, true, method, realcoeff) -> sc.weak.real
# 
# sc.strong.binned$method <- factor(sc.strong.binned$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# sc.weak.binned$method <- factor(sc.weak.binned$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# sc.strong.real$method <- factor(sc.strong.real$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# sc.weak.real$method <- factor(sc.weak.real$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
# 
# 
# plot.dnds.strong <- ggplot(dnds.strong, aes(x = true, y = dnds)) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(~method)
# plot.dnds.weak <- ggplot(dnds.weak, aes(x = true, y = dnds)) + geom_point() + geom_abline(slope=1, intercept=0) + facet_grid(~method)
# 
# plot.sc.strong.binned <- ggplot(sc.strong.binned, aes(x = binnedcoeff)) + geom_density(alpha = 0.8,fill="blue", color="blue",aes(x=true)) + geom_density(alpha = 0.7, fill = "orange", color="orange") + facet_grid(~method)
# plot.sc.weak.binned <- ggplot(sc.weak.binned, aes(x = binnedcoeff)) + geom_density(fill="blue", color="blue", aes(x=true)) + geom_density(alpha = 0.7, fill = "orange", color="orange") + facet_grid(~method)
# 
# 
# strong.grid <- plot_grid(plot.dnds.strong, plot.sc.strong.binned, nrow=2)
# weak.grid   <- plot_grid(plot.dnds.weak, plot.sc.weak.binned, nrow=2)
# 
# 
# 
# 
# 
# 
# methods <- c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01")
# for (name in sim.dndsasets){
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
