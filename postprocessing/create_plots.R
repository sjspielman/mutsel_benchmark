require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(grid)


repr_sim <- "1IBS_A" # Selected because this dataset has the most codon columns (291)

result_directory <- "../results/summarized_results/"
true_directory   <- "../scripts/simulation/true_simulation_parameters/"
maintext_plot_directory   <- "plots/maintext/"
si_plot_directory   <- "plots/SI/"

jsd.results <- read.csv(paste0(result_directory, "simulation_jsd.csv"))
sim.dnds <- read.csv(paste0(result_directory, "simulation_derived_dnds.csv"))
sim.dnds %>% spread(method,dnds) %>% gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, site, true, dnds, method)-> spread.dnds
methods_levels <- c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")
methods_labels <- c("Unpenalized", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel")
sim.datasets <- unique(sim.dnds$dataset)


##########################################################################################
################################ Main text plots #########################################
##########################################################################################



##### selection coefficients across methods ####
theme_set(theme_cowplot() + theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=10), axis.title = element_text(size=13), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 11), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))
strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
true.sc <- read_csv(paste0(true_directory, repr_sim, "_delstrong_true_selcoeffs.csv"))
strong.sc.full %>% filter(method != "true") -> strong.sc.full
strong.sc.full$true.binned <- true.sc$binnedcoeff
strong.sc.full$true.real <- true.sc$realcoeff
strong.sc.full$method <- factor(strong.sc.full$method, levels = methods_levels, labels = methods_labels)

figsc <- ggplot(strong.sc.full, aes(x = binnedcoeff)) + geom_density(aes(x = true.binned), fill = "grey40") + geom_density(fill = "yellow", alpha = 0.4) + facet_grid(~method) + xlab("Scaled Selection Coefficients") + ylab("Density")
save_plot(paste0(maintext_plot_directory, "sc_across_methods_raw.pdf"), figsc, base_width=12, base_height=2.5)



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

save_plot(paste0(maintext_plot_directory, "jsd_maintext.pdf"), fig2, base_width = 7, base_height = 6) 

#### SI weak jsd boxplots
mean.jsd.weak.boxplots <- ggplot(jsd.weak.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
save_plot(paste0(maintext_plot_directory, "jsd_weak_SI.pdf"), mean.jsd.weak.boxplots, base_width = 6.5, base_height = 4) 




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
save_plot(paste0(maintext_plot_directory, "repr_dnds_scatter.pdf"), fig3, base_width=10, base_height=2)


## Figure 4A
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14)))
sim.dnds.r.b %>% ggplot(aes(x = method, y = r2)) + geom_boxplot() + xlab("Inference Method") + ylab("Variance Explained") -> boxplot.sim.r2
## Figure 4B
sim.dnds.r.b %>%  ggplot(aes(x = method, y = b)) + geom_boxplot() + xlab("Inference Method") + ylab("Estimator Bias") + geom_hline(yintercept=0 ) -> boxplot.sim.b
## Full Figure 4
fig4 <- plot_grid(boxplot.sim.r2, boxplot.sim.b, nrow=2, labels=c("A", "B"), label_size=17, scale=0.925)
save_plot(paste0(maintext_plot_directory, "r2_bias.pdf"), fig4, base_width=7, base_height=6)


#### dnds, jsd scatter and slope 
# function to return pvalue from an lm object
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

true.jsd <- sim.dnds %>% filter(method=="true") %>% mutate(truednds = dnds) %>% select(-method,-dnds) %>% left_join(jsd.results)
strong <- true.jsd %>% filter(del == "strong", dataset == repr_sim)
strong$method <- factor(strong$method, levels = methods_levels, labels = methods_labels)
fig5a <- ggplot(strong, aes(x = truednds, y = jsd)) + geom_point(size=1) + facet_grid(~method) + geom_smooth(method="lm", color="red") + scale_y_continuous(limits=c(0, 0.42)) + xlab("True dN/dS") + ylab("Site JSD") + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10), axis.title = element_text(size = 13), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=10))

alpha <- 0.01/length(methods_levels) #  Bonferroni correction.
true.jsd %>% filter(del == "strong") %>% group_by(dataset,method) %>% do(fit=lm(jsd~truednds, data=.)) %>% mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < alpha) %>% select(-fit) -> jsd.true.slope.p
jsd.true.slope.p$method <- factor(jsd.true.slope.p$method, levels = methods_levels, labels = methods_labels)
fig5b <- ggplot(jsd.true.slope.p, aes(x = method, y = slope, shape=sig)) + geom_point(position = position_jitter(w = 0.28)) + scale_shape_manual(values=c(1,19)) + geom_hline(yintercept=0) + theme(legend.position="none") + xlab("Inference Method") + ylab("Slope")+ theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 10), axis.title = element_text(size = 13))

fig5 <- plot_grid(fig5a, fig5b, nrow=2, labels=c("A", "B"), label_size=17)
save_plot(paste0(maintext_plot_directory, "jsd_dnds.pdf"), fig5, base_width=9, base_height=4)



## weak and strong dnds, selection coefficients
theme_set(theme_cowplot() + theme(axis.text = element_text(size=11), axis.title = element_text(size=12), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 10), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))

strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
weak.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delweak_selection_coefficients.csv"))
strong.sc <- filter(strong.sc.full, method %in% c("true", "nopenal", "phylobayes"))
weak.sc <- filter(weak.sc.full, method %in% c("true", "nopenal", "phylobayes"))
sc <- rbind(strong.sc, weak.sc)
sc$method <- factor(sc$method, levels=c("true", "nopenal", "phylobayes"))
spread.dnds.named <- spread.dnds
spread.dnds.named$del <- factor(spread.dnds.named$del, levels=c("strong", "weak"), labels = c("Strongly deleterious", "Weakly deleterious"))

fig6a <- spread.dnds.named %>% filter(dataset == repr_sim, method == "nopenal") %>% ggplot(aes(x = true, y = dnds)) + geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + facet_grid(~del) + xlab("True dN/dS") + ylab("swMutSel dN/dS") + scale_x_continuous(limits = c(0,0.8)) + scale_y_continuous(limits=c(0, 0.8))
fig6b <- sc %>% filter(method %in% c("true","nopenal")) %>% ggplot(aes(x = binnedcoeff, fill = method)) + geom_density() + facet_grid(~del) + scale_fill_manual(values =c("grey40", rgb(1, 1, 0, 0.4))) + ylab("Density") + xlab("Scaled Selection Coefficients") + theme(strip.text = element_blank())
fig6c <- spread.dnds.named %>% filter(dataset == repr_sim, method == "phylobayes") %>% ggplot(aes(x = true, y = dnds)) + geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + facet_grid(~del) + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + scale_x_continuous(limits = c(0,0.8)) + scale_y_continuous(limits=c(0, 0.8))
fig6d <- sc %>% filter(method %in% c("true","phylobayes")) %>% ggplot(aes(x = binnedcoeff, fill = method)) + geom_density() + facet_grid(~del) + scale_fill_manual(values = c("grey40", rgb(1, 1, 0, 0.4))) + ylab("Density") + xlab("Scaled Selection Coefficients") + theme(strip.text = element_blank())

fig6 <- plot_grid(fig6a, fig6c, fig6b, fig6d, nrow=2, labels=c("A", "B", "C", "D"), label_size = 15, scale=0.95)
save_plot(paste0(maintext_plot_directory, "dnds_sc_weakstrong_raw.pdf"), fig6, base_width=10, base_height=5)






##########################################################################################
################################### SI plots #############################################
##########################################################################################

# dN/dS scatterplots for other datasets
theme_set(theme_cowplot() + theme(plot.title = element_text(vjust = -8, hjust=1.08), plot.margin = unit(c(0.1,2.2,0.1,0.1),"cm"), axis.text = element_text(size = 11), axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 11), panel.border = element_rect(size = 0.5), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), strip.text = element_text(size=13)))
scatters_strong <- list()
scatters_weak <- list()
i <- 1
for (d in sim.datasets){
    if (d == repr_sim) next
    subdat <- spread.dnds %>% filter(dataset == d)
    subdat$method <- factor(subdat$method, levels = methods_levels, labels = methods_labels)

    pstrong <- subdat %>% filter(del == "strong") %>% ggplot(aes(x = true, y = dnds)) + geom_point(size=1.5) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("Predicted dN/dS") + scale_y_continuous(limits=c(0,0.85)) + scale_x_continuous(limits=c(0,0.85)) + facet_grid(~method) + ggtitle(d)
    pweak   <- subdat %>% filter(del == "weak") %>% ggplot(aes(x = true, y = dnds)) + geom_point(size=1.5) + geom_abline(slope = 1, intercept = 0, color="red") + xlab("True dN/dS") + ylab("Predicted dN/dS") + scale_y_continuous(limits=c(0,0.85)) + scale_x_continuous(limits=c(0,0.85)) + facet_grid(~method) + ggtitle(d)

    scatters_strong[[i]] <- pstrong
    scatters_weak[[i]] <- pweak
    i <- i+1
}
grid_strong <- plot_grid(plotlist = scatters_strong, nrow=10, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), scale=0.98)
save_plot(paste0(plot_directory, "strong_dnds_scatter_SI.pdf"), grid_strong, base_width = 14, base_height=20)
grid_weak <- plot_grid(plotlist = scatters_weak, nrow=10, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), scale=0.98)
save_plot(paste0(si_plot_directory, "weak_dnds_scatter_SI.pdf"), grid_weak, base_width = 14, base_height=20)


# selection coefficient distributions for all datasets across all methods
theme_set(theme_cowplot() + theme(plot.title = element_text(vjust = -10, hjust=1.085, size=8), plot.margin = unit(c(0.1,2.6,0.1,0.1),"cm"), axis.text = element_text(size=11), axis.title = element_text(size=12), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 11), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))
sc_strong <- list()
sc_weak <- list()
i <- 1
for (d in sim.datasets){
    print(d)
    strong.sc <- read_csv(paste0(result_directory, d, "_delstrong_selection_coefficients.csv"))
    strong.sc <- strong.sc %>% filter(method != "true")
    strong.sc$method <- factor(strong.sc$method, levels = methods_levels, labels = methods_labels)
    strong.sc.true <- read_csv(paste0(true_directory, d, "_delstrong_true_selcoeffs.csv"))

    weak.sc <- read_csv(paste0(result_directory, d, "_delweak_selection_coefficients.csv"))
    weak.sc <- weak.sc %>% filter(method != "true")
    weak.sc$method <- factor(weak.sc$method, levels = methods_levels, labels = methods_labels)
    weak.sc.true <- read_csv(paste0(true_directory, d, "_delweak_true_selcoeffs.csv"))
    
    pstrong <- ggplot(NULL) + geom_density(data = strong.sc.true, aes(x = binnedcoeff), fill = "grey40") + geom_density(data = strong.sc, aes(x = binnedcoeff), fill = "yellow", alpha = 0.4) + facet_grid(~method) + xlab("Scaled Selection Coefficients") + ylab("Density") + ggtitle(d)
    pweak <- ggplot(NULL) + geom_density(data = weak.sc.true, aes(x = binnedcoeff), fill = "grey40") + geom_density(data = weak.sc, aes(x = binnedcoeff), fill = "yellow", alpha = 0.4) + facet_grid(~method) + xlab("Scaled Selection Coefficients") + ylab("Density") + ggtitle(d) 

    sc_strong[[i]] <- pstrong
    sc_weak[[i]] <- pweak
    i <- i+1
}
grid_strong <- plot_grid(plotlist = sc_strong, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))
save_plot(paste0(si_plot_directory, "selcoeffs_strong_SI.pdf"), grid_strong, base_width = 12, base_height=18)
grid_weak <- plot_grid(plotlist = sc_weak, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))
save_plot(paste0(si_plot_directory, "selcoeffs_weak_SI.pdf"), grid_weak, base_width = 12, base_height=18)




# r2, bias boxplots for weakly deleterious datasets
theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14)))

spread.dnds %>% filter(del == "weak") %>% group_by(dataset, method) %>% do(rraw = cor(.$true, .$dnds), braw = glm(dnds ~ offset(true), dat=.)) %>% mutate(r = rraw[1], r2 = r^2, b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) %>% ungroup() -> sim.dnds.r.b
sim.dnds.r.b$method <- factor(sim.dnds.r.b$method, levels = methods_levels, labels = methods_labels)
sim.dnds.r.b %>% ggplot(aes(x = method, y = r2)) + geom_boxplot() + xlab("Inference Method") + ylab("Variance Explained") -> boxplot.sim.r2
sim.dnds.r.b %>%  ggplot(aes(x = method, y = b)) + geom_boxplot() + xlab("Inference Method") + ylab("Estimator Bias") + geom_hline(yintercept=0 ) -> boxplot.sim.b
r2.b.weak <- plot_grid(boxplot.sim.r2, boxplot.sim.b, nrow=2, labels=c("A", "B"), label_size=17, scale=0.925)
save_plot(paste0(si_plot_directory, "weak_r2_bias_SI.pdf"), r2.b.weak, base_width=7, base_height=6)




# JSD boxplots for weakly deleterious datasets
theme_set(theme_cowplot() + theme(axis.text = element_text(size=11), axis.title = element_text(size=12), panel.border = element_rect(size = 0.5)))

jsd.results %>% filter(del == "weak") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.weak.summary
jsd.weak.summary$method <- factor(jsd.weak.summary$method, levels = methods_levels, labels = methods_labels)
mean.jsd.weak.boxplots <- ggplot(jsd.weak.summary, aes(x = method, y = meanjsd)) + geom_boxplot() + xlab("Inference Method") + ylab("Average JSD") + scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))
save_plot(paste0(si_plot_directory, "jsd_weak_SI.pdf"), mean.jsd.weak.boxplots, base_width = 7, base_height = 4) 



# Show the behavior of fitness values for a representative dataset
theme_set(theme_cowplot() + theme(axis.text = element_text(size=11), axis.title = element_text(size=12), panel.border = element_rect(size = 0.5), strip.text = element_text(size = 10), panel.margin = unit(0.75, "lines"), strip.background = element_rect(fill="white"), legend.position="none"))

dat <- read.csv(paste0(repr_sim, "_delweak_rescaled_fitness_across_methods.txt"))
dat$truefit <- rep(dat$fitness[dat$method == "true"], 9)
dat %>% filter(method!="true") -> dat
dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
weak.fit.scatter <- ggplot(dat, aes(x = truefit, y = fitness))+geom_point(alpha=0.3) +facet_grid(~method)+geom_abline(slope=1, intercept=0, color="red", size=1) + xlab("True fitness") + ylab("Inferred fitness")
save_plot(paste0(si_plot_directory, "weak_fitness_scatter_SI.pdf"), weak.fit.scatter, base_width = 8, base_height=3)


dat <- read.csv(paste0(repr_sim, "_delstrong_rescaled_fitness_across_methods.txt"))
dat$truefit <- rep(dat$fitness[dat$method == "true"], 9)
dat %>% filter(method!="true") -> dat
dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
strong.fit.scatter <- ggplot(dat, aes(x = truefit, y = fitness))+geom_point(alpha=0.3) +facet_grid(~method)+geom_abline(slope=1, intercept=0, color="red", size=1) + xlab("True fitness") + ylab("Inferred fitness")
save_plot(paste0(si_plot_directory, "strong_fitness_scatter_SI.pdf"), strong.fit.scatter, base_width = 8, base_height=3)



