# SJS
# Creates MS figures (output to either figures/main_text or figures/SI, depending on where figure is located in MS)

require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(grid)


theme_set(theme_cowplot() + theme(panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.25, "cm"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=14)))


result_directory <- "dataframes/"
true_directory   <- "../simulation/true_simulation_parameters/"
maintext_plot_directory   <- "NEWfigures/"
si_plot_directory   <- "NEWfigures/"

yeast <- read.csv(paste0(result_directory, "yeast_results.csv"))
dms <- read.csv(paste0(result_directory, "dms_results.csv"))

methods_levels <- c("nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes") # "mvn1","d1.0", 
methods_labels <- c("Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel")
yeast <- yeast %>% filter(method %in% methods_levels)
yeast$method <- factor(yeast$method, levels = methods_levels, labels = methods_labels)
yeast$bl <- factor(yeast$bl, levels = c(0.01, 0.5))
dms <- dms %>% filter(method %in% methods_levels)
dms$method <- factor(dms$method, levels = methods_levels, labels = methods_labels)
dms$bl <- factor(dms$bl, levels = c(0.01, 0.5))

primary.data <- yeast %>% filter(del == "strong") %>% select(-del) %>% left_join(dms)
bl_colors <- c("red", "blue")


alpha <- 0.01 # Significance
corrected.alpha <- alpha/length(methods_levels) #Bonferroni significance
repr_sim <- "1R6M_A"
datasets <- c(as.character(unique(yeast$dataset)), "NP", "HA")

summarize_dnds_entropy <- function(dat){
  dat %>% do(rraw.dnds = cor(.$true.dnds, .$dnds), 
     braw.dnds = glm(dnds ~ offset(true.dnds), dat=.),
     rraw.h = cor(.$true.entropy, .$entropy), 
     braw.h = glm(entropy ~ offset(true.entropy), dat=.)) %>% 
    mutate(r2.dnds = rraw.dnds[[1]]^2, 
           b.dnds = summary(braw.dnds)$coeff[1], 
           b.dnds.pvalue = summary(braw.dnds)$coeff[4],
           b.dnds.sig = b.dnds.pvalue<corrected.alpha,
           r2.entropy = rraw.h[[1]]^2, 
           b.entropy = summary(braw.h)$coeff[1], 
           b.entropy.pvalue = summary(braw.h)$coeff[4],
           b.entropy.sig = b.entropy.pvalue<corrected.alpha) %>% 
    select(-rraw.dnds, -rraw.h, -braw.dnds, -braw.h, -b.dnds.pvalue, -b.entropy.pvalue) -> dat.sum
  
  dat.sum
}

yeast.dnds.entropy.stats <- yeast %>% group_by(dataset, del, bl, method) %>% summarize_dnds_entropy()
dms.dnds.entropy.stats <- dms %>% group_by(dataset, bl, method) %>% summarize_dnds_entropy()
primary.dnds.entropy.stats <- primary.data %>% group_by(dataset, bl, method) %>% summarize_dnds_entropy()

# function to return pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



##########################################################################################
################################ Main text plots #########################################
##########################################################################################
print("Creating main text plots")



##################################################################################################################################
############### Figure 2: Boxplots of JSD and absolute sum differences for representative dataset, all datasets #################
##################################################################################################################################
print("Figure 2")

yeast %>% filter(del == "strong") %>% left_join(dms) %>% group_by(dataset, method, bl) %>% summarize(meanjsd = mean(jsd), meandiffsum = mean(diffsum)) -> jsd.diffsum


bl.legend.raw <- jsd.diffsum %>% ggplot(aes(x = method, y = meanjsd, fill = bl)) + 
  geom_boxplot() +
  scale_fill_manual(values = bl_colors, name = "Branch length") + 
  theme(legend.position = "bottom", legend.title = element_text(size=13), legend.text = element_text(size=12))
grobs <-  ggplotGrob(bl.legend.raw)$grobs
bl.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


yeast %>% filter(dataset == repr_sim, del == "strong") %>%
  ggplot(aes(x = method, y = jsd, fill = as.factor(bl))) + 
  geom_boxplot(position = position_dodge(0.9), width=0.5) +
  xlab("Inference Method") + ylab("Site JSD") + 
  theme(legend.position = "none") -> box.a

yeast %>% filter(dataset == repr_sim, del == "strong", bl == 0.01) %>%
  ggplot(aes(x = method, y = jsd)) + 
  geom_boxplot(position = position_dodge(0.9)) +
  xlab("Inference Method") + ylab("Site JSD") + 
  theme(legend.position = "none") -> box.b

jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = meanjsd)) + 
  geom_jitter(width = 0.5) + 
  xlab("Inference Method") + ylab("Average JSD") +
  scale_y_continuous(limits=c(0.13,0.5))  + 
  theme(legend.position = "none") -> box.c

jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = method, y = meanjsd)) + 
  geom_jitter(width = 0.5) + 
  xlab("Inference Method") + ylab("Average JSD") +
  scale_y_continuous(limits=c(0.13,0.5))  + 
  theme(legend.position = "none") -> box.d


jsd.boxplots <- plot_grid(box.a, box.b, box.c, box.d, nrow=2, labels=c("A", "B", "C", "D"))


jsd.diffsum %>%
  ggplot(aes(x = method, y = meanjsd)) + 
  geom_jitter(width = 0.5, size=2) + facet_grid(~bl)+
  xlab("Inference Method") + ylab("Average JSD") +
  scale_y_continuous(limits=c(0.13,0.5))  + 
  theme(legend.position = "none") -> mean.jsd.boxplots

jsd.diffsum %>%
  ggplot(aes(x = method, y = meanjsd, fill = as.factor(bl))) + 
  geom_boxplot(position = position_dodge(0)) +
  scale_fill_manual(values = bl_colors, name = "Branch length") + 
  xlab("Inference Method") + ylab("Average JSD") +
  scale_y_continuous(limits=c(0.13,0.5))  + 
  theme(legend.position = "none") -> mean.jsd.boxplots



fig2 <- plot_grid(repr.jsd.boxplots, mean.jsd.boxplots, bl.legend, nrow=3, labels=c("A", "B"), rel_heights=c(1,1,0.13))
save_plot(paste0(maintext_plot_directory, "jsd_boxplots.pdf"), fig2, base_width = 8, base_height = 5) 

yeast %>% filter(dataset == repr_sim, del == "strong") %>%
  ggplot(aes(x = method, y = diffsum, fill = as.factor(bl))) + 
  geom_boxplot(position = position_dodge(0.9)) +
  scale_fill_manual(values = bl_colors, name = "Branch length") + 
  xlab("Inference Method") + ylab("Site Difference")  + 
  theme(legend.position = "none") -> repr.diffsum.boxplots

yeast.jsd.diffsum %>% filter(del == "strong") %>%
  ggplot(aes(x = method, y = meandiffsum, fill = as.factor(bl))) + 
  geom_boxplot(position = position_dodge(0.9)) + 
  scale_fill_manual(values = bl_colors, name = "Branch length") + 
  xlab("Inference Method") + ylab("Average Difference")  + 
  theme(legend.position = "none") -> mean.diffsum.boxplots


fig.diffsum <- plot_grid(repr.diffsum.boxplots, mean.diffsum.boxplots, bl.legend, nrow=3, labels=c("A", "B"), rel_heights=c(1,1,0.13))
save_plot(paste0(maintext_plot_directory, "diffsum_boxplots.pdf"), fig.diffsum, base_width = 7, base_height = 5.5) 







##########################################################################################
##### Figure 3: True vs. predicted dN/dS: repr scatter, correlation and bias jitters #####
##########################################################################################
print("Figure 3")





repr.scatter.dnds <- yeast %>% filter(dataset == repr_sim, del == "strong") %>%
  ggplot(aes(x = true.dnds, y = dnds)) + 
  geom_point(size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
  xlab("True dN/dS") + ylab("Predicted dN/dS") + 
  scale_y_continuous(limits=c(0,0.88)) + scale_x_continuous(limits=c(0,0.88)) + 
  facet_grid(bl~method)
save_plot(paste0(maintext_plot_directory, "repr_scatter_dnds.pdf"),repr.scatter.dnds, base_width=11, base_height=4)






jitter.bl.legend.raw <- yeast.dnds.entropy.stats %>%
  ggplot(aes(x = method, y = r.dnds, color = as.factor(bl))) + 
  geom_jitter() + scale_color_manual(values=bl_colors, name = "Branch Length") + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12))
grobs <-  ggplotGrob(jitter.bl.legend.raw)$grobs
jitter.bl.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


jitter.r.dnds <- yeast.dnds.entropy.stats %>% filter(del == "strong") %>%
  ggplot(aes(x = method, y = r.dnds, color = as.factor(bl))) + 
  geom_jitter(width = 0.5, size=2) + scale_color_manual(values=bl_colors) + 
  xlab("Inference Method") + ylab("Pearson Correlation") + theme(legend.position = "none", axis.title.y = element_text(size=11.5))


jitter.b.dnds <- yeast.dnds.entropy.stats %>% filter(del == "strong") %>%
  ggplot(aes(x = method, y = b.dnds, color = as.factor(bl), shape = b.dnds.sig)) + 
  geom_jitter(width = 0.5, size=2) + 
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values=bl_colors) + 
  xlab("Inference Method") + ylab("Estimator Bias") + 
  geom_hline(yintercept=0 ) + theme(legend.position = "none", axis.title.y = element_text(size=13))

r.b.dnds <- plot_grid(jitter.r.dnds, jitter.b.dnds, jitter.bl.legend, nrow=3, labels=c("A", "B"), rel_heights=c(1,1,0.08))
save_plot(paste0(maintext_plot_directory, "r_bias_dnds.pdf"), r.b.dnds, base_width=8, base_height=5)


jitter.r.entropy <- yeast.dnds.entropy.stats %>% filter(del == "strong") %>%
  ggplot(aes(x = method, y = r.entropy, color = as.factor(bl))) + 
  geom_jitter(width = 0.5, size=2) + scale_color_manual(values=bl_colors) + 
  xlab("Inference Method") + ylab("Pearson Correlation") + theme(legend.position = "none", axis.title.y = element_text(size=11.5))
jitter.b.entropy <- yeast.dnds.entropy.stats %>% filter(del == "strong") %>%
  ggplot(aes(x = method, y = b.entropy, color = as.factor(bl), shape = b.entropy.sig)) + 
  geom_jitter(width = 0.5, size=2) + 
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values=bl_colors) + 
  xlab("Inference Method") + ylab("Estimator Bias") + 
  geom_hline(yintercept=0 ) + theme(legend.position = "none", axis.title.y = element_text(size=13))
r.b.entropy <- plot_grid(jitter.r.dnds, jitter.b.dnds, jitter.bl.legend, nrow=3, labels=c("A", "B"), rel_heights=c(1,1,0.08))
save_plot(paste0(si_plot_directory, "r_bias_entropy.pdf"), r.b.dnds, base_width=8, base_height=5)





###### RESULTS DIFFER BETWEEN BRANCH LENGTHS FOR THIS SECTION, SO THINK A BIT MORE ABOUT WHAT TO DO HERE #######
##### More specifically, pattern is the same but relationship to 0 is different. This is fully ok with me, because dN/dS is fairly useless at bl0.01.

##########################################################################################
##### Figure 4: JSD regressed on true dN/dS for representative, and slopes for all #######
##########################################################################################
print("Figure 4")

dnds.on.jsd.bl0.5 <- yeast %>% filter(del == "strong", dataset == repr_sim, bl == 0.5) %>%
  ggplot(aes(x = true.dnds, y = jsd)) + 
  geom_point(size=1) + facet_grid(~method) + 
  geom_smooth(method="lm", color="red") + 
  scale_y_continuous(limits=c(0, 0.6)) + 
  xlab("True dN/dS") + ylab("Site JSD")

#dnds.on.jsd.bl0.01 <- yeast %>% filter(del == "strong", dataset == repr_sim, bl == 0.01) %>%
#  ggplot(aes(x = true.dnds, y = jsd)) + 
#  geom_point(size=1) + facet_grid(~method) + 
#  geom_smooth(method="lm", color="red") + 
#  xlab("True dN/dS") + ylab("Site JSD")

jsd.true.slope.p <- yeast %>% filter(del == "strong") %>% 
  group_by(dataset,method,bl) %>% 
  do(fit=lm(jsd~true.dnds, data=.)) %>% 
  mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>% 
  select(-fit)

dnds.jsd.jitter.bl0.5 <- jsd.true.slope.p %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = slope, shape=sig)) + 
  geom_point(position = position_jitter(w = 0.4)) + 
  scale_shape_manual(values=c(1,19)) + 
  geom_hline(yintercept=0) + 
  theme(legend.position="none") + 
  xlab("Inference Method") + ylab("Slope")

#dnds.jsd.jitter.bl0.01 <- jsd.true.slope.p %>% filter(bl == 0.01) %>%
#  ggplot(aes(x = method, y = slope, shape=sig)) + 
#  geom_point(position = position_jitter(w = 0.4)) + 
#  scale_shape_manual(values=c(1,19)) + 
#  geom_hline(yintercept=0) + 
#  theme(legend.position="none") + 
#  xlab("Inference Method") + ylab("Slope")

fig4 <- plot_grid(dnds.on.jsd.bl0.5, dnds.jsd.jitter.bl0.5, nrow=2, labels=c("A", "B"))
save_plot(paste0(maintext_plot_directory, "jsd_dnds.pdf"), fig4, base_width=9.5, base_height=4)




##########################################################################################
##### Figure 5: weak vs. strong scatterplots and selection coefficient distributions #####
##########################################################################################
print("Figure 5")

strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
weak.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delweak_selection_coefficients.csv"))
strong.sc <- filter(strong.sc.full, method %in% c("true", "nopenal", "mvn1", "d0.01", "phylobayes"))
weak.sc <- filter(weak.sc.full, method %in% c("true", "nopenal", "mvn1", "d0.01", "phylobayes"))
sc <- rbind(strong.sc, weak.sc)
sc$method <- factor(sc$method, levels=c("true", "nopenal", "mvn1", "d0.01", "phylobayes"))
spread.dnds.named <- spread.dnds
spread.dnds.named$del <- factor(spread.dnds.named$del, levels=c("strong", "weak"), labels = c("Strongly deleterious", "Weakly deleterious"))

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.2, 0.1, 0.1, 0.2),"cm"), 
                                  panel.margin = unit(0.25, "cm"),
                                  axis.text = element_text(size=10), 
                                  axis.title = element_text(size=10.5, face = "bold"), 
                                  strip.background = element_rect(fill="white"),
                                  strip.text = element_text(size = 9)))

fig5a <- spread.dnds.named %>% filter(dataset == repr_sim, method == "nopenal") %>% 
  ggplot(aes(x = true, y = dnds)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
  facet_grid(~del) + xlab("True dN/dS") + ylab("Unpenalized dN/dS") +   
  scale_x_continuous(limits = c(0,0.9)) + 
  scale_y_continuous(limits = c(0,0.9))

fig5b <- spread.dnds.named %>% filter(dataset == repr_sim, method == "d0.01") %>% 
  ggplot(aes(x = true, y = dnds)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
  facet_grid(~del) + xlab("True dN/dS") + ylab("d0.01 dN/dS") + 
  scale_x_continuous(limits = c(0,0.9)) + 
  scale_y_continuous(limits = c(0,0.9))


fig5c <- spread.dnds.named %>% filter(dataset == repr_sim, method == "mvn1") %>% 
  ggplot(aes(x = true, y = dnds)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
  facet_grid(~del) + xlab("True dN/dS") + ylab("mvn1 dN/dS") + 
  scale_x_continuous(limits = c(0,0.9)) + 
  scale_y_continuous(limits = c(0,0.9))


fig5d <- spread.dnds.named %>% filter(dataset == repr_sim, method == "phylobayes") %>% 
  ggplot(aes(x = true, y = dnds)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
  facet_grid(~del) + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + 
  scale_x_continuous(limits = c(0,0.9)) + 
  scale_y_continuous(limits = c(0,0.9))


fig5e <- sc %>% filter(method %in% c("true","nopenal")) %>% 
  ggplot(aes(x = binnedcoeff, fill = method)) + 
  geom_density() + facet_grid(~del) + 
  scale_fill_manual(name = "", labels = c("True", "Unpenalized"), values =c("grey40", rgb(1, 1, 0, 0.4))) + 
  ylab("Density") + xlab("Scaled Selection Coefficients") + 
  theme(strip.text = element_blank(), 
        legend.margin=unit(-0.05,"cm"),
        legend.position = "bottom", 
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm"))

fig5f <- sc %>% filter(method %in% c("true","d0.01")) %>% 
  ggplot(aes(x = binnedcoeff, fill = method)) + 
  geom_density() + facet_grid(~del) + 
  scale_fill_manual(name = "", labels = c("True", "d0.01"), values = c("grey40", rgb(1, 1, 0, 0.4))) + 
  ylab("Density") + xlab("Scaled Selection Coefficients") + 
  theme(strip.text = element_blank(), 
        legend.margin=unit(-0.05,"cm"),
        legend.position = "bottom", 
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))

fig5g <- sc %>% filter(method %in% c("true","mvn1")) %>% 
  ggplot(aes(x = binnedcoeff, fill = method)) + 
  geom_density() + facet_grid(~del) + 
  scale_fill_manual(name = "", labels = c("True", "mvn1"), values = c("grey40", rgb(1, 1, 0, 0.4))) + 
  ylab("Density") + xlab("Scaled Selection Coefficients") + 
  theme(strip.text = element_blank(), 
        legend.margin=unit(-0.05,"cm"),
        legend.position = "bottom", 
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))


fig5h <- sc %>% filter(method %in% c("true","phylobayes")) %>% 
  ggplot(aes(x = binnedcoeff, fill = method)) + 
  geom_density() + facet_grid(~del) + 
  scale_fill_manual(name = "", labels = c("True", "pbMutSel"), values = c("grey40", rgb(1, 1, 0, 0.4))) + 
  ylab("Density") + xlab("Scaled Selection Coefficients") + 
  theme(strip.text = element_blank(), 
        legend.margin=unit(-0.05,"cm"),
        legend.position = "bottom", 
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))




fig5 <- plot_grid(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f, fig5g, fig5h, nrow=2, labels=c("A", "B", "C", "D", "E", "F", "G", "H"), vjust=1, scale=0.98)
save_plot(paste0(maintext_plot_directory, "dnds_sc_weakstrong.pdf"), fig5, base_width=13.25, base_height=4.25)





##########################################################################################
############# Figure 6: weak vs. strong dN/dS correlation and estimator bias #############
##########################################################################################
print("Figure 6")

all.corrs <- all.corrs.estbias %>% select(-b, -b.pvalue, -sig.bias) %>% spread(del, r)
all.corrs$method <- factor(all.corrs$method, levels = methods_levels, labels = methods_labels)
all.estbias <- all.corrs.estbias %>% select(-r, -b.pvalue, -sig.bias) %>% spread(del, b)
all.estbias$method <- factor(all.estbias$method, levels = methods_levels, labels = methods_labels)



theme_set(theme_cowplot() + theme(plot.margin = unit(c(1.0, 0.1, 0.1, 0.2),"cm"), 
                                  axis.text = element_text(size=11), 
                                  axis.title = element_text(size = 11, face="bold"), 
                                  legend.title = element_text(size = 10), 
                                  legend.text = element_text(size=9)))

p1 <- ggplot(all.corrs, aes(x = strong, y = weak, color = method)) + 
  geom_point() + geom_abline(slope=1, intercept=0) + 
  xlab("Strongly deleterious correlation") + 
  ylab("Weakly deleterious correlation  ") + 
  scale_x_continuous(limits=c(0.7, 1.0)) + 
  scale_y_continuous(limits=c(0.7, 1.0)) + 
  scale_color_brewer(palette = "Dark2", name= "Inference Method")

p2 <- ggplot(all.estbias, aes(x = strong, y = weak, color = method)) + 
  geom_point() + geom_abline(slope=1, intercept=0) + 
  xlab("Strongly deleterious bias") + 
  ylab("Weakly deleterious bias" )+ 
  scale_y_continuous(limits=c(-0.05, 0.32)) +
  scale_x_continuous(limits=c(-0.05, 0.32)) +
  scale_color_brewer(palette = "Dark2", name= "Inference Method") +
  geom_hline(yintercept=0, color="grey50") +
  geom_vline(xintercept=0, color="grey50") +
  theme(legend.position = "none")

grobs <- ggplotGrob(p1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
prow <- plot_grid(p1 + theme(legend.position = "none"), p2, labels=c("A", "B"), vjust=0.8, scale=0.98)
strong.weak.r.b <- plot_grid(prow, legend, rel_widths = c(2, .5))
save_plot(paste0(maintext_plot_directory, "scatter_strong_weak_rb.pdf"), strong.weak.r.b, base_width=8, base_height=2.75)










##########################################################################################
################################### SI plots #############################################
##########################################################################################
print("Creating SI plots")



##########################################################################################
########### Figures S1, S3: Predicted vs. true dN/dS for strong, weak inferences #########
##########################################################################################
print("Figures S1, S3")

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.1,3,0.1,0.1),"cm"), 
                                  axis.text = element_text(size = 11), 
                                  axis.title.x = element_text(size = 13), 
                                  axis.title.y = element_text(size = 11),
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.25, "cm"),
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=14, face="bold")))
scatters_strong <- list()
scatters_weak <- list()
i <- 1
for (d in datasets){
  splitd <- gsub("_", ", chain ", d)
  subdat <- spread.dnds %>% filter(dataset == d)
  subdat$method <- factor(subdat$method, levels = methods_levels, labels = methods_labels)
  
  
  pstrong.raw <- subdat %>% filter(del == "strong") %>% ggplot(aes(x = true, y = dnds)) + 
    geom_point(size=1.5) + 
    geom_abline(slope = 1, intercept = 0, color="red") + 
    xlab("True dN/dS") + ylab("Predicted dN/dS") + 
    scale_y_continuous(limits=c(0,0.85)) + 
    scale_x_continuous(limits=c(0,0.85)) + 
    facet_grid(~method)
  
  
  pweak.raw   <- subdat %>% filter(del == "weak") %>% ggplot(aes(x = true, y = dnds)) + 
    geom_point(size=1.5) + 
    geom_abline(slope = 1, intercept = 0, color="red") +
    xlab("True dN/dS") + ylab("Predicted dN/dS") + 
    scale_y_continuous(limits=c(0,0.85)) + 
    scale_x_continuous(limits=c(0,0.85)) + 
    facet_grid(~method)
  
  if (i != 1){
    pstrong.raw <- pstrong.raw + theme(strip.background = element_blank(), strip.text = element_blank())
    pweak.raw <- pweak.raw + theme(strip.background = element_blank(), strip.text = element_blank())
  }
  pstrong <- ggdraw(pstrong.raw) + draw_label(splitd, x = 0.95, y = 0.55, size=12)
  pweak <- ggdraw(pweak.raw) + draw_label(splitd, x = 0.95, y = 0.55, size=12)
  
  
  scatters_strong[[i]] <- pstrong
  scatters_weak[[i]] <- pweak
  i <- i+1
}
grid_strong <- plot_grid(plotlist = scatters_strong, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), vjust = c(2.5, rep(0.25,10)), rel_heights = c(0.105, rep(0.0895,10)))
save_plot(paste0(si_plot_directory, "strong_dnds_scatter_SI.pdf"), grid_strong, base_width=14, base_height=20)
grid_weak <- plot_grid(plotlist = scatters_weak, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), vjust = c(2.5, rep(0.25,10)), rel_heights = c(0.105, rep(0.0895,10)))
save_plot(paste0(si_plot_directory, "weak_dnds_scatter_SI.pdf"), grid_weak, base_width=14, base_height=20)





##########################################################################################
#### Figures S2, S4: Predicted and true S distributions for strong, weak inferences #####
##########################################################################################
print("Figures S2, S4")

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.1,3,0.1,0.1),"cm"), 
                                  axis.text = element_text(size=11), 
                                  axis.title = element_text(size=12), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.25, "cm"),
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=12, face="bold"),
                                  legend.position="none"))                                   

sc_strong <- list()
sc_weak <- list()
i <- 1
for (d in datasets){
  print(d)
  splitd <- gsub("_", ", chain ", d)
  strong.sc <- read_csv(paste0(result_directory, d, "_delstrong_selection_coefficients.csv"))
  true.strong <- strong.sc$binnedcoeff[strong.sc$method == "true"]
  strong.sc <- strong.sc %>% filter(method != "true")
  strong.sc$method <- factor(strong.sc$method, levels = methods_levels, labels = methods_labels)
  strong.sc$truebinned <- rep(true.strong, 8)
  
  weak.sc <- read_csv(paste0(result_directory, d, "_delweak_selection_coefficients.csv"))
  true.weak <- weak.sc$binnedcoeff[weak.sc$method == "true"]
  weak.sc <- weak.sc %>% filter(method != "true")
  weak.sc$method <- factor(weak.sc$method, levels = methods_levels, labels = methods_labels)
  weak.sc$truebinned <- rep(true.weak, 8)
  
  
  pstrong.raw <- ggplot(strong.sc) + geom_density(aes(x = truebinned), fill = "grey40") + 
    geom_density(aes(x = binnedcoeff), fill = "yellow", alpha = 0.4) + 
    xlab("Scaled Selection Coefficients") + ylab("Density") +
    facet_grid(~method) 
  
  pweak.raw <- ggplot(weak.sc) + geom_density(aes(x = truebinned), fill = "grey40") + 
    geom_density(aes(x = binnedcoeff), fill = "yellow", alpha = 0.4) + 
    xlab("Scaled Selection Coefficients") + ylab("Density") +
    facet_grid(~method) 
  
  if (i != 1){
    pstrong.raw <- pstrong.raw + theme(strip.background = element_blank(), strip.text = element_blank()) 
    pweak.raw <- pweak.raw + theme(strip.background = element_blank(), strip.text = element_blank()) 
  }
  pstrong <- ggdraw(pstrong.raw) + draw_label(splitd, x = 0.95, y = 0.5, size=12)
  pweak <- ggdraw(pweak.raw) + draw_label(splitd, x = 0.95, y = 0.55, size=12)
  
  
  sc_strong[[i]] <- pstrong
  sc_weak[[i]] <- pweak
  i <- i+1
}
grid_strong <- plot_grid(plotlist = sc_strong, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), vjust = c(2.5, rep(0.25,10)), rel_heights = c(0.105, rep(0.0895,10)))
save_plot(paste0(si_plot_directory, "strong_selcoeffs_SI.pdf"), grid_strong, base_width = 12, base_height=18)
grid_weak <- plot_grid(plotlist = sc_weak, nrow=11, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), vjust = c(2.5, rep(0.25,10)), rel_heights = c(0.105, rep(0.0895,10)))
save_plot(paste0(si_plot_directory, "weak_selcoeffs_SI.pdf"), grid_weak, base_width = 12, base_height=18)


# 
# # 
# 
# ##########################################################################################
# ######### Figure S5: Slope of JSD vs. true dN/dS jitter plots for weak simulations #######
# ##########################################################################################
# print("Figure S5")
# 
# true.jsd %>% filter(del == "weak") %>% 
#             group_by(dataset,method) %>% 
#             do(fit=lm(jsd~truednds, data=.)) %>% 
#             mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>% 
#             select(-fit) -> jsd.true.slope.p
# jsd.true.slope.p$method <- factor(jsd.true.slope.p$method, levels = methods_levels, labels = methods_labels)
# 
# p <- ggplot(jsd.true.slope.p, aes(x = method, y = slope, shape=sig)) + 
#             geom_point(position = position_jitter(w = 0.6)) + 
#             scale_shape_manual(values=c(1,19)) + 
#             geom_hline(yintercept=0) + 
#             scale_y_continuous(limits=c(-0.25, 0.1)) +
#             xlab("Inference Method") + ylab("Slope") + 
#             theme(legend.position = "none", 
#                   axis.text.x = element_text(size = 11), 
#                   axis.text.y = element_text(size = 10), 
#                   axis.title = element_text(size = 13))
# 
# save_plot(paste0(si_plot_directory, "weak_jsd_dnds_SI.pdf"), p, base_width=9, base_height=2)
# 
# 
