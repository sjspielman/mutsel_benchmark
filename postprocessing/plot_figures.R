# SJS
# Creates MS figures (output to either figures/main_text or figures/SI, depending on where figure is located in MS)

require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(grid)

result_directory <- "../results/summarized_results/"
true_directory   <- "../simulation/true_simulation_parameters/"
maintext_plot_directory   <- "figures/maintext/"
si_plot_directory   <- "figures/SI/"

jsd <- read.csv(paste0(result_directory, "jsd.csv"))
dnds <- read.csv(paste0(result_directory, "dnds.csv"))

methods_levels <- c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")
methods_labels <- c("Unpenalized", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "pbMutSel")
alpha <- 0.01 # Significance
corrected.alpha <- alpha/length(methods_levels) #Bonferroni significance
repr_sim <- "1R6M_A"
datasets <- unique(dnds$dataset)


dnds %>% spread(method,dnds) %>% 
         gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
         select(dataset, del, site, true, dnds, method) -> spread.dnds
true.jsd <- dnds %>% filter(method=="true") %>% 
                     mutate(truednds = dnds) %>% 
                     select(-method,-dnds) %>% left_join(jsd)

spread.dnds %>% group_by(dataset, del, method) %>%
                do(rraw = cor(.$true, .$dnds), braw = glm(dnds ~ offset(true), dat=.)) %>% 
                mutate(r = rraw[[1]], b = summary(braw)$coeff[1], b.pvalue = summary(braw)$coeff[4], sig.bias = b.pvalue<corrected.alpha) %>% 
                select(-rraw, -braw) -> all.corrs.estbias



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

##########################################################################################
############## Figure 1: Selection coefficients for representative dataset ###############
##########################################################################################
print("Figure 1")

strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
truebinned <- strong.sc.full$binnedcoeff[strong.sc.full$method == "true"]
strong.sc.full %>% filter(method != "true") -> strong.sc.full
strong.sc.full$truebinned <- truebinned 
strong.sc.full$method <- factor(strong.sc.full$method, levels = methods_levels, labels = methods_labels)

theme_set(theme_cowplot() + theme(axis.text.y = element_text(size=12), 
                                  axis.text.x = element_text(size=10), 
                                  axis.title = element_text(size=13, face="bold"), 
                                  panel.border = element_rect(size = 0.5), 
                                  strip.text = element_text(size = 11), 
                                  panel.margin = unit(0.75, "lines"),
                                  strip.background = element_rect(fill="white")))

plot_for_legend <- strong.sc.full %>% filter(method %in% c("mvn1", "mvn100")) %>%
                      ggplot(aes(x = binnedcoeff, fill=method)) + geom_density() + 
                      scale_fill_manual(name = "", labels = c("True", "Inferred"), values =c("grey40", rgb(1, 1, 0, 0.4))) + 
                      theme(legend.position="bottom", 
                            legend.margin = unit(0, "cm"), 
                            legend.text = element_text(size=9), 
                            legend.key.size = unit(0.4, "cm")) 

grobs <- ggplotGrob(plot_for_legend)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

p <-ggplot(strong.sc.full, aes(x = binnedcoeff)) + 
           geom_density(aes(x = truebinned), fill = "grey40") + 
           geom_density(fill = "yellow", alpha = 0.4) + 
           facet_grid(~method) + xlab("Scaled Selection Coefficients") + ylab("Density") + 
           theme(plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
           
figsc_with_legend <- plot_grid(p, legend, nrow=2, rel_heights=c(1,0.1))
save_plot(paste0(maintext_plot_directory, "sc_across_methods.pdf"), figsc_with_legend, base_width=12, base_height=2.5)





##########################################################################################
########## Figure 2: Boxplots of JSD for representative dataset, all datasets ############
##########################################################################################
print("Figure 2")

jsd %>% filter(del == "strong") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.strong.summary
jsd.strong.summary$method <- factor(jsd.strong.summary$method, levels = methods_levels, labels = methods_labels)
jsd %>% filter(del == "weak") %>% group_by(dataset, method) %>% summarize(meanjsd = mean(jsd)) -> jsd.weak.summary
jsd.weak.summary$method <- factor(jsd.weak.summary$method, levels = methods_levels, labels = methods_labels)
jsd %>% filter(del == "strong", dataset == repr_sim) -> jsd.repr.sim
jsd.repr.sim$method <- factor(jsd.repr.sim$method, levels = methods_levels, labels = methods_labels)

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 10), 
                                  axis.text.y = element_text(size = 12), axis.title = element_text(size = 13)))
## Figure 2A
repr.jsd.boxplots <- ggplot(jsd.repr.sim, aes(x = method, y = jsd)) + 
                            geom_boxplot() + 
                            xlab("Inference Method") + ylab("Site JSD") + 
                            scale_y_continuous(limits = c(0, 0.45), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) 
mean.jsd.boxplots <- ggplot(jsd.strong.summary, aes(x = method, y = meanjsd)) + 
                            geom_boxplot() + 
                            xlab("Inference Method") + ylab("Average JSD") + 
                            scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))

fig2 <- plot_grid(repr.jsd.boxplots, mean.jsd.boxplots, nrow=2, labels=c("A", "B"), label_size = 17, scale=0.925)
save_plot(paste0(maintext_plot_directory, "jsd_maintext.pdf"), fig2, base_width = 7, base_height = 6) 




##########################################################################################
####### Figure 3: Scatterplots of true, predicted dN/dS for representative dataset #######
##########################################################################################
print("Figure 3")

dnds.repr.strong <- spread.dnds %>% filter(dataset == repr_sim, del == "strong")
dnds.repr.strong$method <- factor(dnds.repr.strong$method, levels = methods_levels, labels = methods_labels)
dnds.repr.weak <- spread.dnds %>% filter(dataset == repr_sim, del == "weak")
dnds.repr.weak$method <- factor(dnds.repr.weak$method, levels = methods_levels, labels = methods_labels)

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 8.5), 
                                  axis.title = element_text(size = 11, face = "bold"), 
                                  panel.border = element_rect(size = 0.5), 
                                  panel.margin = unit(0.25, "cm"), 
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=10)))

fig3 <- ggplot(dnds.repr.strong, aes(x = true, y = dnds)) + 
        geom_point(size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
        xlab("True dN/dS") + ylab("Predicted dN/dS") + 
        scale_y_continuous(limits=c(0,0.88)) + scale_x_continuous(limits=c(0,0.88)) + 
        facet_grid(~method) 
save_plot(paste0(maintext_plot_directory, "repr_dnds_scatter.pdf"), fig3, base_width=10, base_height=2)




##########################################################################################
####### Figure 4: Jitter plots of correlations, estimator bias for predicted dN/dS  ######
##########################################################################################
print("Figure 4")

theme_set(theme_cowplot() + theme(axis.text.x = element_text(size = 11), 
                                  axis.text.y = element_text(size = 12), 
                                  axis.title = element_text(size = 14, face="bold")))


strong.corrs.estbias <- all.corrs.estbias %>% filter(del == "strong") 
strong.corrs.estbias$method <- factor(strong.corrs.estbias$method, levels = methods_levels, labels = methods_labels)

# all are significant, so no need for shape attribute
jitter.r <- ggplot(strong.corrs.estbias, aes(x = method, y = r)) + 
                    geom_jitter(w = 0.2) + 
                    xlab("Inference Method") + ylab("Pearson Correlation") +
                    scale_y_continuous(limits=c(0.7, 1.0))

jitter.b <- ggplot(strong.corrs.estbias, aes(x = method, y = b, shape = sig.bias)) + 
                    geom_jitter(w = 0.6) + 
                    scale_shape_manual(values=c(1,19)) + 
                    xlab("Inference Method") + ylab("Estimator Bias") + 
                    geom_hline(yintercept=0 ) + theme(legend.position="none")
jitter.r.b <- plot_grid(jitter.r, jitter.b, nrow=2, labels=c("A", "B"), label_size=17, scale=0.925)
save_plot(paste0(maintext_plot_directory, "strong_r_bias.pdf"), jitter.r.b, base_width=7.5, base_height=6)






##########################################################################################
##### Figure 5: JSD regressed on true dN/dS for representative, and slopes for all #######
##########################################################################################
print("Figure 5")

strong <- true.jsd %>% filter(del == "strong", dataset == repr_sim)
strong$method <- factor(strong$method, levels = methods_levels, labels = methods_labels)
theme_set(theme_cowplot() + theme(axis.text = element_text(size = 9),  
                                  axis.title = element_text(size = 13, face = "bold"), 
                                  panel.border = element_rect(size = 0.5),
                                  strip.background = element_rect(fill="white"), 
                                  strip.text = element_text(size=10)))

fig5a <- ggplot(strong, aes(x = truednds, y = jsd)) + 
         geom_point(size=1) + 
         facet_grid(~method) + 
         geom_smooth(method="lm", color="red") + 
         scale_y_continuous(limits=c(0, 0.45)) + 
         xlab("True dN/dS") + ylab("Site JSD")
         
true.jsd %>% filter(del == "strong") %>% 
             group_by(dataset,method) %>% 
             do(fit=lm(jsd~truednds, data=.)) %>% 
             mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>% 
             select(-fit) -> jsd.true.slope.p
jsd.true.slope.p$method <- factor(jsd.true.slope.p$method, levels = methods_levels, labels = methods_labels)
fig5b <- ggplot(jsd.true.slope.p, aes(x = method, y = slope, shape=sig)) + 
                geom_point(position = position_jitter(w = 0.4)) + 
                scale_shape_manual(values=c(1,19)) + 
                geom_hline(yintercept=0) + 
                theme(legend.position="none") + 
                xlab("Inference Method") + ylab("Slope")+ 
                theme(axis.text.x = element_text(size = 11), 
                      axis.text.y = element_text(size = 10))

fig5 <- plot_grid(fig5a, fig5b, nrow=2, labels=c("A", "B"))
save_plot(paste0(maintext_plot_directory, "jsd_dnds.pdf"), fig5, base_width=9.5, base_height=4)






##########################################################################################
##### Figure 6: weak vs. strong scatterplots and selection coefficient distributions #####
##########################################################################################
print("Figure 6")

strong.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delstrong_selection_coefficients.csv"))
weak.sc.full <- read_csv(paste0(result_directory, repr_sim, "_delweak_selection_coefficients.csv"))
strong.sc <- filter(strong.sc.full, method %in% c("true", "nopenal", "phylobayes"))
weak.sc <- filter(weak.sc.full, method %in% c("true", "nopenal", "phylobayes"))
sc <- rbind(strong.sc, weak.sc)
sc$method <- factor(sc$method, levels=c("true", "nopenal", "phylobayes"))
spread.dnds.named <- spread.dnds
spread.dnds.named$del <- factor(spread.dnds.named$del, levels=c("strong", "weak"), labels = c("Strongly deleterious", "Weakly deleterious"))

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.3, 0.1, 0.1, 0.2),"cm"), 
                                  panel.margin = unit(0.25, "cm"),
                                  axis.text = element_text(size=9), 
                                  axis.title = element_text(size=9.5, face = "bold"), 
                                  strip.background = element_rect(fill="white"),
                                  strip.text = element_text(size = 9)))

fig6a <- spread.dnds.named %>% filter(dataset == repr_sim, method == "nopenal") %>% 
         ggplot(aes(x = true, y = dnds)) + 
         geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
         facet_grid(~del) + xlab("True dN/dS") + ylab("swMutSel dN/dS") +   
         scale_x_continuous(limits = c(0,0.9)) + 
         scale_y_continuous(limits = c(0,0.9))
         
fig6b <- spread.dnds.named %>% filter(dataset == repr_sim, method == "phylobayes") %>% 
         ggplot(aes(x = true, y = dnds)) + 
         geom_point() + geom_abline(slope = 1, intercept = 0, color="red") + 
         facet_grid(~del) + xlab("True dN/dS") + ylab("pbMutSel dN/dS") + 
         scale_x_continuous(limits = c(0,0.9)) + 
         scale_y_continuous(limits = c(0,0.9))
                                  
                                  
fig6c <- sc %>% filter(method %in% c("true","nopenal")) %>% 
                ggplot(aes(x = binnedcoeff, fill = method)) + 
                geom_density() + facet_grid(~del) + 
                scale_fill_manual(name = "", labels = c("True", "Inferred"), values =c("grey40", rgb(1, 1, 0, 0.4))) + 
                ylab("Density") + xlab("Scaled Selection Coefficients") + 
                theme(strip.text = element_blank(), 
                      legend.margin=unit(-0.05,"cm"),
                      legend.position = "bottom", 
                      legend.text = element_text(size=8), 
                      legend.key.size = unit(0.3, "cm"))
                                         
fig6d <- sc %>% filter(method %in% c("true","phylobayes")) %>% 
                ggplot(aes(x = binnedcoeff, fill = method)) + 
                geom_density() + facet_grid(~del) + 
                scale_fill_manual(name = "", labels = c("True", "Inferred"), values = c("grey40", rgb(1, 1, 0, 0.4))) + 
                ylab("Density") + xlab("Scaled Selection Coefficients") + 
                theme(strip.text = element_blank(), 
                      legend.margin=unit(-0.05,"cm"),
                      legend.position = "bottom", 
                      legend.text = element_text(size=8),
                      legend.key.size = unit(0.3, "cm"))

fig6 <- plot_grid(fig6a, fig6b, fig6c, fig6d, nrow=2, labels=c("A", "B", "C", "D"), label_size = 13, scale=0.97)
save_plot(paste0(maintext_plot_directory, "dnds_sc_weakstrong.pdf"), fig6, base_width=7, base_height=4)





##########################################################################################
############# Figure 7: weak vs. strong dN/dS correlation and estimator bias #############
##########################################################################################
print("Figure 7")

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

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.1,2.2,0.1,0.1),"cm"), 
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
    pstrong <- ggdraw(pstrong.raw) + draw_label(d, x = 0.97, y = 0.5, size=12)
    pweak <- ggdraw(pweak.raw) + draw_label(d, x = 0.97, y = 0.5, size=12)


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

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.1,2.2,0.1,0.1),"cm"), 
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
    pstrong <- ggdraw(pstrong.raw) + draw_label(d, x = 0.96, y = 0.5, size=12)
    pweak <- ggdraw(pweak.raw) + draw_label(d, x = 0.96, y = 0.5, size=12)

    
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
