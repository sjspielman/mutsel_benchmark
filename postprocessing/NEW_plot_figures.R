# SJS
# Creates MS figures (output to either figures/main_text or figures/SI, depending on where figure is located in MS)

library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(grid)



################### THIS FUNCTION IS FROM THE PACKAGE "SMATR", which cannot be loaded because it breaks my linear models. ################
slope.test <- function( y, x, test.value=1, V=matrix(0,2,2))
{
   
    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)

    resDF <- n - 2
    fCrit <- qf( 1-alpha, 1, resDF )

    dat <- cbind(y[iref], x[iref])
    r.factor <- 1
	  vr <- ( cov(dat) - V )*(n-1)
		vr <- t(dat)%*%dat - V*n
    r <- vr[1,2]/sqrt( vr[1,1]*vr[2,2] )

    bCI     <- matrix( NA, 1, 2 )
    varTest <- matrix( 0, 2, 2 )

     # linear regression code only. SJS
     b            <- vr[1,2]/vr[2,2]
     varRes       <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] )/resDF
     varB         <- varRes/vr[2,2] * r.factor
     bCI[1,1]     <- b - sqrt(varB)*sqrt(fCrit)
     bCI[1,2]     <- b + sqrt(varB)*sqrt(fCrit)
     varTest[1,1] <- vr[1,1] - 2*test.value*vr[1,2] + test.value^2*vr[2,2]
     varTest[1,2] <- vr[1,2] - test.value*vr[2,2]
     varTest[2,2] <- vr[2,2]
     rTest  <- varTest[1,2] / sqrt( varTest[1,1] ) / sqrt( varTest[2,2] )
     F      <- rTest^2/(1 - rTest^2)/r.factor*(n-2)
     pValue <- 1 - pf( F, 1, resDF)
     list( p=pValue, b=b )

}



# function to return pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Linear models on dnds, entropy inferences
summarize_dnds_entropy <- function(dat, siglevel){
  dat %>% group_by(dataset, type, method, bl) %>%
    do(rraw.dnds = cor(.$dnds, .$true.dnds),
       braw.dnds = glm(dnds ~ offset(true.dnds), dat=.),
       rraw.h = cor(.$entropy, .$true.entropy), 
       braw.h = glm(entropy ~ offset(true.entropy), dat=.))   %>%
    mutate(r2.dnds   = rraw.dnds[[1]]^2,
           b.dnds = summary(braw.dnds)$coeff[[1]],
           b.dnds.sig = summary(braw.dnds)$coeff[4]<siglevel,
           r2.entropy = rraw.h[[1]]^2,
           b.entropy = summary(braw.h)$coeff[1],
           b.entropy.sig = summary(braw.h)$coeff[4]<siglevel) %>%
    select(-rraw.dnds, -rraw.h, -braw.dnds, -braw.h) -> dat.sum1

  
  dat %>% group_by(dataset, method) %>% 
    do(x = slope.test(.$dnds, .$true.dnds, test.value=1)) %>%
    mutate(slope.dnds.sig = x$p<siglevel, slope.dnds = x$b) %>%
    select(-x) -> dat.sum2
  
  dat %>% group_by(dataset, method) %>% 
    do(x = slope.test(.$entropy, .$true.entropy, test.value=1)) %>%
    mutate(slope.entropy.sig = x$p<siglevel, slope.entropy = x$b) %>%
    select(-x) -> dat.sum3
  
  part <- left_join(dat.sum2, dat.sum3)
  dat.sum <- left_join(dat.sum1, part)
  
  dat.sum
}

theme_set(theme_cowplot() + theme(panel.margin = unit(0.75, "lines")))

result_directory         <- "dataframes/"
true_directory           <- "../simulation/true_simulation_parameters/"
maintext_plot_directory  <- "NEWfigures/"
si_plot_directory        <- "NEWfigures/"
repr.datasets            <- c("1R6M_A", "LAC")
methods_levels           <- c("nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes") 
methods_labels           <- c("Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel")
alpha                    <- 0.01 # Significance
corrected.alpha          <- alpha/length(methods_levels) #Bonferroni significance


dat <- read.csv(paste0(result_directory, "inference_results.csv"))
datasets <- c(as.character(unique(dat$dataset)))


#dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
#dat$bl <- factor(dat$bl, levels = c(0.01, 0.5))
type.colors <- c("red", "blue")
dat <- dat %>% mutate(type = ifelse(dataset %in% c("NP", "HA", "LAC"), "dms", "yeast"))
dat$type <- factor(dat$type, levels=c("yeast", "dms"), labels = c("Yeast", "DMS"))
dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
sub.dat <- dat %>% filter(true.dnds >= 0.4, true.dnds <= 0.75)



sumstats <- dat %>% summarize_dnds_entropy(corrected.alpha)
sub.sumstats <- sub.dat %>% summarize_dnds_entropy(corrected.alpha)
###################### JSD and diff sum #######################

dat %>% group_by(dataset, type, method, bl) %>% summarize(meanjsd = mean(jsd), meandiffsum = mean(diffsum))-> jsd.diffsum
jsd.diffsum$method <- factor(jsd.diffsum$method, levels = methods_levels, labels = methods_labels)

jsd.diffsum %>% 
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type))) + 
  geom_jitter(width = 0.5) +
  scale_color_manual(values=type.colors, name = "Data type ") + 
  theme(legend.position = "bottom") -> type.legend.raw
grobs <- ggplotGrob(type.legend.raw)$grobs
type.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type))) + 
  geom_jitter(width = 0.5) + ggtitle("BL = 0.5") +
  scale_color_manual(values=type.colors, name = "Data type") + 
  xlab("Inference Method") + ylab("Average JSD") + theme(legend.position = "none", axis.text.x = element_text(size=11)) -> jsd.jitters.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type))) + 
  geom_jitter(width = 0.5) + ggtitle("BL = 0.01") +
  scale_color_manual(values=type.colors, name = "Data type") + 
  xlab("Inference Method") + ylab("Average JSD") + theme(legend.position = "none", axis.text.x = element_text(size=11)) -> jsd.jitters.raw2
jsd.jitters1 <- plot_grid(jsd.jitters.raw1, jsd.jitters.raw2, nrow=1, labels=c("A", "B"), scale=0.98)
jsd.jitters  <- plot_grid(jsd.jitters1, type.legend, ncol=1, rel_heights=c(1, 0.07))
save_plot(paste0(maintext_plot_directory, "jsd_jitters.pdf"), jsd.jitters, base_width = 11, base_height=3)


jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = meandiffsum, color = as.factor(type))) + 
  geom_jitter(width = 0.5) + ggtitle("BL = 0.5") +
  scale_color_manual(values=type.colors, name = "Data type") + 
  xlab("Inference Method") + ylab("Average Difference") + theme(legend.position = "none") -> diffsum.jitters.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = method, y = meandiffsum, color = as.factor(type))) + 
  geom_jitter(width = 0.5) + ggtitle("BL = 0.01") +
  scale_color_manual(values=type.colors, name = "Data type") + 
  xlab("Inference Method") + ylab("Average Difference") + theme(legend.position = "none") -> diffsum.jitters.raw2
diffsum.jitters1 <- plot_grid(diffsum.jitters.raw1, diffsum.jitters.raw2, nrow=1, labels=c("A", "B"), scale=0.98)
diffsum.jitters  <- plot_grid(diffsum.jitters1, type.legend, ncol=1, rel_heights=c(1, 0.07))
save_plot(paste0(maintext_plot_directory, "diffsum_jitters.pdf"), diffsum.jitters, base_width = 11.5, base_height=3.5)


###################### dN/dS and entropy scatterplots #######################

sumstats %>% 
    select(-b.dnds, -b.dnds.sig, -b.entropy, -b.entropy.sig, -slope.dnds.sig, -slope.dnds, -slope.entropy.sig, -slope.entropy) %>%
    filter(dataset %in% repr.datasets) %>% select(-dataset) %>% group_by(type, method) %>%
    mutate(r2.dnds.short = as.numeric(format(r2.dnds, digits=3)), r2.entropy.short = as.numeric(format(r2.entropy, digits=3))) -> scatter.stats
scatter.stats %>%  -> scatter.stats
stats.bl0.5 <- scatter.stats %>% filter(bl == 0.5)
stats.bl0.01 <- scatter.stats %>% filter(bl == 0.01)
dat %>% filter(dataset %in% repr.datasets, bl == 0.5) %>% select(-dataset) -> scatter.repr.bl0.5
dat %>% filter(dataset %in% repr.datasets, bl == 0.01)%>% select(-dataset)  -> scatter.repr.bl0.01

theme_set(theme_cowplot() + theme(strip.text = element_text(size = 13),
                                panel.margin = unit(0.75, "lines"),
                                axis.text = element_text(size=10)))
              
scatter.dnds.bl0.5 <- ggplot(data = NULL) + 
  geom_point(data = scatter.repr.bl0.5, aes(x = true.dnds, y = dnds), size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
  geom_text(data = stats.bl0.5, aes(label=paste0("r^2==",r2.dnds.short)), x = 0.25, y = 0.88, parse=TRUE, size=5) +
  xlab("True dN/dS") + ylab("Predicted dN/dS") + 
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,1)) + 
  facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_dnds_bl0.5.pdf"), scatter.dnds.bl0.5, base_width=15.5, base_height=5.5)

scatter.dnds.bl0.01 <- ggplot(data = NULL) + 
    geom_point(data = scatter.repr.bl0.01, aes(x = true.dnds, y = dnds), size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
    geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.dnds.short)), x = 0.25, y = 1.06, parse=TRUE, size=4.75) +
    xlab("True dN/dS") + ylab("Predicted dN/dS") + 
    scale_y_continuous(limits=c(0,1.1)) + scale_x_continuous(limits=c(0,1.1)) + 
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_dnds_bl0.01.pdf"), scatter.dnds.bl0.01, base_width=15.5, base_height=5.5)


scatter.entropy.bl0.5 <- ggplot(data = NULL) + 
    geom_point(data = scatter.repr.bl0.5, aes(x = true.entropy, y = entropy), size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
    geom_text(data = stats.bl0.5, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.6, y = 2.5, parse=TRUE, size=4) +
    xlab("True Entropy") + ylab("Predicted Entropy") + 
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_entropy_bl0.5.pdf"), scatter.entropy.bl0.5, base_width=15.5, base_height=5.5)

scatter.entropy.bl0.01 <- ggplot(data = NULL) + 
    geom_point(data = scatter.repr.bl0.01, aes(x = true.entropy, y = entropy), size=1) + geom_abline(slope = 1, intercept = 0, color="red") + 
    geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.6, y = 2.9, parse=TRUE, size=4) +
    xlab("True Entropy") + ylab("Predicted Entropy") + 
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_entropy_bl0.01.pdf"), scatter.entropy.bl0.01, base_width=15.5, base_height=5.5)




theme_set(theme_cowplot() + theme(strip.text = element_text(size = 13),
                                  panel.margin = unit(0.75, "lines"),
                                  axis.text = element_text(size=10)))

jitter.r2.dnds <- sumstats %>% # All r2 are significant
  ggplot(aes(x = method, y = r2.dnds, color = as.factor(type))) + 
  geom_jitter(width = 0.6, size=2, alpha=0.7) + 
  scale_color_manual(values=c("red", "blue")) +
  facet_grid(~bl) + theme(legend.position = "none") + 
  xlab("Inference Method") + ylab(expression(r^2)) 
jitter.b.dnds <- sumstats %>%
    ggplot(aes(x = method, y = b.dnds, color = as.factor(type), shape = as.factor(b.dnds.sig))) + 
    geom_jitter(width = 0.6, size=2, alpha=0.7) + 
    scale_color_manual(values=c("red", "blue")) + 
    scale_shape_manual(values=c(1,19)) +
    geom_hline(yintercept=0) + 
    facet_grid(~bl) + theme(legend.position = "none") + 
    xlab("Inference Method") + ylab("Estimator Bias")
jitter.slope.dnds <- sumstats %>%   # All slopes are significant
  ggplot(aes(x = method, y = slope.dnds, color = as.factor(type))) + 
  geom_jitter(width = 0.6, size=2, alpha=0.7) + 
  scale_color_manual(values=c("red", "blue")) +
  geom_hline(yintercept=0) + 
  facet_grid(~bl) + theme(legend.position = "none") + 
  xlab("Inference Method") + ylab("True-Inferred Slope")
r2.b.slope.dnds <- plot_grid(jitter.r2.dnds, jitter.b.dnds, jitter.slope.dnds, nrow=3, labels=c("A", "B", "C"))
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_dnds.pdf"), r2.b.slope.dnds, base_width=8.5, base_height=7.5)

jitter.r2.entropy <- sumstats %>%
  ggplot(aes(x = method, y = r2.entropy, color = as.factor(type))) + 
  geom_jitter(width = 0.6, size=2, alpha=0.7) + 
  scale_color_manual(values=c("red", "blue")) +
  facet_grid(~bl) + theme(legend.position = "none") + 
  xlab("Inference Method") + ylab(expression(r^2)) 
jitter.b.entropy <- sumstats %>%
  ggplot(aes(x = method, y = b.entropy, color = as.factor(type), shape = as.factor(b.entropy.sig))) + 
  geom_jitter(width = 0.6, size=2, alpha=0.7) + 
  scale_color_manual(values=c("red", "blue")) + 
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + 
  facet_grid(~bl) + theme(legend.position = "none") + 
  xlab("Inference Method") + ylab("Estimator Bias")
jitter.slope.entropy <- sumstats %>% 
  ggplot(aes(x = method, y = slope.entropy, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
  geom_jitter(width = 0.6, size=2, alpha=0.7) + 
  scale_color_manual(values=c("red", "blue")) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + 
  facet_grid(~bl) + theme(legend.position = "none") + 
  xlab("Inference Method") + ylab("True-Inferred Slope")
r2.b.slope.entropy <- plot_grid(jitter.r2.entropy, jitter.b.entropy, jitter.slope.entropy, nrow=3, labels=c("A", "B", "C"))
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_entropy.pdf"), r2.b.slope.entropy, base_width=8.5, base_height=7.5)










##########################################################################################
##### Figure 4: JSD regressed on true dN/dS for representative, and slopes for all #######
##########################################################################################
print("Figure 4")
dms.bl0.5 <- dms %>% filter(bl == 0.5)
yeast %>% filter(dataset == repr_sim, bl == 0.5) %>% 
  select(-dataset) %>% 
  mutate(dataset = del) %>% 
  select(-del) %>% 
  rbind(dms.bl0.5) %>%
  ggplot(aes(x = true.dnds, y = jsd)) + 
  geom_point(size=1) + facet_grid(dataset~method) + 
  geom_smooth(method="lm", color="red") + 
  xlab("True dN/dS") + ylab("Site JSD") -> jsd.dnds.reprs
save_plot(paste0(maintext_plot_directory, "jsd_dnds_scatter.pdf"), jsd.dnds.reprs,base_width = 11, base_height=8)



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

fig4 <- plot_grid(, dnds.jsd.jitter.bl0.5, nrow=2, labels=c("A", "B"))
save_plot(paste0(maintext_plot_directory, "jsd_dnds_scatter.pdf"), dnds.on.jsd.bl0.5, base_width=9.5, base_height=4)




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
