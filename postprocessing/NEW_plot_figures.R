o# SJS
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
     list( pValue, b )

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
    mutate(slope.dnds.sig = x[[1]]<siglevel, slope.dnds = x[[2]]) %>%
    select(-x) -> dat.sum2
  dat %>% group_by(dataset, method) %>%
    do(x = slope.test(.$entropy, .$true.entropy, test.value=1)) %>%
    mutate(slope.entropy.sig = x[[1]]<siglevel, slope.entropy = x[[2]]) %>%
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
dms.datasets             <- c("GAL4", "LAC", "HA", "NP")
methods_levels           <- c("nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes")
methods_labels           <- c("Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel")
alpha                    <- 0.01 # Significance
corrected.alpha          <- alpha/length(methods_levels) #Bonferroni significance


dat <- read.csv(paste0(result_directory, "inference_results.csv"))
true.selcoeffs <- read_csv(paste0(result_directory, "true_selection_coefficients.csv"))
datasets <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A", "GAL4", "LAC", "HA", "NP")
sub.datasets <- c("1R6M_A", "GAL4", "LAC", "HA", "NP")


dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
dat$bl <- factor(dat$bl, levels = c(0.01, 0.5))
#dat$datasets <- factor(dat$dataset, levels=datasets)
type.colors <- c("red", "blue")
dat <- dat %>% mutate(type = ifelse(!(dataset %in% dms.datasets), "Yeast", "DMS"))
dat$type <- factor(dat$type, levels=c("Yeast", "DMS"))

sumstats <- dat %>% summarize_dnds_entropy(corrected.alpha)


##################### Distributions of selection coefficients and selection pressures ( meaning dN/dS, entropy) ######################
true.selcoeffs %>% filter(dataset %in% sub.datasets) -> sub.selcoeffs
sub.selcoeffs$dataset <- factor(sub.selcoeffs$dataset, levels=sub.datasets)
dat %>% filter(dataset %in% sub.datasets, method == "mvn10") -> sub.dat
sub.dat$dataset <- factor(sub.dat$dataset, levels=sub.datasets)

compare.sc.densities <- sub.selcoeffs %>% ggplot(aes(x = binnedcoeffs)) + geom_density() + facet_grid(~dataset, scales="free_y") + scale_x_continuous(limits=c(-10,10)) + xlab("True S")
compare.truednds     <- sub.dat %>% ggplot(aes(x = dataset, y = true.dnds)) + geom_boxplot() + xlab("Dataset") + ylab("True dN/dS")
compare.trueh        <- sub.dat %>% ggplot(aes(x = dataset, y = true.entropy)) + geom_boxplot() + xlab("Dataset") + ylab("True Entropy")
compare.true.dnds_entropy <- plot_grid(compare.truednds,compare.trueh, nrow=1, labels="AUTO")
save_plot(paste0(maintext_plot_directory, "compare_sc.pdf"), compare.sc.densities, base_width=9, base_height=3)
save_plot(paste0(si_plot_directory, "compare_true_dnds_entropy.pdf"), compare.true.dnds_entropy, base_width=8, base_height=3)


###################### JSD and diff sum #######################

dat %>% group_by(dataset, type, method, bl) %>% summarize(meanjsd = mean(jsd), meandiffsum = mean(diffsum))-> jsd.diffsum

dat %>% filter(true.dnds >= 0.25) %>%
sub.dat %>% group_by(dataset, type, method, bl) %>% summarize(meanjsd = mean(jsd), meandiffsum = mean(diffsum))-> sub.jsd.diffsum

jsd.diffsum %>%
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type), group = dataset)) +
  geom_point(size=2) + geom_line() +
  scale_color_manual(values=type.colors, name = "Data type ") +
  theme(legend.position = "bottom") -> type.legend.raw
grobs <- ggplotGrob(type.legend.raw)$grobs
type.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type), group = dataset)) +
  geom_point(size=1.5) + geom_line() + ggtitle("BL = 0.5") +
  scale_color_manual(values=type.colors, name = "Data type") +
  xlab("Inference Method") + ylab("Average JSD") + theme(legend.position = "none", axis.text.x = element_text(size=11)) -> jsd.jitters.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = method, y = meanjsd, color = as.factor(type), group = dataset)) +
  geom_point(size=1.5) + geom_line() + ggtitle("BL = 0.01") +
  scale_color_manual(values=type.colors, name = "Data type") +
  xlab("Inference Method") + ylab("Average JSD") + theme(legend.position = "none", axis.text.x = element_text(size=11)) -> jsd.jitters.raw2
jsd.jitters1 <- plot_grid(jsd.jitters.raw1, jsd.jitters.raw2, nrow=1, labels=c("A", "B"), scale=0.98)
jsd.jitters  <- plot_grid(jsd.jitters1, type.legend, ncol=1, rel_heights=c(1, 0.09))
save_plot(paste0(maintext_plot_directory, "jsd_lineplot.pdf"), jsd.jitters, base_width = 11, base_height=3)


jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = meandiffsum, color = as.factor(type), group = dataset)) +
  geom_point(size=1.5) + geom_line() + ggtitle("BL = 0.5") +
  scale_color_manual(values=type.colors, name = "Data type") +
  xlab("Inference Method") + ylab("Average Difference") + theme(legend.position = "none") -> diffsum.jitters.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = method, y = meandiffsum, color = as.factor(type), group = dataset)) +
  geom_point(size=1.5) + geom_line() + ggtitle("BL = 0.01") +
  scale_color_manual(values=type.colors, name = "Data type") +
  xlab("Inference Method") + ylab("Average Difference") + theme(legend.position = "none") -> diffsum.jitters.raw2
diffsum.jitters1 <- plot_grid(diffsum.jitters.raw1, diffsum.jitters.raw2, nrow=1, labels=c("A", "B"), scale=0.98)
diffsum.jitters  <- plot_grid(diffsum.jitters1, type.legend, ncol=1, rel_heights=c(1, 0.09))
save_plot(paste0(maintext_plot_directory, "diffsum_lineplot.pdf"), diffsum.jitters, base_width = 11.5, base_height=3.5)


###################### dN/dS and entropy scatterplots #######################

sumstats %>%
    select(-b.dnds, -b.dnds.sig, -b.entropy, -b.entropy.sig, -slope.dnds.sig, -slope.dnds, -slope.entropy.sig, -slope.entropy) %>%
    filter(dataset %in% repr.datasets) %>% select(-dataset) %>% group_by(type, method) %>%
    mutate(r2.dnds.short = as.numeric(format(r2.dnds, digits=3)), r2.entropy.short = as.numeric(format(r2.entropy, digits=3))) -> scatter.stats
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
    geom_point(data = scatter.repr.bl0.01, aes(x = true.dnds, y = dnds), size=1) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.dnds.short)), x = 0.25, y = 1.06, parse=TRUE, size=4.75) +
    xlab("True dN/dS") + ylab("Predicted dN/dS") +
    scale_y_continuous(limits=c(0,1.1)) + scale_x_continuous(limits=c(0,1.1)) +
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_dnds_bl0.01.pdf"), scatter.dnds.bl0.01, base_width=15.5, base_height=5.5)


scatter.entropy.bl0.5 <- ggplot(data = NULL) +
    geom_point(data = scatter.repr.bl0.5, aes(x = true.entropy, y = entropy), size=1) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    geom_text(data = stats.bl0.5, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.6, y = 2.75, parse=TRUE, size=4) +
    xlab("True Entropy") + ylab("Predicted Entropy") +
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_entropy_bl0.5.pdf"), scatter.entropy.bl0.5, base_width=15.5, base_height=5.5)

scatter.entropy.bl0.01 <- ggplot(data = NULL) +
    geom_point(data = scatter.repr.bl0.01, aes(x = true.entropy, y = entropy), size=1) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.6, y =3.1, parse=TRUE, size=4) +
    xlab("True Entropy") + ylab("Predicted Entropy") +
    coord_cartesian(ylim=c(0,3.2)) +
    facet_grid(type~method)
save_plot(paste0(maintext_plot_directory, "scatter_entropy_bl0.01.pdf"), scatter.entropy.bl0.01, base_width=15.5, base_height=5.5)


theme_set(theme_cowplot() + theme(axis.text.x = element_text(size=10),
                                  axis.title.y = element_text(size=14)))

r2.dnds.bl0.5 <- sumstats %>% filter(bl==0.5) %>%
  ggplot(aes(x = method, y = r2.dnds, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  theme(legend.position = "none") +
  xlab("Inference Method") +ylab("Variance Explained")
r2.dnds.bl0.01 <- sumstats %>% filter(bl==0.01) %>%
  ggplot(aes(x = method, y = r2.dnds, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  theme(legend.position = "none") +
  xlab("Inference Method") + ylab("Variance Explained")

b.dnds.bl0.5 <- sumstats %>% filter(bl==0.5) %>%
  ggplot(aes(x = method, y = b.dnds, group = dataset, color = as.factor(type), shape = as.factor(b.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) +
  theme(legend.position = "none") +
  xlab("Inference Method") + ylab("Estimator Bias")
b.dnds.bl0.01 <- sumstats %>% filter(bl==0.01) %>%
  ggplot(aes(x = method, y = b.dnds, group = dataset, color = as.factor(type), shape = as.factor(b.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) +
  theme(legend.position = "none") +
  xlab("Inference Method") + ylab("Estimator Bias")

slope.dnds.bl0.5 <- sumstats %>% filter(bl==0.5) %>%
  ggplot(aes(x = method, y = slope.dnds, group = dataset, color = as.factor(type), shape = as.factor(slope.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=1) +
  theme(legend.position = "none") +
  xlab("Inference Method") + ylab("True-Inferred Slope")
slope.dnds.bl0.01 <- sumstats %>% filter(bl==0.01) %>%
  ggplot(aes(x = method, y = slope.dnds, group = dataset, color = as.factor(type), shape = as.factor(slope.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=1) +
  theme(legend.position = "none") +
  xlab("Inference Method") + ylab("True-Inferred Slope")

r2.b.slope.dnds.bl0.5 <- plot_grid(r2.dnds.bl0.5, b.dnds.bl0.5, slope.dnds.bl0.5, nrow=3, labels=c("A", "B", "C"), align="hv")
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_dnds_bl0.5.pdf"), r2.b.slope.dnds.bl0.5, base_width=5, base_height=8)
r2.b.slope.dnds.bl0.01 <- plot_grid(r2.dnds.bl0.01, b.dnds.bl0.01, slope.dnds.bl0.01, nrow=3, labels=c("A", "B", "C"), align="hv")
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_dnds_bl0.01.pdf"), r2.b.slope.dnds.bl0.01, base_width=5, base_height=8)



r2.entropy <- sumstats %>%
  ggplot(aes(x = method, y = r2.entropy, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  facet_grid(~bl) + theme(legend.position = "none") +
  xlab("Inference Method") + ylab("Variance Explained")
b.entropy <- sumstats %>%
  ggplot(aes(x = method, y = b.entropy, group = dataset, color = as.factor(type), shape = as.factor(b.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) +
  facet_grid(~bl) + theme(legend.position = "none") +
  xlab("Inference Method") + ylab("Estimator Bias")
slope.entropy <- sumstats %>%
  ggplot(aes(x = method, y = slope.entropy, group = dataset, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=1) +
  facet_grid(~bl) + theme(legend.position = "none") +
  xlab("Inference Method") + ylab("True-Inferred Slope")
r2.b.slope.entropy <- plot_grid(r2.entropy, b.entropy, slope.entropy, nrow=3, labels=c("A", "B", "C"))
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_entropy.pdf"), r2.b.slope.entropy, base_width=9, base_height=7.5)





#############################################################################################
############################# Selection coefficients histograms #############################

theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.2, 2, 0.2, 0.2),"cm")))
grid_list <- c()
i <- 1
for (d in sub.datasets){
  sc <- true.selcoeffs %>% filter(dataset == d) %>% mutate(method = "True")
  for (m in methods_levels){
    sc.temp <- read.csv(paste0("dataframes/",d,"_bl0.5_", m, "_selcoeffs.csv"))
    temp <- data.frame(dataset = d, binnedcoeffs = sc.temp$binnedcoeff, method = m)
    sc <- rbind(sc, temp)
  }
  sc$method <- factor(sc$method, levels = c("True", "nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes"), labels = c("True", "Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel"))
  x <- ggplot(sc, aes(x = binnedcoeffs)) + geom_histogram(fill = "white", color = "black") + facet_grid(~method) + xlab("Selection Coefficients") + ylab("Count")
  x2 <- ggdraw(x) + draw_label(d, x = 0.965, y = 0.5, size=14, fontface = "bold")
  grid_list[[i]] <- x2
  i <- i + 1
}

grid <- plot_grid(plotlist = grid_list, nrow=5)
save_plot(paste0(maintext_plot_directory, "yeast_dms_sc_grid.pdf"), grid, base_width = 12.5, base_height = 11)



##########################################################################################
##### Figure 4: JSD regressed on true dN/dS for representative, and slopes for all #######
##########################################################################################
dat %>% filter(dataset %in% repr.datasets, bl == 0.5) %>%
  ggplot(aes(x = true.dnds, y = jsd)) +
  geom_point(size=1) + facet_grid(dataset~method) +
  geom_smooth(method="lm", color="red") +
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1), labels=c('0.0', '0.25', '0.5', '0.75', '1.0'))+
  xlab("True dN/dS") + ylab("Site JSD") -> jsd.dnds.reprs.bl0.5
dat %>% filter(dataset %in% repr.datasets, bl == 0.01) %>%
  ggplot(aes(x = true.dnds, y = jsd)) +
  geom_point(size=1) + facet_grid(dataset~method) +
  geom_smooth(method="lm", color="red") +
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1), labels=c('0.0', '0.25', '0.5', '0.75', '1.0'))+
  xlab("True dN/dS") + ylab("Site JSD") -> jsd.dnds.reprs.bl0.01
save_plot(paste0(maintext_plot_directory, "jsd_dnds_scatter_bl0.5.pdf"), jsd.dnds.reprs.bl0.5, base_width = 13, base_height=3)
save_plot(paste0(maintext_plot_directory, "jsd_dnds_scatter_bl0.01.pdf"), jsd.dnds.reprs.bl0.01, base_width = 13, base_height=3)

dat %>% filter(dataset %in% repr.datasets, bl == 0.5) %>%
  ggplot(aes(x = true.entropy, y = jsd)) +
  geom_point(size=1) + facet_grid(dataset~method) +
  geom_smooth(method="lm", color="red") +
  xlab("True Entropy") + ylab("Site JSD") + ggtitle("BL = 0.5")-> jsd.entropy.reprs.bl0.5
dat %>% filter(dataset %in% repr.datasets, bl == 0.01) %>%
  ggplot(aes(x = true.entropy, y = jsd)) +
  geom_point(size=1) + facet_grid(dataset~method) +
  geom_smooth(method="lm", color="red") +
  xlab("True Entropy") + ylab("Site JSD") + ggtitle("BL = 0.01")-> jsd.entropy.reprs.bl0.01
jsd.entropy.scatters <- plot_grid(jsd.entropy.reprs.bl0.5, jsd.entropy.reprs.bl0.01, nrow=2, labels="AUTO")
save_plot(paste0(maintext_plot_directory, "jsd_entropy_scatter.pdf"), jsd.entropy.scatters, base_width = 13, base_height=7)


jsd.true.slope.p <- dat %>% group_by(dataset,method,bl,type) %>%
  do(fit=lm(jsd~true.dnds, data=.)) %>%
  mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>%
  select(-fit)

dnds.jsd.bl0.5 <- jsd.true.slope.p %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  scale_shape_manual(values=c(1,19)) +
  scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +
  theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope")

subset.jsd.true.slope.p <- dat %>% filter(true.dnds >= 0.3 & true.dnds <= 0.6) %>%
  group_by(dataset,method,bl,type) %>%
  do(fit=lm(jsd~true.dnds, data=.)) %>%
  mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>%
  select(-fit)

subset.dnds.jsd.bl0.5 <- subset.jsd.true.slope.p %>% filter(bl == 0.5) %>%
  ggplot(aes(x = method, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  scale_shape_manual(values=c(1,19)) +
  scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +
  theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope")

jsd.true.slope <- plot_grid(dnds.jsd.bl0.5, subset.dnds.jsd.bl0.5, nrow=1, labels="AUTO")
save_plot(paste0(maintext_plot_directory, "jsd_truednds_slopes.pdf"), jsd.true.slope, base_width = 10, base_height = 3)
#dnds.jsd.jitter.bl0.01 <- jsd.true.slope.p %>% filter(bl == 0.01) %>%
#  ggplot(aes(x = method, y = slope, shape=sig)) +
#  geom_point(position = position_jitter(w = 0.4)) +
#  scale_shape_manual(values=c(1,19)) +
#  geom_hline(yintercept=0) +
#  theme(legend.position="none") +
#  xlab("Inference Method") + ylab("Slope")

fig4 <- plot_grid(, dnds.jsd.jitter.bl0.5, nrow=2, labels=c("A", "B"))
save_plot(paste0(maintext_plot_directory, "jsd_dnds_scatter.pdf"), dnds.on.jsd.bl0.5, base_width=9.5, base_height=4)
