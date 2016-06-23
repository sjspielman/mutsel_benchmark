# SJS
# Creates MS figures and perform some stats

library(cowplot)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(grid)
library(multcomp)
library(lme4)
library(lmerTest)

# Function to test whether a slope for simple linear regression is significantly different from 1
test.slope <- function(y, x, compare=1)
{
  df <- length(y) - 2
  fit <- lm(y~x)
  b <- summary(fit)$coefficients[2]
  b.se <- summary(fit)$coefficients[4]
  t <- (b - compare)/b.se
  pvalue <- (1 - pt(t, df))*2
  list(pvalue,b)
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
    do(rraw.dnds  = lm(.$dnds ~ .$true.dnds, dat=.),
       braw.dnds  = glm(dnds ~ offset(true.dnds), dat=.),
       rraw.h     = lm(.$entropy ~ .$true.entropy, dat=.),
       braw.h     = glm(entropy ~ offset(true.entropy), dat=.))   %>%
    mutate(r2.dnds        = summary(rraw.dnds)$r.squared,
           r2.dnds.sig    = lmp(rraw.dnds)<siglevel,
           b.dnds         = summary(braw.dnds)$coeff[[1]],
           b.dnds.sig     = summary(braw.dnds)$coeff[4]<siglevel,
           r2.entropy     = summary(rraw.h)$r.squared,
           r2.entropy.sig = lmp(rraw.h)<siglevel,
           b.entropy      = summary(braw.h)$coeff[1],
           b.entropy.sig  = summary(braw.h)$coeff[4]<siglevel) %>%
    dplyr::select(-rraw.dnds, -rraw.h, -braw.dnds, -braw.h) -> dat.sum1
  dat %>% group_by(dataset, method, bl) %>%
    do(x = test.slope(.$dnds, .$true.dnds)) %>%
    mutate(slope.dnds.sig = x[[1]]<siglevel, slope.dnds = x[[2]]) %>%
    dplyr::select(-x) -> dat.sum2
  dat %>% group_by(dataset, method, bl) %>%
    do(x = test.slope(.$entropy, .$true.entropy)) %>%
    mutate(slope.entropy.sig = x[[1]]<siglevel, slope.entropy = x[[2]]) %>%
    dplyr::select(-x) -> dat.sum3
  part <- left_join(dat.sum2, dat.sum3)
  dat.sum <- left_join(dat.sum1, part)
  dat.sum
}

theme_set(theme_cowplot() + theme(panel.margin = unit(0.9, "lines")))

result_directory         <- "dataframes/"
true_directory           <- "../simulation/true_simulation_parameters/"
maintext_plot_directory  <- "maintext_figures/"
si_plot_directory        <- "SI_figures/"
repr.datasets            <- c("1R6M_A", "NP")
dms.datasets             <- c("Gal4", "LAC", "HA", "NP")
methods_levels           <- c("nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes")
methods_labels           <- c("Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel")
alpha                    <- 0.01 # Significance
corrected.alpha          <- alpha/length(methods_levels) #Bonferroni significance


dat <- read.csv(paste0(result_directory, "inference_results.csv"))
true.selcoeffs <- read_csv(paste0(result_directory, "true_selection_coefficients.csv"))
datasets <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A", "Gal4", "LAC", "HA", "NP")
sub.datasets <- c("1R6M_A", "Gal4", "LAC", "HA", "NP")


dat <- mutate(dat, type = ifelse(!(dataset %in% dms.datasets), "Natural", "DMS"))
dat$method <- factor(dat$method, levels = methods_levels, labels = methods_labels)
dat$bl <- factor(dat$bl, levels = c(0.5, 0.01))
dat$type <- factor(dat$type, levels=c("Natural", "DMS"))
dat$dataset <- factor(dat$dataset, levels=datasets)
type.colors <- c("red", "blue")
sumstats <- dat %>% summarize_dnds_entropy(corrected.alpha)

#### Shared "type" legend
fake <- data.frame(x = seq(1,4), y = seq(1,4), type = c("Natural", "Natural", "DMS", "DMS"))
fake$type <- factor(fake$type, levels=c("Natural", "DMS"))
fake %>%
  ggplot(aes(x = y, y = y, color = as.factor(type))) +
  geom_point(size=2.5) + geom_line(size=0.75) +
  scale_color_manual(values = type.colors, name = "Data type ") +
  theme(legend.position = "bottom") -> type.legend.raw
grobs <- ggplotGrob(type.legend.raw)$grobs
type.legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]




##################### Distributions of true selection coefficients and selection pressures ( meaning dN/dS, entropy) ######################
true.selcoeffs %>% filter(dataset %in% sub.datasets) -> sub.selcoeffs
sub.selcoeffs <- mutate(sub.selcoeffs, type = ifelse(!(dataset %in% dms.datasets), "Natural", "DMS"))
sub.selcoeffs$type <- factor(sub.selcoeffs$type, levels=c("Natural", "DMS"))
sub.selcoeffs$dataset <- factor(sub.selcoeffs$dataset, levels=sub.datasets)

dat %>% filter(dataset %in% sub.datasets, method == "mvn10") -> sub.dat
sub.dat$dataset <- factor(sub.dat$dataset, levels=sub.datasets)

compare.truesc <- sub.selcoeffs %>% ggplot(aes(x = binnedcoeffs)) + geom_density(aes(fill = type), alpha=0.7) + facet_grid(~dataset, scales="free_y") + scale_x_continuous(limits=c(-10,10)) + ylab("Density")  + xlab("True S") + scale_fill_manual(name = "Data Type", values = type.colors) + theme(legend.title = element_text(size=12), legend.text = element_text(size=10), legend.key.size = unit(0.4, "cm"))
save_plot(paste0(maintext_plot_directory, "small_selcoeffs.pdf"), compare.truesc, base_width=10, base_height=2)



compare.truednds <- dat %>% ggplot(aes(x = dataset, y = true.dnds)) + geom_boxplot(aes(fill = type)) + xlab("Dataset") + ylab("True dN/dS") + scale_fill_manual(name = "Data Type", values = type.colors) + scale_x_discrete(limits=datasets) +
theme(legend.position = "none")
compare.trueh <- dat %>% ggplot(aes(x = dataset, y = true.entropy)) + geom_boxplot(aes(fill = type)) + xlab("Dataset") + ylab("True Entropy") + scale_fill_manual(name = "Data Type", values = type.colors) +  scale_x_discrete(limits=datasets) + theme(legend.position = "none")
compare.dnds.h <- plot_grid(compare.truednds, compare.trueh, ncol = 1, labels = "AUTO", align = "hv")
compare.grid.legend <- plot_grid(compare.dnds.h, type.legend, ncol = 1, rel_heights = c(1,0.05))
save_plot(paste0(si_plot_directory, "full_boxplot_dnds_entropy.pdf"), compare.grid.legend, base_width=12, base_height=5.5)

###################### JSD and diff sum #######################

dat %>% group_by(dataset, type, method, bl) %>% summarize(meanjsd = mean(jsd), meandiffsum = mean(diffsum)) %>% mutate(x=as.numeric(method)) -> jsd.diffsum

### JSD statistics:
dat %>% filter(bl == 0.5) -> dat1
dat %>% filter(bl == 0.01) -> dat2

fit1 <- lmer(jsd ~ method + type + (1|dataset), data = dat1)
fit1.mc.method <- glht(fit1, linfct=mcp(method='Tukey')) # JSD among methods
summary(fit1.mc.method)


jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = meanjsd, color = type, group = dataset)) +
  geom_point(size=2) + geom_line(size=0.75, alpha=0.5) + ggtitle("BL = 0.5") +
  geom_text_repel(data = subset(jsd.diffsum, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_color_manual(values=type.colors, name = "Data type") + ylab("Average JSD")+
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(jsd.diffsum$x)), labels = levels(jsd.diffsum$method),limits = c(1, max(jsd.diffsum$x) + 1)) +
  theme(legend.position = "none") -> jsd.lineplot.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = meanjsd, color = as.factor(type), group = dataset)) +
  geom_point(size=2 ) + geom_line(size=0.75, alpha=0.5) + ggtitle("BL = 0.01") +
  geom_text_repel(data = subset(jsd.diffsum, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_color_manual(values=type.colors, name = "Data type") + ylab("Average JSD")+
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(jsd.diffsum$x)), labels = levels(jsd.diffsum$method),limits = c(1, max(jsd.diffsum$x) + 1))  +
  theme(legend.position = "none") -> jsd.lineplot.raw2
jsd.lineplot1 <- plot_grid(jsd.lineplot.raw1, jsd.lineplot.raw2, nrow=1, labels="AUTO")
jsd.lineplot  <- plot_grid(jsd.lineplot1, type.legend, ncol=1, rel_heights=c(1, 0.09))
save_plot(paste0(maintext_plot_directory, "jsd_lineplot.pdf"), jsd.lineplot, base_width = 12, base_height=3.5)

jsd.diffsum %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = meandiffsum, color = as.factor(type), group = dataset)) +
  geom_point(size=2 ) + geom_line(size=0.75, alpha=0.5) + ggtitle("BL = 0.5") +
  geom_text_repel(data = subset(jsd.diffsum, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_color_manual(values=type.colors, name = "Data type") + ylab("Average Sum of Differences")+
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(jsd.diffsum$x)), labels = levels(jsd.diffsum$method),limits = c(1, max(jsd.diffsum$x) + 1)) +
  theme(legend.position = "none") -> diffsum.lineplot.raw1
jsd.diffsum %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = meandiffsum, color = as.factor(type), group = dataset)) +
  geom_point(size=2 ) + geom_line(size=0.75, alpha=0.5) + ggtitle("BL = 0.01") +
  geom_text_repel(data = subset(jsd.diffsum, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_color_manual(values=type.colors, name = "Data type") + ylab("Average Sum of Differences")+
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(jsd.diffsum$x)), labels = levels(jsd.diffsum$method),limits = c(1, max(jsd.diffsum$x) + 1))  +
  theme(legend.position = "none") -> diffsum.lineplot.raw2
diffsum.lineplot1 <- plot_grid(diffsum.lineplot.raw1, diffsum.lineplot.raw2, nrow=1, labels="AUTO")
diffsum.lineplot  <- plot_grid(diffsum.lineplot1, type.legend, ncol=1, rel_heights=c(1, 0.09))
save_plot(paste0(si_plot_directory, "diffsum_lineplot.pdf"), diffsum.lineplot, base_width = 12, base_height=3.5)


###################### dN/dS and entropy scatterplots #######################

sumstats %>%
    dplyr::select(-b.dnds, -b.dnds.sig, -b.entropy, -b.entropy.sig, -slope.dnds.sig, -slope.dnds, -slope.entropy.sig, -slope.entropy) %>%
    filter(dataset %in% repr.datasets) %>% group_by(dataset, method) %>%
    mutate(r2.dnds.short = as.numeric(round(r2.dnds, 2)), r2.entropy.short = as.numeric(round(r2.entropy, 2))) -> scatter.stats
stats.bl0.5 <- scatter.stats %>% filter(bl == 0.5)
stats.bl0.01 <- scatter.stats %>% filter(bl == 0.01)
dat %>% filter(dataset %in% repr.datasets, bl == 0.5)  -> scatter.repr.bl0.5
dat %>% filter(dataset %in% repr.datasets, bl == 0.01) -> scatter.repr.bl0.01
r2.size <- 3.5

scatter.theme <- theme(panel.margin = unit(0.5, "lines"),
                       legend.position = "none",
                       axis.text  = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       strip.text = element_text(size = 10.5),
                       strip.background = element_rect(size = 0.5))

scatter.dnds.bl0.5 <- ggplot(data = NULL) +
  geom_point(data = scatter.repr.bl0.5, aes(x = true.dnds, y = dnds, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  geom_text(data = stats.bl0.5, aes(label=paste0("r^2==",r2.dnds.short)), x = 0.25, y = 0.95, parse=TRUE, size=r2.size) +
  xlab("True dN/dS") + ylab("Predicted dN/dS") + facet_grid(dataset~method)+ scale_color_manual(values=type.colors, name = "Dataset") + scatter.theme + scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1)) + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))
scatter.entropy.bl0.5 <- ggplot(data = NULL) +
    geom_point(data = scatter.repr.bl0.5, aes(x = true.entropy, y = entropy, color = as.factor(type)), size=1) +
    geom_abline(slope = 1, intercept = 0, size=0.75) +
    geom_text(data = stats.bl0.5, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.75, y = 2.75, parse=TRUE, size=r2.size) +
    xlab("True Entropy") + ylab("Predicted Entropy") + scatter.theme +
    facet_grid(dataset~method) + scale_color_manual(values=type.colors, name = "Dataset") + coord_cartesian(xlim=c(0,3.1), ylim=c(0,3.1))
#scatter.bl0.5.grid <- plot_grid(scatter.dnds.bl0.5, scatter.entropy.bl0.5, nrow=2, labels="AUTO", align ="v")
#save_plot(paste0(maintext_plot_directory, "scatter_dnds_entropy_bl0.5.pdf"), scatter.bl0.5.grid, base_width=10, base_height=5.75)

scatter.dnds.bl0.01 <- ggplot(data = NULL) +
  geom_point(data = scatter.repr.bl0.01, aes(x = true.dnds, y = dnds, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.dnds.short)), x = 0.3, y = 1.12, parse=TRUE, size=r2.size) +
  xlab("True dN/dS") + ylab("Predicted dN/dS") +
  scale_x_continuous(limits=c(0,1.1), breaks=c(0,0.5,1)) + scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.5,1)) +
  facet_grid(dataset~method) + scale_color_manual(values=type.colors, name = "Dataset")+ scatter.theme
scatter.entropy.bl0.01 <- ggplot(data = NULL) +
  geom_point(data = scatter.repr.bl0.01, aes(x = true.entropy, y = entropy, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  geom_text(data = stats.bl0.01, aes(label=paste0("r^2==",r2.entropy.short)), x = 0.67, y =3, parse=TRUE, size=r2.size-0.25) +
  xlab("True Entropy") + ylab("Predicted Entropy") + scatter.theme +
  coord_cartesian(xlim=c(0,3.1), ylim=c(0,3.1)) + scale_color_manual(values=type.colors, name = "Dataset") + facet_grid(dataset~method)
#scatter.bl0.01.grid <- plot_grid(scatter.dnds.bl0.01, scatter.entropy.bl0.01, nrow=2, labels="AUTO")
#save_plot(paste0(si_plot_directory, "scatter_dnds_entropy_bl0.01.pdf"), scatter.bl0.01.grid, base_width=11.5,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                base_height=6)


scatter.grid <- plot_grid(scatter.dnds.bl0.5, scatter.entropy.bl0.5, scatter.dnds.bl0.01, scatter.entropy.bl0.01, nrow=2, labels="AUTO", align ="v")
save_plot(paste0(maintext_plot_directory, "scatter_dnds_entropy_full.pdf"), scatter.grid, base_width=15.5, base_height=5.5)











##########################################################################################
################################# r^2, bias, and slope plots #############################
##########################################################################################

this.theme <- theme_cowplot() + theme(axis.text.x = element_text(size=9.5),
                                      axis.text.y = element_text(size = 10.5),
                                      axis.title  = element_text(size=12),
                                      legend.position = "none")
sumstats2 <- mutate(sumstats, x=as.numeric(method)) # manually create x values
sumstats2$type <- factor(sumstats2$type, levels=c("Natural", "DMS"))

r2.dnds.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = r2.dnds, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x = 0.5, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + this.theme + scale_y_continuous(limits=c(0.1, 1)) +
  scale_color_manual(values=type.colors)  + ylab("Variance Explained") + ggtitle("dN/dS\n")
r2.h.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = r2.entropy, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x=0.5, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + scale_y_continuous(name = "Variance Explained", limits=c(0.1, 1)) + ggtitle("Entropy\n")+ scale_color_manual(values=type.colors) + this.theme

b.dnds.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = b.dnds, group = dataset, color = as.factor(type), shape = as.factor(b.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x = 0.5, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + this.theme + ylab("Estimator Bias")
b.h.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = b.entropy, group = dataset, color = as.factor(type), shape = as.factor(b.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x = 0.5, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + this.theme  + ylab("Estimator Bias")

slope.dnds.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = slope.dnds, group = dataset, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x = 0.5, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + scale_y_continuous(name = "True-Inferred Slope", limits=c(0.2, 1.1)) + geom_hline(yintercept=1) + this.theme
slope.h.bl0.5 <- sumstats2 %>% filter(bl==0.5) %>%
  ggplot(aes(x = x, y = slope.entropy, group = dataset, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.5 & type == "DMS"), nudge_x = 0.5, x=6, size=3, aes(label = dataset)) + scale_y_continuous(name = "True-Inferred Slope", limits=c(0.2, 1.1)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + geom_hline(yintercept=1) + this.theme
dnds.grid.bl0.5 <- plot_grid(r2.dnds.bl0.5, b.dnds.bl0.5, slope.dnds.bl0.5, nrow=3, labels=c("A", "B", "C"), align="hv")
entropy.grid.bl0.5 <- plot_grid(r2.h.bl0.5, b.h.bl0.5, slope.h.bl0.5, nrow=3, labels=c("D", "E", "F"), align="hv")
stats.grid.bl0.5 <- plot_grid(dnds.grid.bl0.5, entropy.grid.bl0.5, ncol=2, align = "h")
save_plot(paste0(maintext_plot_directory, "r2_bias_slope_bl0.5.pdf"), stats.grid.bl0.5, base_width=9.75, base_height=7)






r2.dnds.bl0.01 <- sumstats2 %>% filter(bl==0.01) %>%
  ggplot(aes(x = x, y = r2.dnds, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.01 & type == "DMS"), nudge_x = 0.01, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + this.theme + scale_y_continuous(limits=c(0, 0.7)) +
  scale_color_manual(values=type.colors) + ggtitle("dN/dS\n") + ylab("Variance Explained")
r2.h.bl0.01 <- sumstats2 %>% filter(bl==0.01) %>%
  ggplot(aes(x = x, y = r2.entropy, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.01 & type == "DMS"), nudge_x=0.01, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) + scale_y_continuous(limits=c(0, 0.7)) +
  scale_color_manual(values=type.colors) + this.theme  + ggtitle("Entropy\n") + ylab("Variance Explained")

b.dnds.bl0.01 <- sumstats2 %>% filter(bl==0.01) %>%
  ggplot(aes(x = x, y = b.dnds, group = dataset, color = as.factor(type), shape = as.factor(b.dnds.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.01 & type == "DMS"), nudge_x = 0.01, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + this.theme + ylab("Estimator Bias")
b.h.bl0.01 <- sumstats2 %>% filter(bl==0.01) %>%
  ggplot(aes(x = x, y = b.entropy, group = dataset, color = as.factor(type), shape = as.factor(b.entropy.sig))) +
  geom_point(size=2) + geom_line(alpha=0.5) +
  geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & bl==0.01 & type == "DMS"), nudge_x = 0.01, x=6, size=3, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.6)) +
  scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept=0) + this.theme  + ylab("Estimator Bias")




















##########################################################################################
#########  JSD regressed on true dN/dS, H for representative, and slopes for all #########
##########################################################################################
full <- dat %>% group_by(dataset,method,bl,type) %>%
  do(fit=lm(jsd~true.dnds, data=.)) %>%
  mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>%
  dplyr::select(-fit)
subset <- dat %>% filter(true.dnds >= 0.3 & true.dnds <= 0.6) %>%
  group_by(dataset,method,bl,type) %>%
  do(fit=lm(jsd~true.dnds, data=.)) %>%
  mutate(slope = round(fit[[1]][[2]],3), pvalue = lmp(fit), sig = pvalue < corrected.alpha) %>%
  dplyr::select(-fit)

full2 <- mutate(full, x=as.numeric(method))
subset2 <- mutate(subset, x=as.numeric(method))

theme_set(theme_cowplot() + theme(panel.margin = unit(0.6, "lines")))

dat %>% filter(dataset %in% repr.datasets, bl == 0.5) %>%
  ggplot(aes(x = true.dnds, y = jsd, color = type)) +
  geom_point(size=1) + facet_grid(dataset~method) +
  geom_smooth(method="lm", color = "black") + scale_color_manual(values = type.colors) +
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1), labels=c('0.0', '0.25', '0.5', '0.75', '1.0'))+
  xlab("True dN/dS") + ylab("Site JSD") + theme(legend.position = "none") -> jsd.dnds.reprs.bl0.5
dnds.jsd.bl0.5 <- full2 %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(full2, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(full2$x)), labels = levels(full2$method),limits = c(1, max(full2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) + theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("Full")+ scale_y_continuous(limits=c(-0.5, 0.35))
subset.dnds.jsd.bl0.5 <- subset2 %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(subset2, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.25, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(subset2$x)), labels = levels(subset2$method),limits = c(1, max(subset2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +  theme(legend.position="none") + scale_y_continuous(limits=c(-0.5, 0.35))+
  xlab("Inference Method") + ylab("Slope") + ggtitle("Subset")
jsd.dnds.true.slope <- ggdraw() + draw_plot(jsd.dnds.reprs.bl0.5, 0, 0.5, 1, 0.5) +
                                  draw_plot(dnds.jsd.bl0.5, 0, 0, 0.5, 0.5) +
                                  draw_plot(subset.dnds.jsd.bl0.5, 0.5, 0, 0.5, 0.5) +
                                  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

save_plot(paste0(maintext_plot_directory, "jsd_truednds_slopes.pdf"), jsd.dnds.true.slope , base_width = 11.5, base_height = 5)




entropy.jsd.bl0.5 <- full2 %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(full2, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(full2$x)), labels = levels(full2$method),limits = c(1, max(full2$x) + 0.35)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) + theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("Full")
subset.entropy.jsd.bl0.5 <- subset2 %>% filter(bl == 0.5) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(subset2, method == "pbMutSel" & bl == 0.5 & type == "DMS"), nudge_x = 0.25, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(subset2$x)), labels = levels(subset2$method),limits = c(1, max(subset2$x) + 0.35)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +  theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("Subset")

jsd.true.slope <- plot_grid(entropy.jsd.bl0.5, subset.entropy.jsd.bl0.5, nrow=2, labels="AUTO", align = "hv")
save_plot(paste0(si_plot_directory, "jsd_trueentropy_slopes.pdf"), jsd.true.slope, base_width = 7, base_height = 6)




dnds.jsd.bl0.01 <- full2 %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(full2, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(full2$x)), labels = levels(full2$method),limits = c(1, max(full2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) + theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("BL = 0.01")
subset.dnds.jsd.bl0.01 <- subset2 %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(subset2, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.25, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(subset2$x)), labels = levels(subset2$method),limits = c(1, max(subset2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +  theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("Subset")
jsd.true.slope <- plot_grid(dnds.jsd.bl0.01, subset.dnds.jsd.bl0.01, nrow=2, labels="AUTO", align = "hv")
save_plot(paste0(si_plot_directory, "jsd_truednds_slopes_bl0.01.pdf"), jsd.true.slope, base_width = 7, base_height = 6)



entropy.jsd.bl0.01 <- full2 %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(full2, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.5, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(full2$x)), labels = levels(full2$method),limits = c(1, max(full2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) + theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("BL = 0.01")
subset.entropy.jsd.bl0.01 <- subset2 %>% filter(bl == 0.01) %>%
  ggplot(aes(x = x, y = slope, shape=sig, group = dataset, color = as.factor(type))) +
  geom_point(size=2) + geom_line(alpha = 0.5) +
  geom_text_repel(data = subset(subset2, method == "pbMutSel" & bl == 0.01 & type == "DMS"), nudge_x = 0.25, x=6.1, size=4, aes(label = dataset)) +
  scale_x_continuous(name = "Inference Method", breaks = sort(unique(subset2$x)), labels = levels(subset2$method),limits = c(1, max(subset2$x) + 0.5)) +
  scale_shape_manual(values=c(1,19)) + scale_color_manual(values = type.colors) +
  geom_hline(yintercept=0) +  theme(legend.position="none") +
  xlab("Inference Method") + ylab("Slope") + ggtitle("Subset")
jsd.true.slope <- plot_grid(entropy.jsd.bl0.01, subset.entropy.jsd.bl0.01, nrow=2, labels="AUTO", align = "hv")
save_plot(paste0(si_plot_directory, "jsd_trueentropy_slopes_bl0.01.pdf"), jsd.true.slope, base_width = 7, base_height = 6)








#############################################################################################
##################################### All the scatterplots ##################################

dat %>% filter(!(dataset %in% repr.datasets), bl == 0.5)  -> scatter.bl0.5
dat %>% filter(!(dataset %in% repr.datasets), bl == 0.01) -> scatter.bl0.01

scatter.theme <- theme_cowplot() + theme(legend.position = "none",
                                         strip.text = element_text(size=10),
                                         axis.text  = element_text(size=10))

scatter.dnds.bl0.5 <- ggplot(scatter.bl0.5) +
  geom_point(aes(x = true.dnds, y = dnds, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  xlab("True dN/dS") + ylab("Predicted dN/dS") + facet_grid(dataset~method)+ scale_color_manual(values=type.colors, name = "Dataset") + scatter.theme + scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1)) + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))
save_plot(paste0(si_plot_directory, "scatter_dnds_bl0.5_full.pdf"), scatter.dnds.bl0.5 , base_width = 8, base_height = 10.5)

scatter.dnds.bl0.01 <- ggplot(scatter.bl0.01) +
  geom_point(aes(x = true.dnds, y = dnds, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  xlab("True dN/dS") + ylab("Predicted dN/dS") + facet_grid(dataset~method)+ scale_color_manual(values=type.colors, name = "Dataset") + scatter.theme + scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1)) + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))
save_plot(paste0(si_plot_directory, "scatter_dnds_bl0.01_full.pdf"), scatter.dnds.bl0.01, base_width = 8, base_height = 10.5)

scatter.entropy.bl0.5 <- ggplot(scatter.bl0.5) +
  geom_point(aes(x = true.entropy, y = entropy, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  xlab("True Entropy") + ylab("Predicted Entropy") + facet_grid(dataset~method)+ scale_color_manual(values=type.colors, name = "Dataset") + scatter.theme
save_plot(paste0(si_plot_directory, "scatter_entropy_bl0.5_full.pdf"), scatter.entropy.bl0.5 , base_width = 8, base_height = 10.5)

scatter.entropy.bl0.01 <- ggplot(scatter.bl0.01) +
  geom_point(aes(x = true.entropy, y = entropy, color = as.factor(type)), size=1) +
  geom_abline(slope = 1, intercept = 0, size=0.75) +
  xlab("True Entropy") + ylab("Predicted Entropy") + facet_grid(dataset~method)+ scale_color_manual(values=type.colors, name = "Dataset") + scatter.theme
save_plot(paste0(si_plot_directory, "scatter_entropy_bl0.01_full.pdf"), scatter.entropy.bl0.01, base_width = 8, base_height = 10.5)












#############################################################################################
############################# Selection coefficients histograms #############################

# Plot a subset of histograms, place the remaining structural in SI
theme_set(theme_cowplot() + theme(plot.margin = unit(c(0.2, 2, 0.2, 0.2),"cm"),
                                  legend.position = "none",
                                  strip.text = element_text(size=10)))
grid_list <- c()
i <- 1
for (d in sub.datasets){
  sc <- true.selcoeffs %>% filter(dataset == d) %>% mutate(method = "True")
  for (m in methods_levels){
    sc.temp <- read.csv(paste0("dataframes/",d,"_bl0.5_", m, "_selcoeffs.csv"))
    temp <- data.frame(dataset = d, binnedcoeffs = sc.temp$binnedcoeff, method = m)
    sc <- rbind(sc, temp)
  }
  if(d %in% dms.datasets){
    col <- "blue"
  }
  else
  {
  col <- "red"
  }
  sc$method <- factor(sc$method, levels = c("True", "nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes"), labels = c("True", "Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel"))
  x <- ggplot(sc, aes(x = binnedcoeffs)) + geom_histogram(fill = col, color = col) + facet_grid(~method) + xlab("Selection Coefficients") + ylab("Count")
  x2 <- ggdraw(x) + draw_label(d, x = 0.955, y = 0.5, size=14, fontface = "bold")
  grid_list[[i]] <- x2
  i <- i + 1
}
grid <- plot_grid(plotlist = grid_list, nrow=5)
save_plot(paste0(maintext_plot_directory, "selcoeff_histograms_subset.pdf"), grid, base_width = 10.5, base_height = 9)


# bl0.5, everything except sub
theme_set(theme_cowplot() + theme(plot.margin = unit(c(0, 1.5, 0, 0.01),"cm"),  # trbl
                                  axis.text = element_text(size = 7),
                                  axis.title = element_text(size = 8),
                                  strip.text = element_blank(),
                                  strip.background = element_blank()))
label_size = 8
grid_list <- c()
i <- 1
for (d in datasets){
  if (d %in% sub.datasets) next
  sc <- true.selcoeffs %>% filter(dataset == d) %>% mutate(method = "True")
  for (m in methods_levels){
    sc.temp <- read.csv(paste0("dataframes/",d,"_bl0.5_", m, "_selcoeffs.csv"))
    temp <- data.frame(dataset = d, binnedcoeffs = sc.temp$binnedcoeff, method = m)
    sc <- rbind(sc, temp)
  }
  sc$method <- factor(sc$method, levels = c("True", "nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes"), labels = c("True", "Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel"))
  x <- ggplot(sc, aes(x = binnedcoeffs)) + geom_histogram(fill = "grey40", color = "grey40") + facet_grid(~method) + xlab("") + ylab("Count")
  x2 <- ggdraw(x) + draw_label(d, x = 0.96, y = 0.6, size=label_size, fontface = "bold")
  grid_list[[i]] <- x2
  i <- i + 1
}
grid <- plot_grid(plotlist = grid_list, nrow=i-1)
title <- ggdraw() + draw_label(paste("          True", "Unpenalized", "mvn100  ", "mvn10   ", "d0.01   ", "d0.1   ", "pbMutSel      ", sep = "              "), size = label_size, fontface="bold")
base <- ggdraw() + draw_label("Selection Coefficients", size = label_size)
grid.with.title <- plot_grid(title, grid, base, ncol=1, rel_heights=c(0.03, 1, 0.02))
save_plot(paste0(si_plot_directory, "selcoeff_histograms_bl0.5.pdf"), grid.with.title, base_width = 7, base_height = 9)


# bl0.01, all
dms.grid_list <- c()
natural.grid_list <- c()
i <- 1
j <- 1
for (d in datasets){
  sc <- true.selcoeffs %>% filter(dataset == d) %>% mutate(method = "True")
  for (m in methods_levels){
    sc.temp <- read.csv(paste0("dataframes/",d,"_bl0.01_", m, "_selcoeffs.csv"))
    temp <- data.frame(dataset = d, binnedcoeffs = sc.temp$binnedcoeff, method = m)
    sc <- rbind(sc, temp)
  }
  sc$method <- factor(sc$method, levels = c("True", "nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes"), labels = c("True", "Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel"))
  x <- ggplot(sc, aes(x = binnedcoeffs)) + geom_histogram(fill = "grey40", color = "grey40") + facet_grid(~method) + ylab("Count") + xlab("")
  x2 <- ggdraw(x) + draw_label(d, x = 0.96, y = 0.7, size=label_size, fontface = "bold")

  if (d %in% dms.datasets){
    dms.grid_list[[i]] <- x2
    i <- i + 1
  }else{
    natural.grid_list[[j]] <- x2
    j <- j + 1
  }
}
dms.grid <- plot_grid(plotlist = dms.grid_list, nrow=i-1)
natural.grid <- plot_grid(plotlist = natural.grid_list, nrow = j-1)
title <- ggdraw() + draw_label(paste("          True", "Unpenalized", "mvn100  ", "mvn10   ", "d0.01   ", "d0.1   ", "pbMutSel      ", sep = "              "), size = label_size, fontface="bold")
base <- ggdraw() + draw_label("Selection Coefficients", size = label_size)
dms.grid.with.title <- plot_grid(title, dms.grid, base, ncol=1, rel_heights=c(0.06, 1, 0.06))
save_plot(paste0(si_plot_directory, "selcoeff_histograms_dms_bl0.01.pdf"), dms.grid.with.title, base_width = 7, base_height = 4)
natural.grid.with.title <- plot_grid(title, natural.grid, base, ncol=1, rel_heights=c(0.03, 1, 0.02))
save_plot(paste0(si_plot_directory, "selcoeff_histograms_natural_bl0.01.pdf"), natural.grid.with.title, base_width = 7, base_height = 9)















































# r2.dnds.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = r2.dnds, shape = r2.entropy.sig, group = dataset, color = type)) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6., size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.7)) + facet_grid(~bl, scales = "free_y") +
#   scale_color_manual(values=type.colors)  + scale_shape_manual(values=c(1,19)) + theme(legend.position = "none") + ylab("Variance Explained") + this.theme
# b.dnds.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = b.dnds, group = dataset, color = type, shape = as.factor(b.dnds.sig))) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6, size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.7)) +
#   scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) + facet_grid(~bl, scales = "free_y") +
#   geom_hline(yintercept=0) + theme(legend.position = "none") + ylab("Estimator Bias") + this.theme
# slope.dnds.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = slope.dnds, group = dataset, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) + facet_grid(~bl, scales = "free_y") +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6, size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.7)) +
#   geom_hline(yintercept=0.75) + theme(legend.position = "none") + ylab("True-Inferred Slope") + this.theme
# dnds.grid <- plot_grid(r2.dnds.bl0.5, b.dnds.bl0.5, slope.dnds.bl0.5, nrow=3, labels=c("A", "B", "C"), align="hv")
# dnds.grid.legend <- plot_grid(dnds.grid, type.legend, ncol=1, rel_heights=c(1,0.05))
# save_plot(paste0(maintext_plot_directory, "r2_bias_slope_dnds.pdf"), dnds.grid.legend, base_width=6, base_height=8)
#
#
#
#
# r2.entropy.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = r2.entropy, shape = r2.entropy.sig, group = dataset, color = type)) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6., size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.75)) + facet_grid(~bl, scales = "free_y") +
#   scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) +
#   theme(legend.position = "none") + ylab("Variance Explained")
# b.entropy.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = b.entropy, group = dataset, color = type, shape = as.factor(b.entropy.sig))) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6, size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.75)) +
#   scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) + facet_grid(~bl, scales = "free_y") +
#   geom_hline(yintercept=0) + theme(legend.position = "none") + ylab("Estimator Bias")
# slope.entropy.bl0.5 <- sumstats2 %>%
#   ggplot(aes(x = x, y = slope.entropy, group = dataset, color = as.factor(type), shape = as.factor(slope.entropy.sig))) +
#   geom_point(size=2) + geom_line(size = 0.75, alpha=0.5) +
#   scale_color_manual(values=type.colors) + scale_shape_manual(values=c(1,19)) + facet_grid(~bl, scales = "free_y") +
#   geom_text_repel(data = subset(sumstats2, method == "pbMutSel" & type == "DMS"), nudge_x = 0.6, x=6, size=3, aes(label = dataset)) +
#   scale_x_continuous(name = "Inference Method", breaks = sort(unique(sumstats2$x)), labels = levels(sumstats2$method),limits = c(1, max(sumstats2$x) + 0.75)) +
#   geom_hline(yintercept=0.75) + theme(legend.position = "none") + ylab("True-Inferred Slope")
# entropy.grid <- plot_grid(r2.entropy.bl0.5, b.entropy.bl0.5, slope.entropy.bl0.5, nrow=3, labels=c("A", "B", "C"), align="hv")
# entropy.grid.legend <- plot_grid(entropy.grid, type.legend, ncol=1, rel_heights=c(1,0.05))
# save_plot(paste0(si_plot_directory, "r2_bias_slope_entropy.pdf"), entropy.grid.legend, base_width=9.25, base_height=8)
