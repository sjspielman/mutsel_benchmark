require(cowplot)
require(dplyr)
require(tidyr)
require(readr)

sim.dat <- read.csv("../../results/summarized_results/simulation_derived_dnds.csv")
sim.dat %>% spread(method,dnds) %>% gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% select(dataset, del, site, true, dnds, method)-> dnds
sim.dat$method <- factor(sim.dat$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))

dnds.strong <- dnds %>% filter(dataset == "1B4T_A", del == "strong")
dnds.strong$method <- factor(dnds.strong$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))
dnds.weak <- dnds %>% filter(dataset == "1B4T_A", del == "weak")
dnds.weak$method <- factor(dnds.weak$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))


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
