require(cowplot)
require(dplyr)
require(tidyr)
require(readr)

dnds_raw <- read.csv("../../results/summarized_results/empirical_derived_dnds.csv")
dnds_raw %>% spread(method,dnds) -> dnds

for (name in emp.datasets){
    print(name)
    sc <- read_csv(paste0("../../results/", name, "_selection_coefficients.csv"))
    sc$method <- factor(sc$method, levels=c("nopenal", "phylobayes", "mvn1", "mvn10", "mvn100", "d1.0", "d0.1", "d0.01"))

    a <- dnds %>% filter(dataset == name) %>% ggplot(aes(x = nopenal, y = phylobayes)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0)
    b <- sc %>% filter(method %in% c("nopenal", "phylobayes")) %>% ggplot(aes(x=binnedcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
    p <- plot_grid(a,b,nrow=2)
    print(p)
    readline()
}    

    
#     a <- dnds %>% filter(dataset == dataset) %>% ggplot(aes(x = nopenal, y = d0.01)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d0.01")
#     b <- dnds %>% filter(dataset == dataset) %>%ggplot(aes(x = nopenal, y = d0.1)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d0.1")
#     cc <- dnds %>% filter(dataset == dataset) %>%ggplot(aes(x = nopenal, y = d1.0)) + geom_point(alpha=0.7) + geom_abline(slope=1, intercept=0) + ggtitle("d1.0")
#     
# 
#     d <- sc %>% filter(method %in% c("nopenal", "d0.01")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
#     e <- sc %>% filter(method %in% c("nopenal", "d0.1")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
#     f <- sc %>% filter(method %in% c("nopenal", "d1.0")) %>% ggplot(aes(x=realcoeff, fill=method)) + geom_density(alpha = 0.5) + scale_fill_manual(values=c("navy", "yellow"))
    
    plotgrid <- plot_grid(a,b,cc,d,e,f,nrow=2)
    print(plotgrid)
    readline()
}




# Some possibly useful code, REQUIRES PHYLOBAYES TO BE DONE! # 
extract_corr_bias <- function(df, type)
{
    df %>% spread(method,dnds) %>% 
           group_by(dataset) %>% 
           do(rraw = cor(.$nopenal, .$phylobayes), braw = lm(phylobayes ~ offset(nopenal), data=.)) %>% 
           mutate(r = rraw[1], b = summary(braw)$coeff[1]) %>% select(-rraw, -braw) %>% arrange(b,r) -> newdf
    newdf$type <- type
    newdf
}


dnds.sim <- read.csv("../../results/summarized_results/simulation_derived_dnds.csv")
dnds.sim %>% filter(del == "weak") %>% extract_corr_bias("weak") -> sim.weak.r.b
dnds.sim %>% filter(del == "strong") %>% extract_corr_bias("strong") -> sim.strong.r.b
dnds.emp <- read.csv("../../results/summarized_results/empirical_derived_dnds.csv")
emp.r.b <- extract_corr_bias(dnds.emp, "empirical")
full.r.b <- rbind(sim.weak.r.b, sim.strong.r.b, emp.r.b)
full.r.b$type <- factor(full.r.b$type, levels=c("strong", "weak", "empirical"), labels=c("Simulated, strong", "Simulated, weak", "Empirical"))
rplot_swpb <- ggplot(full.r.b, aes(x = type, y = r, group = type, fill=type)) + geom_boxplot() + scale_fill_manual(values=c("orangered2", "royalblue3")) + xlab("Data type") + ylab("swMutSel-pbMutSel Correlation")
bplot_swpb <- ggplot(full.r.b, aes(x = type, y = b, group = type, fill=type)) + geom_boxplot() + scale_fill_manual(values=c("orangered2", "royalblue3")) + xlab("Data type") + ylab("Average swMutSel > pbMutSel")
rb <- plot_grid(rplot, bplot, nrow=1)

