# SJS
# Statistical analysis

require(lme4)
require(multcomp)
require(dplyr)
require(tidyr)
require(readr)

result_directory <- "../results/summarized_results/"
true_directory   <- "../scripts/simulation/true_simulation_parameters/"
jsd.dat  <- read.csv(paste0(result_directory, "simulation_jsd.csv"))
sim.dnds <- read.csv(paste0(result_directory, "simulation_derived_dnds.csv"))



# Correlation between weak and strong true dN/dS
sim.dnds %>% filter(method == "true") %>% group_by(dataset) %>% spread(del,dnds) %>% ungroup() %>% group_by(dataset) %>% do(rraw = cor(.$strong, .$weak)) %>% mutate(r = rraw[[1]]) %>% ungroup() %>% summarize(meanr = mean(r))
# 0.999979

fit1 <- lmer(jsd ~ method + del + (1|dataset), data = jsd.dat)
fit1.mc.method <- glht(fit1, linfct=mcp(method='Tukey')) # JSD among methods
summary(fit1.mc.method)
# Linear Hypotheses:
#                             Estimate Std. Error z value Pr(>|z|)    
# d0.1 - d0.01 == 0          0.0336028  0.0013270  25.322  < 0.001 ***
# d1.0 - d0.01 == 0          0.0966083  0.0013270  72.801  < 0.001 ***
# mvn1 - d0.01 == 0          0.0281787  0.0013270  21.235  < 0.001 ***
# mvn10 - d0.01 == 0        -0.0161990  0.0013270 -12.207  < 0.001 ***
# mvn100 - d0.01 == 0       -0.0149410  0.0013270 -11.259  < 0.001 ***
# nopenal - d0.01 == 0      -0.0148239  0.0013270 -11.171  < 0.001 ***
# phylobayes - d0.01 == 0    0.0538093  0.0013748  39.140  < 0.001 ***
# d1.0 - d0.1 == 0           0.0630056  0.0013270  47.479  < 0.001 ***
# mvn1 - d0.1 == 0          -0.0054241  0.0013270  -4.087  0.00116 **    
# mvn10 - d0.1 == 0         -0.0498018  0.0013270 -37.529  < 0.001 ***
# mvn100 - d0.1 == 0        -0.0485438  0.0013270 -36.581  < 0.001 ***
# nopenal - d0.1 == 0       -0.0484267  0.0013270 -36.493  < 0.001 ***
# phylobayes - d0.1 == 0     0.0202065  0.0013748  14.698  < 0.001 ***
# mvn1 - d1.0 == 0          -0.0684296  0.0013270 -51.567  < 0.001 ***
# mvn10 - d1.0 == 0         -0.1128073  0.0013270 -85.009  < 0.001 ***
# mvn100 - d1.0 == 0        -0.1115493  0.0013270 -84.061  < 0.001 ***
# nopenal - d1.0 == 0       -0.1114322  0.0013270 -83.972  < 0.001 ***
# phylobayes - d1.0 == 0    -0.0427990  0.0013748 -31.131  < 0.001 ***
# mvn10 - mvn1 == 0         -0.0443777  0.0013270 -33.442  < 0.001 ***
# mvn100 - mvn1 == 0        -0.0431197  0.0013270 -32.494  < 0.001 ***
# nopenal - mvn1 == 0       -0.0430026  0.0013270 -32.406  < 0.001 ***
# phylobayes - mvn1 == 0     0.0256306  0.0013748  18.643  < 0.001 ***
# mvn100 - mvn10 == 0        0.0012580  0.0013270   0.948  0.98125         !! mvn100 = mvn10
# nopenal - mvn10 == 0       0.0013751  0.0013270   1.036  0.96895         !! nopenal = mvn10
# phylobayes - mvn10 == 0    0.0700083  0.0013748  50.923  < 0.001 ***
# nopenal - mvn100 == 0      0.0001171  0.0013270   0.088  1.00000         !! nopenal = mvn100
# phylobayes - mvn100 == 0   0.0687503  0.0013748  50.008  < 0.001 ***
# phylobayes - nopenal == 0  0.0686332  0.0013748  49.923  < 0.001 ***


fit1.mc.del <- glht(fit1, linfct=mcp(del='Tukey'))   # Effect of weak/strong deleteriousness on JSD, controlling for method and dataset
summary(fit1.mc.del)
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)
# weak - strong == 0 0.0003191  0.0006696   0.477    0.634                 !! strong = weak




sim.dnds %>% spread(method,dnds) %>% 
             gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
             select(dataset, del, site, true, dnds, method) %>% group_by(dataset, method, del) %>%
             do(rraw = cor(.$true, .$dnds)) %>% mutate(r2 = rraw[[1]]^2) %>% select(-rraw) -> corrs

fit2 <- lmer(r2 ~ method + del + (1|dataset), data = corrs)
fit2.mc.method <- glht(fit2, linfct=mcp(method='Tukey')) # dN/dS r^2 among methods
summary(fit2.mc.method)
# Linear Hypotheses:
#                             Estimate Std. Error z value Pr(>|z|)    
# d0.1 - d0.01 == 0         -1.283e-01  1.022e-02 -12.552  < 0.001 ***
# d1.0 - d0.01 == 0         -2.812e-01  1.022e-02 -27.518  < 0.001 ***
# mvn1 - d0.01 == 0         -3.097e-02  1.022e-02  -3.030  0.05039 .      !! This is cusp.
# mvn10 - d0.01 == 0         9.046e-03  1.022e-02   0.885  0.98745        !! mvn10 = d0.01
# mvn100 - d0.01 == 0        8.444e-03  1.022e-02   0.826  0.99168        !! mvn100 = d0.01
# nopenal - d0.01 == 0       8.428e-03  1.022e-02   0.825  0.99179        !! noepenal = d0.01
# phylobayes - d0.01 == 0   -1.847e-01  1.022e-02 -18.074  < 0.001 ***
# d1.0 - d0.1 == 0          -1.529e-01  1.022e-02 -14.965  < 0.001 ***
# mvn1 - d0.1 == 0           9.732e-02  1.022e-02   9.523  < 0.001 ***
# mvn10 - d0.1 == 0          1.373e-01  1.022e-02  13.438  < 0.001 ***
# mvn100 - d0.1 == 0         1.367e-01  1.022e-02  13.379  < 0.001 ***
# nopenal - d0.1 == 0        1.367e-01  1.022e-02  13.377  < 0.001 ***
# phylobayes - d0.1 == 0    -5.644e-02  1.022e-02  -5.522  < 0.001 ***
# mvn1 - d1.0 == 0           2.503e-01  1.022e-02  24.488  < 0.001 ***
# mvn10 - d1.0 == 0          2.903e-01  1.022e-02  28.403  < 0.001 ***
# mvn100 - d1.0 == 0         2.897e-01  1.022e-02  28.344  < 0.001 ***
# nopenal - d1.0 == 0        2.897e-01  1.022e-02  28.342  < 0.001 ***
# phylobayes - d1.0 == 0     9.651e-02  1.022e-02   9.443  < 0.001 ***
# mvn10 - mvn1 == 0          4.001e-02  1.022e-02   3.915  0.00233 ** 
# mvn100 - mvn1 == 0         3.941e-02  1.022e-02   3.856  0.00295 ** 
# nopenal - mvn1 == 0        3.939e-02  1.022e-02   3.855  0.00299 ** 
# phylobayes - mvn1 == 0    -1.538e-01  1.022e-02 -15.044  < 0.001 ***
# mvn100 - mvn10 == 0       -6.014e-04  1.022e-02  -0.059  1.00000       !! mvn100 == mvn10
# nopenal - mvn10 == 0      -6.181e-04  1.022e-02  -0.060  1.00000       !! nopenal == mvn10
# phylobayes - mvn10 == 0   -1.938e-01  1.022e-02 -18.959  < 0.001 ***
# nopenal - mvn100 == 0     -1.674e-05  1.022e-02  -0.002  1.00000       !! nopenal == mvn100
# phylobayes - mvn100 == 0  -1.932e-01  1.022e-02 -18.901  < 0.001 ***
# phylobayes - nopenal == 0 -1.932e-01  1.022e-02 -18.899  < 0.001 ***


fit2.mc.del <- glht(fit2, linfct=mcp(del='Tukey')) # Effect of weak/strong deleteriousness on dN/dS correlations, controlling for method and dataset
summary(fit2.mc.del)
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)  
# weak - strong == 0 -0.008472   0.005110  -1.658   0.0973 .             !! strong = weak




##########################################################################################
################### K-S tests between selection coefficients #############################
datasets <- unique(sim.dnds$dataset)
methods <- c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")
n <- length(datasets)*length(methods)

##### K-S Tests between strong, weak. ######
ksp.sw <- data.frame(dataset = character(), method = character(), pvalue = numeric())
for (d in datasets){
    print(d)
    strong.sc <- read_csv(paste0(result_directory, d, "_delstrong_selection_coefficients.csv"))
    weak.sc <- read_csv(paste0(result_directory, d, "_delweak_selection_coefficients.csv"))
    
    for (m in methods){
        print(m)
        strong <- strong.sc %>% filter(method == m)
        weak <- weak.sc %>% filter(method == m)
        ks <- ks.test(strong$realcoeff, weak$realcoeff)
        corrected.p <- ks$p.value * n
        if (corrected.p >1) corrected.p = 1
        temp <- data.frame(dataset = d, method = m, pvalue = corrected.p)
        ksp.sw <- rbind(ksp.sw, temp)
}}
ksp.sw %>% filter(pvalue <= 0.01)



##### K-S Tests between strong and true. All P<2.2e-16 ######
ksp.st <- data.frame(dataset = character(), method = character(), pvalue = numeric())
for (d in datasets){
    print(d)
    strong.sc <- read_csv(paste0(result_directory, d, "_delstrong_selection_coefficients.csv"))    
    true <- strong.sc %>% filter(method == "true")
    
    for (m in methods){
        print(m)
        strong <- strong.sc %>% filter(method == m)    
        ks <- ks.test(strong$binnedcoeff, true$binnedcoeff)
        corrected.p <- ks$p.value * n
        if (corrected.p >1) corrected.p = 1
        temp <- data.frame(dataset = d, method = m, pvalue = corrected.p)
        ksp.st <- rbind(ksp.st, temp)
}}




##### K-S Tests between weak and true. All P<2.2e-16 ######
ksp.wt <- data.frame(dataset = character(), method = character(), pvalue = numeric())
for (d in datasets){
    print(d)
    weak.sc <- read_csv(paste0(result_directory, d, "_delweak_selection_coefficients.csv"))    
    true <- weak.sc %>% filter(method == "true")
    
    for (m in methods){
        print(m)
        weak <- weak.sc %>% filter(method == m)
        ks <- ks.test(weak$binnedcoeff, true$binnedcoeff)
        corrected.p <- ks$p.value * n
        if (corrected.p >1) corrected.p = 1
        temp <- data.frame(dataset = d, method = m, pvalue = corrected.p)
        ksp.wt <- rbind(ksp.wt, temp)
}}


