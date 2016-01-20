# SJS
# Statistical analyses, in no particular order..... :)

require(lme4)
require(lmerTest)
require(multcomp)
require(dplyr)
require(tidyr)
require(readr)

result_directory <- "../results/summarized_results/"
true_directory   <- "../simulation/true_simulation_parameters/"
jsd.dat  <- read.csv(paste0(result_directory, "jsd.csv"))
sim.dnds <- read.csv(paste0(result_directory, "dnds.csv"))
sim.datasets   <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A")
deltypes       <- c("strong", "weak")

# Proportion of highly changes across data 
prop <- data.frame("dataset" = character(), "del" = character(), "prop" = numeric())
for (d in sim.datasets){
    for (del in deltypes){
        
        dat <- read_csv(paste0(true_directory, d, "_del", del, "_true_selcoeffs.csv"))
        total <- nrow(dat)
        big <- sum(dat$binnedcoeff == 10 | dat$binnedcoeff == -10)
        propor <- big/total
        temp <- data.frame("dataset" = d, "del" = del, "prop" = propor)
        prop <- rbind(temp,prop)
    }
}
prop %>% group_by(del) %>% summarize(minprop = min(prop), maxprop = max(prop), meanprop = mean(prop), sdprop = sd(prop))
#      del   minprop   maxprop  meanprop     sdprop
#   (fctr)     (dbl)     (dbl)     (dbl)      (dbl)
# 1   weak 0.0000000 0.0000000 0.0000000 0.00000000
# 2 strong 0.3892646 0.4302033 0.4059725 0.01346862



# Correlation between weak and strong true dN/dS
#sim.dnds %>% filter(method == "true") %>% group_by(dataset) %>% spread(del,dnds) %>% ungroup() %>% group_by(dataset) %>% do(rraw = cor(.$strong, .$weak)) %>% mutate(r = rraw[[1]]) %>% ungroup() %>% summarize(meanr = mean(r), sdr = sd(r)) 
#      meanr        sdr
#      (dbl)      (dbl)
#1 0.9464474 0.01203134


# All methods display highly significant differences between JSD distributions between strong,weak
for (m in unique(jsd.dat$method)){
    jsd.dat %>% filter(method == m) -> jsd.temp
    fit <- lmer(jsd ~ del + (1|dataset), data = jsd.temp)
    print(m)
    print(summary(fit))
    readline()
}
#method    weak-strong    pvalue      
#nopenal	1.483e-02	<2e-16
#d0.01	-6.324e-03	2.82e-05
#d0.1	-4.546e-02	<2e-16
#d1.0	-5.876e-02	<2e-16
#mvn1	-5.506e-02	<2e-16
#mvn10	6.916e-03	3.89e-05
#mvn100	1.452e-02	4.44e-16
#phylobayes	-4.532e-02	<2e-16



fit1 <- lmer(jsd ~ method + del + (1|dataset), data = jsd.dat)
fit1.mc.method <- glht(fit1, linfct=mcp(method='Tukey')) # JSD among methods
summary(fit1.mc.method)
# Linear Hypotheses:
#                             Estimate Std. Error z value Pr(>|z|)    
# d0.1 - d0.01 == 0          0.0154740  0.0012432  12.447   <0.001 ***
# d1.0 - d0.01 == 0          0.0718530  0.0012432  57.796   <0.001 ***
# mvn1 - d0.01 == 0          0.0057558  0.0012432   4.630   <0.001 ***
# mvn10 - d0.01 == 0        -0.0103028  0.0012432  -8.287   <0.001 ***
# mvn100 - d0.01 == 0       -0.0069815  0.0012432  -5.616   <0.001 ***
# nopenal - d0.01 == 0      -0.0068210  0.0012432  -5.487   <0.001 ***
# phylobayes - d0.01 == 0    0.0333722  0.0012432  26.843   <0.001 ***
# d1.0 - d0.1 == 0           0.0563790  0.0012432  45.349   <0.001 ***
# mvn1 - d0.1 == 0          -0.0097182  0.0012432  -7.817   <0.001 ***
# mvn10 - d0.1 == 0         -0.0257768  0.0012432 -20.734   <0.001 ***
# mvn100 - d0.1 == 0        -0.0224556  0.0012432 -18.062   <0.001 ***
# nopenal - d0.1 == 0       -0.0222950  0.0012432 -17.933   <0.001 ***
# phylobayes - d0.1 == 0     0.0178982  0.0012432  14.397   <0.001 ***
# mvn1 - d1.0 == 0          -0.0660972  0.0012432 -53.166   <0.001 ***
# mvn10 - d1.0 == 0         -0.0821558  0.0012432 -66.083   <0.001 ***
# mvn100 - d1.0 == 0        -0.0788345  0.0012432 -63.412   <0.001 ***
# nopenal - d1.0 == 0       -0.0786740  0.0012432 -63.283   <0.001 ***
# phylobayes - d1.0 == 0    -0.0384808  0.0012432 -30.953   <0.001 ***
# mvn10 - mvn1 == 0         -0.0160587  0.0012432 -12.917   <0.001 ***
# mvn100 - mvn1 == 0        -0.0127374  0.0012432 -10.246   <0.001 ***
# nopenal - mvn1 == 0       -0.0125769  0.0012432 -10.116   <0.001 ***
# phylobayes - mvn1 == 0     0.0276164  0.0012432  22.214   <0.001 ***
# mvn100 - mvn10 == 0        0.0033213  0.0012432   2.672   0.1315       mvn100 = mvn10
# nopenal - mvn10 == 0       0.0034818  0.0012432   2.801   0.0949 .     nopenal = mvn10
# phylobayes - mvn10 == 0    0.0436750  0.0012432  35.131   <0.001 ***
# nopenal - mvn100 == 0      0.0001605  0.0012432   0.129   1.0000       nopenal = mvn100
# phylobayes - mvn100 == 0   0.0403537  0.0012432  32.459   <0.001 ***
# phylobayes - nopenal == 0  0.0401932  0.0012432  32.330   <0.001 ***

fit1.mc.del <- glht(fit1, linfct=mcp(del='Tukey'))   # Effect of weak/strong deleteriousness on JSD, controlling for method and dataset
summary(fit1.mc.del)
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)
# weak - strong == 0 -0.0218335  0.0006216  -35.12   <2e-16 ***  # weakly deleterious has slightly lower JSD then highly deleterious




sim.dnds %>% spread(method,dnds) %>% filter(del == "strong") %>% 
             gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
             select(dataset, del, site, true, dnds, method) %>% group_by(dataset, method) %>%
             do(rraw = cor(.$true, .$dnds)) %>% mutate(r2 = rraw[[1]]^2) %>% select(-rraw) -> strong.corrs

sim.dnds %>% spread(method,dnds) %>% filter(del == "weak") %>% 
             gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
             select(dataset, del, site, true, dnds, method) %>% group_by(dataset, method) %>%
             do(rraw = cor(.$true, .$dnds)) %>% mutate(r2 = rraw[[1]]^2) %>% select(-rraw) -> weak.corrs
             
sim.dnds %>% spread(method,dnds) %>%
             gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
             select(dataset, del, site, true, dnds, method) %>% group_by(dataset, del, method) %>%
             do(rraw = cor(.$true, .$dnds)) %>% mutate(r2 = rraw[[1]]^2) %>% select(-rraw) -> all.corrs

sim.dnds %>% spread(method,dnds) %>%
             gather(method, dnds, d0.01, d0.1, d1.0, mvn1, mvn10, mvn100, nopenal, phylobayes) %>% 
             select(dataset, del, site, true, dnds, method) %>% group_by(dataset, del, method) %>%
             do(braw = glm(dnds ~ offset(true), dat=.)) %>% mutate(b = summary(braw)$coeff[1]) %>% select(-braw) -> all.estbias

fit2 <- lmer(r2 ~ method  + (1|dataset), data = strong.corrs)
fit2.mc.method <- glht(fit2, linfct=mcp(method='Tukey')) # dN/dS r^2 among methods
summary(fit2.mc.method)
# Linear Hypotheses:
#                             Estimate Std. Error z value Pr(>|z|)    
# d0.1 - d0.01 == 0         -1.260e-01  1.471e-02  -8.561   <0.001 ***
# d1.0 - d0.01 == 0         -2.799e-01  1.471e-02 -19.024   <0.001 ***
# mvn1 - d0.01 == 0         -2.911e-02  1.471e-02  -1.979   0.4964     mvn1    = d0.01
# mvn10 - d0.01 == 0         8.989e-03  1.471e-02   0.611   0.9988     mvn10   = d0.01
# mvn100 - d0.01 == 0        8.063e-03  1.471e-02   0.548   0.9994     mvn100  = d0.01
# nopenal - d0.01 == 0       8.053e-03  1.471e-02   0.547   0.9994     nopenal = d0.01
# phylobayes - d0.01 == 0   -1.777e-01  1.471e-02 -12.078   <0.001 ***
# d1.0 - d0.1 == 0          -1.539e-01  1.471e-02 -10.463   <0.001 ***
# mvn1 - d0.1 == 0           9.684e-02  1.471e-02   6.582   <0.001 ***
# mvn10 - d0.1 == 0          1.349e-01  1.471e-02   9.172   <0.001 ***
# mvn100 - d0.1 == 0         1.340e-01  1.471e-02   9.109   <0.001 ***
# nopenal - d0.1 == 0        1.340e-01  1.471e-02   9.108   <0.001 ***
# phylobayes - d0.1 == 0    -5.175e-02  1.471e-02  -3.517   0.0101 *  
# mvn1 - d1.0 == 0           2.508e-01  1.471e-02  17.045   <0.001 ***
# mvn10 - d1.0 == 0          2.889e-01  1.471e-02  19.635   <0.001 ***
# mvn100 - d1.0 == 0         2.880e-01  1.471e-02  19.572   <0.001 ***
# nopenal - d1.0 == 0        2.880e-01  1.471e-02  19.571   <0.001 ***
# phylobayes - d1.0 == 0     1.022e-01  1.471e-02   6.946   <0.001 ***
# mvn10 - mvn1 == 0          3.810e-02  1.471e-02   2.590   0.1594      mvn10  = mvn1
# mvn100 - mvn1 == 0         3.717e-02  1.471e-02   2.527   0.1839      mvn100 = mvn1
# nopenal - mvn1 == 0        3.717e-02  1.471e-02   2.526   0.1846      nopenal = mvn1
# phylobayes - mvn1 == 0    -1.486e-01  1.471e-02 -10.099   <0.001 ***
# mvn100 - mvn10 == 0       -9.268e-04  1.471e-02  -0.063   1.0000      mvn100 = mvn10
# nopenal - mvn10 == 0      -9.366e-04  1.471e-02  -0.064   1.0000      nopenal = mvn10
# phylobayes - mvn10 == 0   -1.867e-01  1.471e-02 -12.689   <0.001 ***
# nopenal - mvn100 == 0     -9.804e-06  1.471e-02  -0.001   1.0000      nopenal = mvn100
# phylobayes - mvn100 == 0  -1.858e-01  1.471e-02 -12.626   <0.001 ***
# phylobayes - nopenal == 0 -1.858e-01  1.471e-02 -12.625   <0.001 ***


fit2 <- lmer(r2 ~ method  + del +  (1|dataset), data = all.corrs)
fit2.mc.method <- glht(fit2, linfct=mcp(del='Tukey')) # dnds correlations between strong, weak deleterious (controlling for method and dataset) are the same.
summary(fit2.mc.method)
#Linear Hypotheses:
#                    Estimate Std. Error z value Pr(>|z|)
#weak - strong == 0 -0.003217   0.006420  -0.501    0.616


fit2 <- lmer(b ~ method * del +  (1|dataset), data = all.estbias)
summary(fit2)
#                           Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)              -4.479e-03  2.814e-03  3.071e+01  -1.592  0.12167    
# methodd0.1                6.161e-02  2.728e-03  1.501e+02  22.580  < 2e-16 ***
# methodd1.0                2.638e-01  2.728e-03  1.501e+02  96.700  < 2e-16 ***
# methodmvn1                1.408e-01  2.728e-03  1.501e+02  51.613  < 2e-16 ***
# methodmvn10              -7.543e-03  2.728e-03  1.501e+02  -2.765  0.00641 ** 
# methodmvn100             -1.247e-02  2.728e-03  1.501e+02  -4.569 1.02e-05 ***
# methodnopenal            -1.266e-02  2.728e-03  1.501e+02  -4.640 7.52e-06 ***
# methodphylobayes          9.050e-02  2.728e-03  1.501e+02  33.168  < 2e-16 ***
# delweak                  -2.159e-02  2.728e-03  1.501e+02  -7.912 5.13e-13 ***
# methodd0.1:delweak       -1.825e-02  3.858e-03  1.501e+02  -4.729 5.15e-06 ***
# methodd1.0:delweak       -2.264e-02  3.858e-03  1.501e+02  -5.867 2.73e-08 ***
# methodmvn1:delweak        3.408e-04  3.858e-03  1.501e+02   0.088  0.92975    
# methodmvn10:delweak       7.291e-03  3.858e-03  1.501e+02   1.890  0.06075 .  
# methodmvn100:delweak      6.508e-03  3.858e-03  1.501e+02   1.687  0.09375 .  
# methodnopenal:delweak     6.474e-03  3.858e-03  1.501e+02   1.678  0.09545 .  
# methodphylobayes:delweak  1.896e-02  3.858e-03  1.501e+02   4.913 2.32e-06 ***

pvalues <- c()
diffmean <- c()
for (m in c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")){
    all.estbias %>% filter(method == m) %>% spread(del, b) -> temp.estbias
    p <- t.test(abs(temp.estbias$strong), abs(temp.estbias$weak), paired=T)$p.value * 8
    dmean <- t.test(abs(temp.estbias$strong), abs(temp.estbias$weak), paired=T)$estimate[[1]]
    pvalues <- c(pvalues,p)
    diffmean <- c(dmean, diffmean)
}
# all.estbias %>% group_by(method,del) %>% summarize(mean(b))
#        method    del      mean(b)
#        (fctr) (fctr)        (dbl)
# 1       d0.01 strong -0.004478664
# 2       d0.01   weak -0.026066687
# 3        d0.1 strong  0.057129532
# 4        d0.1   weak  0.017293899
# 5        d1.0 strong  0.259357311
# 6        d1.0   weak  0.215130947
# 7        mvn1 strong  0.136343745
# 8        mvn1   weak  0.115096483
# 9       mvn10 strong -0.012022196
# 10      mvn10   weak -0.026319318
# 11     mvn100 strong -0.016943627
# 12     mvn100   weak -0.032023686
# 13    nopenal strong -0.017139386
# 14    nopenal   weak -0.032253245
# 15 phylobayes strong  0.086017574
# 16 phylobayes   weak  0.083388372





##########################################################################################
################### K-S tests between selection coefficients #############################
datasets <- unique(sim.dnds$dataset)
methods <- c("nopenal", "mvn100", "mvn10", "mvn1", "d0.01", "d0.1", "d1.0", "phylobayes")
n <- length(datasets)*length(methods)


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
ksp.st %>% filter(pvalue >= 0.01)
# [1] dataset method  pvalue 
# <0 rows> (or 0-length row.names)


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
ksp.wt %>% filter(pvalue >= 0.01)
# [1] dataset method  pvalue 
# <0 rows> (or 0-length row.names)

