# SJS
# Statistical analysis

require(lme4)
require(multcomp)
require(dplyr)

result_directory <- "../../results/summarized_results/"
jsd.dat  <- read.csv(paste0(result_directory, "simulation_jsd.csv"))
sim.dnds <- read.csv(paste0(result_directory, "simulation_derived_dnds.csv"))


fit1 <- lmer(jsd ~ method + del + (1|dataset), data = jsd.dat)
fit1.mc.method <- glht(fit1, linfct=mcp(method='Tukey'))
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
# mvn100 - mvn10 == 0        0.0012580  0.0013270   0.948  0.98125         !!mvn100 = mvn10
# nopenal - mvn10 == 0       0.0013751  0.0013270   1.036  0.96895         !!nopenal = mvn10
# phylobayes - mvn10 == 0    0.0700083  0.0013748  50.923  < 0.001 ***
# nopenal - mvn100 == 0      0.0001171  0.0013270   0.088  1.00000         !!nopenal = mvn100
# phylobayes - mvn100 == 0   0.0687503  0.0013748  50.008  < 0.001 ***
# phylobayes - nopenal == 0  0.0686332  0.0013748  49.923  < 0.001 ***


fit1.mc.del <- glht(fit1, linfct=mcp(del='Tukey'))
summary(fit1.mc.del)
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)
# weak - strong == 0 0.0003191  0.0006696   0.477    0.634                 !!strong = weak


