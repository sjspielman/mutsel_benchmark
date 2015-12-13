require(cowplot)
# path <- "../results/raw_results/simulation/derived_dnds_coeffs/"
# names   <- c("1B4T_A", "1G58_B", "1GV3_A", "1HUR_A", "1IBS_A", "1PV1_A", "1QMV_A", "1R6M_A", "1V9S_B", "1W7W_B", "1X1O_B", "1YPI_A", "1ZNN_A", "2A84_A", "2BCG_Y", "2CFE_A", "2CJM_C", "2CNV_A", "2FLI_A", "2G0N_B")
# for (name in names){
# true <- read.csv(paste0(path, name, "_true_selcoeffs.csv"))
# nopenal <- read.csv(paste0(path, name, "_simulated_nopenal_selcoeffs.csv"))
# pbmutsel <- read.csv(paste0(path, name, "_simulated_phylobayes_selcoeffs.csv"))
# mvn10 <- read.csv(paste0(path, name, "_simulated_mvn10_selcoeffs.csv"))
# mvn100 <- read.csv(paste0(path, name, "_simulated_mvn100_selcoeffs.csv"))
# mvn1000 <- read.csv(paste0(path, name, "_simulated_mvn1000_selcoeffs.csv"))
# d0.01 <- read.csv(paste0(path, name, "_simulated_d0.01_selcoeffs.csv"))
# d0.1 <- read.csv(paste0(path, name, "_simulated_d0.1_selcoeffs.csv"))
# d1.0 <- read.csv(paste0(path, name, "_simulated_d1.0_selcoeffs.csv"))
# 
# bw=0.75
# a <- ggplot(true, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("true")
# b <- ggplot(nopenal, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("nopenal")
# cc <- ggplot(pbmutsel, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("pbmutsel")
# d <- ggplot(mvn10, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("mvn10")
# e <- ggplot(mvn100, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("mvn100")
# f <- ggplot(mvn1000, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("mvn1000")
# g <- ggplot(d0.01, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("d0.01")
# h <- ggplot(d0.1, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("d0.1")
# i <- ggplot(d1.0, aes(x = V1)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-10,10)) + ggtitle("d1.0")
# p<-plot_grid(a,b,cc,d,e,f,g,h,i,nrow=3)
# print(p)
# readline()}


path <- "../results/raw_results/simulated_from_empirical/derived_dnds_coeffs/"
path2 <- "../results/raw_results/empirical/derived_dnds_coeffs/"
names  <- c("PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")

for (name in names){
true <- read.csv(paste0(path2, name, "_nopenal_selcoeffs.csv"))
nopenal <- read.csv(paste0(path, name, "_simulated_nopenal_selcoeffs.csv"))
d1.0 <- read.csv(paste0(path, name, "_simulated_d1.0_selcoeffs.csv"))

bw=0.75
a <- ggplot(true, aes(x = realcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("true, raw")
b <- ggplot(true, aes(x = binnedcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("true, binned")
cc <- ggplot(nopenal, aes(x = realcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("nopenal, raw")
d <- ggplot(nopenal, aes(x = binnedcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("nopenal, binned")
e <- ggplot(d1.0, aes(x = realcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("d1.0, raw")
f <- ggplot(d1.0, aes(x = binnedcoeff)) + geom_histogram(binwidth = bw) + scale_x_continuous(limits=c(-11,11)) + ggtitle("d1.0, binned")

p<-plot_grid(a,b,cc,d,e,f,nrow=3)
print(p)
readline()}