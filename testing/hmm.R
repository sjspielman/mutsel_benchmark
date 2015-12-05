#require(cowplot)
#require(dplyr)
names <- c("PF00188", "PF00226", "PF00512", "PF03144", "PF04542")
for (name in names){

empslacraw <- read.table(paste0(name,"_SLAC.txt"), header=T)
simslacraw <- read.table(paste0(name,"_simulated_SLAC.txt"), header=T)
empslac <- empslacraw$dN/mean(empslacraw$dS)
simslac <- simslacraw$dN/mean(simslacraw$dS)
simsw <- read.table(paste0(name,"_nopenal_dnds.txt"))
empsw <- read.table(paste0(name,"_simulated_dnds.txt"))

pdf(paste0(name, ".pdf"), width=9, height=3)
par(mfrow=c(1,4))
f <- lm(empsw$V1~empslac)
plot(empslac,empsw$V1, pch=20, xlab = "empirical slac", ylab="empirical swmutsel",main=f$coefficients[[2]])
abline(0,1)
abline(f, col="red")

f <- lm(simsw$V1~simslac)
plot(simslac,simsw$V1, pch=20, xlab = "simulated slac", ylab="simulated swmutsel",main=f$coefficients[[2]])
abline(0,1)
abline(f, col="red")

f <- lm(empslac~simslac)
plot(simslac,empslac, pch=20, xlab = "simulated slac", ylab="empirical slac",main=f$coefficients[[2]])
abline(0,1)
abline(f, col="red")

f <- lm(empsw$V1~simsw$V1)
plot(simsw$V1,empsw$V1, pch=20, xlab = "simulated swmutsel", ylab="empirical swmutsel",main=f$coefficients[[2]])
abline(0,1)
abline(f, col="red")

dev.off()
}
