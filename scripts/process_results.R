# SJS
# Script to merge results into a single data frame, sans phylobayes. Specifically for simulated data.
require(dplyr)

DATADIR="/Users/sjspielman/Research/mutsel_bench/results/"
meths    <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000")
mutypes  <- c("est", "fixed")
datasets <- c("n11_bl0.16")
numcol = 185
dat <- data.frame("dataset" = factor(), 
                  "site"    = numeric(), 
                  "true"    = numeric(),
                  "slac"    = numeric(),
                  "dnds"    = numeric(),
                  "method"  = numeric(),
                  "mutype"  = factor(),
                  "predc"   = numeric(),
                  "predaa"  = numeric())

dftrue <- read.table(paste(DATADIR,"codon_freq_lib_info.txt", sep=""), header=T)
alltrue <- dftrue$dnds
for (dataset in datasets){
    slac <- read.table(paste(DATADIR, dataset, "_slac.txt", sep=""), header=T)
    slac_dnds <- slac$dN/slac$dS
    
    pred <- read.csv(paste(DATADIR, dataset, "_predicted_dnds.csv", sep="")) 
    pred_codon <- pred$dnds_codon
    pred_aa    <- pred$dnds_aa
    
    for (meth in meths){        
        for (mu in mutypes){
            d <- read.table(paste(DATADIR, dataset, "_swmutsel_", meth, "_", mu, "mu_dnds.txt", sep=""))
            dnds <- d$V1

            temp <- data.frame("dataset" = rep(dataset, numcol), "site" = 1:numcol, "true" = alltrue, "slac" = slac_dnds, "dnds" = dnds , "method" = rep(meth, numcol), "mutype" = rep(mu, numcol), "pred_codon" = pred_codon, "pred_aa" = pred_aa)
            dat <- rbind(dat, temp)
        
        }
    }
}
dat %>% spread(method, dnds) %>% na.omit() -> dat

#dat.raw %>% na.omit() %>% mutate(pos = (pvalue <= alpha & dnds > 1)) %>% group_by(true, ntaxa, bl) %>% mutate(percent_pos = sum(pos)/numcol, treelen = (2*(ntaxa-1)*bl)) -> dat
#write.csv(dat, paste(DATADIR, "homog_results.csv", sep=""), row.names=FALSE, quote=FALSE)






plot_dnds <- function(df, mu_type, reference_values, label)
{
    for (data in unique(df$dataset))
    {
        df %>% filter(dataset == data) -> df2
        outfile = paste("plots/", data, "_dnds_", label, ".pdf", sep="")
        pdf(outfile, width = 9, height=5)

        par(mfrow=c(2,4))
        for (m in unique(dat$method))
        {
            df2 %>% filter(method == m, mutype == mu_type) -> df3
            r <- cor(log(reference_values), log(df3$dnds))
            plot(reference_values, df3$dnds, log='xy', cex=0.75, main=paste(m, round(r,4), sep="; r="), cex.main = 0.95, xlab = "Reference dN/dS", ylab = "Predicted from MutSel dN/dS", pch=20, xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
            abline(0,1,col='red', lwd=2)           
            
        }   
        dev.off()
    }
}

plot_dnds(dat, "fixed", alltrue, "true_fixedmu")
plot_dnds(dat, "est", alltrue, "true_estmu")

plot_dnds(dat, "fixed", slac_dnds, "slac_fixedmu")
plot_dnds(dat, "est", slac_dnds, "slac_estmu")

plot_dnds(dat, "fixed", pred_aa, "predaa_fixedmu")
plot_dnds(dat, "est", pred_aa, "predaa_estmu")

r <- corr(log(alltrue), log(slac_dnds))
pdf("plots/n11_bl0.16_slac_true.pdf", width =5, height=5)
plot(alltrue, slac_dnds, xlab="True dN/dS", ylab = "SLAC dN/dS", main = round(r,4), cex = 0.75, pch=20, log="xy", xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
abline(0,1,col="red")
dev.off()

r <- corr(log(predaa), log(slac_dnds))
pdf("plots/n11_bl0.16_slac_predaa.pdf", width =5, height=5)
plot(slac_dnds, pred_aa, xlab="SLAC dN/dS", ylab = "Predicted from AA freqs dN/dS", main = round(r,4), cex = 0.75, pch=20, log="xy", xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
abline(0,1,col="red")
dev.off()




