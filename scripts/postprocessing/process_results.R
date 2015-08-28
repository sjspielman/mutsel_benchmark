# SJS
# Script to merge results into a single data frame.
require(dplyr)
require(tidyr)

datadir_parent="/Users/sjspielman/Research/mutsel_bench/results/"
meths   <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000")# ,"pb")   # files are a single column of values
datasets <- c("PF00126", "PF01336", "PF04055", "PF00593", "PF01926", "PF07715", "PF01266", "PF02518", "camelid", "amine", "pb2", "1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1RII_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated")

results <- data.frame("dataset" = factor(), 
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "method"  = factor())


for (dataset in datasets){
    
    if (grepl("_simulated",dataset)){
        type <- "simulation"
    }
    else{
        type <- "empirical"
    }
    datadir  <- paste0(datadir_parent,type,"/")

    # dN/dS from FEL1 
    dat <- read.csv(paste0(datadir, dataset, "_fel1.txt"))
    datdnds <- dat$dN.dS
    l <- nrow(dat)
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = "fel1")
    results <- rbind(results, temp)


    # dN/dS predicted from MutSel inference
    for (meth in meths){    
        dat <- read.table(paste0(datadir, dataset, "_", meth, "_dnds.txt"))
        datdnds <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = meth)
        results <- rbind(results, temp)
    }
    
    
    # dN/dS predicted directly from alignments
    dat <- read.csv(paste0(datadir, dataset, "_dnds_from_alignment.txt"))
    dndscodon <- dat$dnds_codon
    dndsaa    <- dat$dnds_aa
    l <- nrow(dat)
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = dndscodon , "method" = "aln_codon")
    results <- rbind(results, temp)
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = dndsaa , "method" = "aln_aa")
    results <- rbind(results, temp)
    
    
    # True simulated dN/dS AS APPLICABLE
    if (type == "simulation")
    {
        newname <- strsplit(dataset, "_simulated")[[1]][1]
        true <- read.csv(paste0(datadir, newname, "_true_dnds.txt"))
        temp <- data.frame("dataset" = dataset, "site" = true$site, "dnds" = true$dnds , "method" = "true")
        results <- rbind(results, temp)    
    }    
        
        
}


write.csv(results, paste0(datadir_parent, "results.csv", sep=""), row.names=FALSE, quote=FALSE)


#dat.raw %>% na.omit() %>% mutate(pos = (pvalue <= alpha & dnds > 1)) %>% group_by(true, ntaxa, bl) %>% mutate(percent_pos = sum(pos)/numcol, treelen = (2*(ntaxa-1)*bl)) -> dat
#write.csv(dat, paste(DATADIR, "homog_results.csv", sep=""), row.names=FALSE, quote=FALSE)



# 
# 
# 
# plot_dnds <- function(df, mu_type, reference_values, label)
# {
#     for (data in unique(df$dataset))
#     {
#         df %>% filter(dataset == data) -> df2
#         outfile = paste("plots/", data, "_dnds_", label, ".pdf", sep="")
#         pdf(outfile, width = 9, height=5)
# 
#         par(mfrow=c(2,4))
#         for (m in unique(dat$method))
#         {
#             df2 %>% filter(method == m, mutype == mu_type) -> df3
#             r <- cor(log(reference_values), log(df3$dnds))
#             plot(reference_values, df3$dnds, log='xy', cex=0.75, main=paste(m, round(r,4), sep="; r="), cex.main = 0.95, xlab = "Reference dN/dS", ylab = "Predicted from MutSel dN/dS", pch=20, xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
#             abline(0,1,col='red', lwd=2)           
#             
#         }   
#         dev.off()
#     }
# }
# 
# plot_dnds(dat, "fixed", alltrue, "true_fixedmu")
# plot_dnds(dat, "est", alltrue, "true_estmu")
# 
# plot_dnds(dat, "fixed", slac_dnds, "slac_fixedmu")
# plot_dnds(dat, "est", slac_dnds, "slac_estmu")
# 
# plot_dnds(dat, "fixed", pred_aa, "predaa_fixedmu")
# plot_dnds(dat, "est", pred_aa, "predaa_estmu")
# 
# r <- corr(log(alltrue), log(slac_dnds))
# pdf("plots/n11_bl0.16_slac_true.pdf", width =5, height=5)
# plot(alltrue, slac_dnds, xlab="True dN/dS", ylab = "SLAC dN/dS", main = round(r,4), cex = 0.75, pch=20, log="xy", xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
# abline(0,1,col="red")
# dev.off()
# 
# r <- corr(log(predaa), log(slac_dnds))
# pdf("plots/n11_bl0.16_slac_predaa.pdf", width =5, height=5)
# plot(slac_dnds, pred_aa, xlab="SLAC dN/dS", ylab = "Predicted from AA freqs dN/dS", main = round(r,4), cex = 0.75, pch=20, log="xy", xlim=c(5e-3, 1), ylim=c(5e-3, 1), frame.plot=F)
# abline(0,1,col="red")
# dev.off()
# 



