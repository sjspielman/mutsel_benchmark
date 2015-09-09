# SJS
# Script to merge results into a single data frame.
require(dplyr)
require(tidyr)

datadir_parent="/Users/sjspielman/Research/mutsel_bench/results/"
meths   <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000")# ,"pb")   # files are a single column of values
datasets <- c("1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated")    #, "amine", "pb2", "PF00126", "PF04055", "PF00593", "PF07715") #, "PF01266", "PF02518", "camelid", "amine", "pb2", 
results <- data.frame("dataset" = factor(), 
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "type"    = factor(),
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
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = "fel1", type = type)
    results <- rbind(results, temp)


    # dN/dS predicted from MutSel inference
    for (meth in meths){    
        dat <- read.table(paste0(datadir, dataset, "_", meth, "_dnds.txt"))
        datdnds <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = meth, type = type)
        results <- rbind(results, temp)
    }
#     
#     
#     # dN/dS predicted directly from alignments
#     dat <- read.csv(paste0(datadir, dataset, "_dnds_from_alignment.txt"))
#     #dndscodon <- dat$dnds_codon
#     empirical    <- dat$dnds_aa
#     l <- nrow(dat)
#     #temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = dndscodon , "method" = "aln_codon")
#     #results <- rbind(results, temp)
#     temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = empirical , "method" = "empirical", type = type)
#     results <- rbind(results, temp)
#     
    
    # True simulated dN/dS AS APPLICABLE
    if (type == "simulation")
    {
        newname <- strsplit(dataset, "_simulated")[[1]][1]
        true <- read.csv(paste0(datadir, newname, "_true_dnds.txt"))
        temp <- data.frame("dataset" = dataset, "site" = true$site, "dnds" = true$dnds , "method" = "true", type = type)
        results <- rbind(results, temp)    
    }    
        
        
}


write.csv(results, paste0(datadir_parent, "results.csv", sep=""), row.names=FALSE, quote=FALSE)
