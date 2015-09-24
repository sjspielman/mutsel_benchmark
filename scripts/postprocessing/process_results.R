# SJS
# Script to merge dN/dS calculations from FEL, swMutSel, and pbMutSel into a single data frame.

require(dplyr)
require(tidyr)

datadir_parent="../../results/"
sim.datasets <- c("1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated")
emp.datasets <- c("amine", "pb2", "PF00126", "PF00593", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")
sim.meths    <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000", "phylobayes")
emp.meths    <- c("nopenal", "phylobayes")



results.dnds <- data.frame("dataset" = factor(), 
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "type"    = factor(),
                      "method"  = factor())


datasets <- c(sim.datasets, emp.datasets)
for (dataset in datasets){
    
    if (grepl("_simulated",dataset)){
        type <- "simulation"
        meths <- sim.meths
        
    }
    else{
        type <- "empirical"
        meths <- emp.meths
    }
    datadir  <- paste0(datadir_parent,type,"/")

    # dN/dS from FEL1 
    dat <- read.csv(paste0(datadir, dataset, "_FEL1.txt"))
    datdnds <- dat$dN.dS
    l <- nrow(dat)
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = "fel1", type = type)
    results.dnds <- rbind(results.dnds, temp)


    # dN/dS predicted from MutSel inference
    for (meth in meths){    
        if (meth == "phylobayes"){
          savemeth <- "pbmutsel"
        }
        else{
          savemeth <- meth
        }
        dat <- read.table(paste0(datadir, dataset, "_", meth, "_dnds.txt"))
        datdnds <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = savemeth, type = type)
        results.dnds <- rbind(results.dnds, temp)
    }
    
    # True simulated dN/dS AS APPLICABLE
    if (type == "simulation")
    {
        newname <- strsplit(dataset, "_simulated")[[1]][1]
        true <- read.csv(paste0(datadir, newname, "_true_dnds.txt"))
        temp <- data.frame("dataset" = dataset, "site" = true$site, "dnds" = true$dnds , "method" = "true", type = type)
        results.dnds <- rbind(results.dnds, temp)    
    }    
        
        
}
write.csv(results.dnds, paste0(datadir_parent, "dnds_results.csv", sep=""), row.names=FALSE, quote=FALSE)













results.jsd <- data.frame("dataset" = factor(), 
                          "site"    = numeric(), 
                          "jsd"    = numeric(),
                          "method"  = factor())


for (dataset in sim.datasets){
    
    newname <- strsplit(dataset, "_simulated")[[1]][1]
    datadir <- paste0(datadir_parent,"simulation/")

    for (meth in sim.meths){    
        if (meth == "phylobayes"){
          savemeth <- "pbmutsel"
        }
        else{
          savemeth <- meth
        }
        dat <- read.table(paste0(datadir, newname, "_", meth, "_jsd.txt"))
        jsd <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = dataset, "site" = 1:l, "jsd" = jsd , "method" = savemeth)
        results.jsd <- rbind(results.jsd, temp)
    }
       
        
}
write.csv(results.jsd, paste0(datadir_parent, "jsd_results.csv", sep=""), row.names=FALSE, quote=FALSE)

