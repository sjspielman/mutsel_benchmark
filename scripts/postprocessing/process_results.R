# SJS
# Script to merge dN/dS calculations from FEL, swMutSel, and pbMutSel into a single data frame.

require(dplyr)
require(tidyr)

results_dir    <- "../../results/"
datadir_parent <- paste0(results_dir,"raw_results/")
sim.datasets   <- c("1B4T_A_simulated", "1G58_B_simulated", "1GV3_A_simulated", "1HUR_A_simulated", "1IBS_A_simulated", "1PV1_A_simulated", "1QMV_A_simulated", "1R6M_A_simulated", "1V9S_B_simulated", "1W7W_B_simulated", "1X1O_B_simulated", "1YPI_A_simulated", "1ZNN_A_simulated", "2A84_A_simulated", "2BCG_Y_simulated", "2CFE_A_simulated", "2CJM_C_simulated", "2CNV_A_simulated", "2FLI_A_simulated", "2G0N_B_simulated")
emp.datasets   <- c("PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")
datasets       <- c(sim.datasets, emp.datasets)
methods        <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000", "phylobayes")


results.dnds <- data.frame("dataset" = factor(), 
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

    # dN/dS from SLAC 
    dat <- read.table(paste0(datadir, "slac/", dataset, "_SLAC.txt"), header=T)
    datdnds <- dat$dN/mean(dat$dS)
    l <- nrow(dat)
    temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = "slac", type = type)
    results.dnds <- rbind(results.dnds, temp)

    # dN/dS predicted from MutSel inference
    for (meth in methods){    
        if (meth == "phylobayes"){
          savemeth <- "pbmutsel"
        }
        else{
          savemeth <- meth
        }
        dat <- read.table(paste0(datadir, "derived_dnds/", dataset, "_", meth, "_dnds.txt"))
        datdnds <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = dataset, "site" = 1:l, "dnds" = datdnds , "method" = savemeth, type = type)
        results.dnds <- rbind(results.dnds, temp)
    }
    
    # True simulated dN/dS AS APPLICABLE
    if (type == "simulation")
    {
        newname <- strsplit(dataset, "_simulated")[[1]][1]
        true <- read.csv(paste0(datadir, "derived_dnds/", newname, "_true_dnds.txt"))
        temp <- data.frame("dataset" = dataset, "site" = true$site, "dnds" = true$dnds , "method" = "true", type = type)
        results.dnds <- rbind(results.dnds, temp)    
    }    
}

write.csv(results.dnds, paste0(results_dir, "dnds_results.csv", sep=""), row.names=FALSE, quote=FALSE)





results.jsd <- data.frame("dataset" = factor(), 
                          "site"    = numeric(), 
                          "jsd"    = numeric(),
                          "method"  = factor())


for (dataset in sim.datasets){
    
    newname <- strsplit(dataset, "_simulated")[[1]][1]
    datadir <- paste0(datadir_parent,"simulation/jsd/")

    for (meth in methods){    
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
write.csv(results.jsd, paste0(results_dir, "jsd_results.csv", sep=""), row.names=FALSE, quote=FALSE)

