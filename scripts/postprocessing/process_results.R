# SJS
# Script to merge dN/dS predictions from swMutSel and Phylobayes into a single data frame, one for simulated and empirical data.

require(dplyr)
require(tidyr)

results_dir    <- "../../results/"
datadir_parent <- paste0(results_dir,"raw_results/")
sim.datasets   <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "1GV3_A", "2BCG_Y", "1IBS_A", "2CFE_A", "1R6M_A", "2FLI_A")
deltypes       <- c("weak", "strong")
emp.datasets   <- c("PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")
methods        <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn1", "mvn10", "mvn100")#  , "phylobayes")


# Simulated datasets #
sim.dat <- data.frame("dataset" = character(), 
                      "del"     = character(),
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "method"  = factor())
datadir <- paste0(datadir_parent, "simulation/")
truedir <- "../simulation/true_simulation_parameters/"
for (name in sim.datasets)
{    
    for (deltype in deltypes)
    { 
        true <- read.csv(paste0(truedir, name, "_del", deltype, "_true_dnds.csv"))

        for (meth in c("true", methods))
        {  
        
            if (meth == "true")
            {
                dat <- read.csv(paste0(truedir, name, "_del", deltype, "_true_dnds.csv"))
                datdnds <- dat$dnds
            }
            else
            {
                dat <- read.table(paste0(datadir, "derived_dnds_coeffs/", name, "_del", deltype, "_",  meth, "_dnds.csv"))
                datdnds <- dat$V1
            }
            l <- nrow(dat)
            temp <- data.frame("dataset" = name, "del" = deltype, "site" = 1:l, "dnds" = datdnds , "method" = meth)
            sim.dat <- rbind(sim.dat, temp)
        }
    }
}
write.csv(sim.dat, paste0(results_dir, "simulation_derived_dnds.csv", sep=""), row.names=FALSE, quote=FALSE)



# Empirical datasets #
emp.dat <- data.frame("dataset" = character(), 
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "method"  = factor())
datadir <- paste0(datadir_parent, "empirical/")
for (name in emp.datasets)
{    
    for (meth in methods){    
        dat <- read.table(paste0(datadir, "derived_dnds_coeffs/", name, "_", meth, "_dnds.csv"))
        datdnds <- dat$V1
        l <- nrow(dat)
        temp <- data.frame("dataset" = name, "site" = 1:l, "dnds" = datdnds , "method" = meth)
        emp.dat <- rbind(emp.dat, temp)
    }
}
write.csv(emp.dat, paste0(results_dir, "empirical_derived_dnds.csv", sep=""), row.names=FALSE, quote=FALSE)
