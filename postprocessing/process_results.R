# SJS
# Script to merge JSD, dN/dS, and selection coefficient inferences from swMutSel and Phylobayes into concise single data frame(s).

require(dplyr)
require(tidyr)

resdir       <- "../results/"
truedir      <- "../simulation/true_simulation_parameters/"
datadir      <- paste0(resdir, "dnds_coeffs_jsd/")
outdir       <- paste0(resdir, "summarized_results/")

sim.datasets <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A")
deltypes     <- c("weak", "strong")
methods      <- c("nopenal", "d0.01", "d0.1", "d1.0", "mvn1", "mvn10", "mvn100" , "phylobayes")


print("Processing JSD")
jsd <- data.frame("dataset" = character(), 
                  "del"     = character(),
                  "site"    = numeric(), 
                  "jsd"     = numeric(),
                  "method"  = factor())
for (name in sim.datasets)
{    
    for (deltype in deltypes)
    { 
        for (meth in methods)
        {  
            dat <- read.csv(paste0(datadir, name, "_del", deltype, "_",  meth, "_jsd.csv"))
            l <- nrow(dat)
            temp <- data.frame("dataset" = name, "del" = deltype, "site" = dat$site, "jsd" = dat$jsd, "method" = meth)
            jsd <- rbind(jsd, temp)
        }
    }
}
write.csv(jsd, paste0(outdir, "simulation_jsd.csv"), row.names=FALSE, quote=FALSE)






print("Processing dN/dS")
sim.dat <- data.frame("dataset" = character(), 
                      "del"     = character(),
                      "site"    = numeric(), 
                      "dnds"    = numeric(),
                      "method"  = factor())
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
            }
            else
            {
                dat <- read.csv(paste0(datadir, name, "_del", deltype, "_",  meth, "_dnds.csv"))
            }
            temp <- data.frame("dataset" = name, "del" = deltype, "site" = dat$site, "dnds" = dat$dnds , "method" = meth)
            sim.dat <- rbind(sim.dat, temp)
        }
    }
}
write.csv(sim.dat, paste0(outdir, "simulation_derived_dnds.csv"), row.names=FALSE, quote=FALSE)



print("Processing selection coefficients")
for (dataset in sim.datasets)
{
    for (del in deltypes){
        sim.selcoeffs <- data.frame("dataset" = character(), "method" = character(), "dummy" = numeric(), "del" = character(), "binnedcoeff" = numeric(), "realcoeff" = numeric())
        full_dataset <- paste0(dataset, "_del", del)
        for (m in c("true", methods))
        {
            if (m == "true")
            {
                directory <- truedir
            }
            else
            {
                directory <- datadir
            }
            
            dat <- read.csv(paste0(directory, full_dataset, "_", m, "_selcoeffs.csv"))
            temp <- data.frame("dataset" = dataset, "method" = m, "dummy" = 1:nrow(dat), "del" = del, "binnedcoeff" = dat$binnedcoeff, "realcoeff" = dat$realcoeff)
            sim.selcoeffs <- rbind(sim.selcoeffs, temp)
        }
        write.csv(sim.selcoeffs, paste0(outdir, full_dataset, "_selection_coefficients.csv"), row.names=FALSE, quote=FALSE)
    }
}
