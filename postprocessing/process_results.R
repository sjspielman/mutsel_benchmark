# SJS
# Script to merge JSD, dN/dS, and selection coefficient inferences from swMutSel and Phylobayes into concise single data frame(s).

require(dplyr)
require(tidyr)

resdir       <- "dataframes/"
truedir      <- "../simulation/true_simulation_parameters/"
datadir      <- "dataframes/"

yeast.datasets <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A")
branch.lengths <- c("0.01", "0.5")
dms.datasets   <- c("HA", "NP")
deltypes       <- c("weak", "strong")
methods        <- c("nopenal", "d0.01", "d0.1", "mvn10", "mvn100" , "phylobayes") #"d1.0", "mvn1",

print("Processing yeast datasets")
dat.yeast <- data.frame("dataset"      = character(),
                        "del"          = character(),
                        "bl"           = numeric(),
                        "site"         = numeric(),
                        "jsd"          = numeric(),
                        "entropy"      = numeric(),
                        "diffsum"      = numeric(),
                        "dnds"         = numeric(),
                        "method"       = factor(), 
                        "true.dnds"    = numeric(),
                        "true.entropy" = numeric())
for (name in yeast.datasets)
{
    print(name)
    for (deltype in deltypes)
    {
        true <- read.csv(paste0(truedir, name, "_del", deltype, "_true_dnds_entropy.csv"))
        for (bl in branch.lengths)
        {
            for (meth in methods)
            {
                dat <- read.csv(paste0(datadir, name, "_del", deltype, "_bl", bl, "_", meth, "_statistics.csv"))
                l <- nrow(dat)
                temp <- data.frame("dataset" = name, "del" = deltype, "bl" = bl, "site" = dat$site, "jsd" = dat$jsd, "entropy" = dat$entropy, "diffsum" = dat$abs_sum_differences, "dnds" = dat$dnds, "method" = meth, "true.dnds" = true$dnds , "true.entropy" = true$entropy)
                dat.yeast <- rbind(dat.yeast, temp)
            }
        }
    }
}
write.csv(dat.yeast, paste0(datadir, "yeast_results.csv"), row.names=FALSE, quote=FALSE)


print("Processing deep mutational scanning")
dat.dms <- data.frame("dataset" = character(), 
                      "bl"      = numeric(),
                      "site"    = numeric(), 
                      "jsd"     = numeric(),
                      "entropy" = numeric(),
                      "diffsum" = numeric(),
                      "dnds"    = numeric(),
                      "method"  = factor())
for (name in dms.datasets)
{    
    true <- read.csv(paste0(truedir, name, "_true_dnds_entropy.csv"))
    for (bl in branch.lengths)
    { 
        for (meth in methods)
        {  
            dat <- read.csv(paste0(datadir, name, "_bl", bl, "_", meth, "_statistics.csv"))
            l <- nrow(dat)
            temp <- data.frame("dataset" = name, "bl" = bl, "site" = dat$site, "jsd" = dat$jsd, "entropy" = dat$entropy, "diffsum" = dat$abs_sum_differences, "dnds" = dat$dnds, "method" = meth, "true.dnds" = true$dnds, "true.entropy" = true$entropy)
            dat.dms <- rbind(dat.dms, temp)
        }
        
    }
}
write.csv(dat.dms, paste0(datadir, "dms_results.csv"), row.names=FALSE, quote=FALSE)


# 
# 
# print("Processing selection coefficients, yeast")
# for (dataset in yeast.datasets)
# {
#     for (bl in branch.lengths)
#     {
#         for (del in deltypes){
#             sim.selcoeffs <- data.frame("dataset" = character(), "bl" = numeric(), "method" = character(), "dummy" = numeric(), "del" = character(), "binnedcoeff" = numeric(), "realcoeff" = numeric())
#             full_dataset <- paste0(dataset, "_del", del)
#             for (m in c("true", methods))
#             {
#                 if (m == "true")
#                 {
#                     directory <- truedir
#                 }
#                 else
#                 {
#                     directory <- datadir
#                 }
#             
#                 dat <- read.csv(paste0(directory, full_dataset, "_bl", bl, "_", m, "_selcoeffs.csv"))
#                 temp <- data.frame("dataset" = dataset, "bl" = bl, "method" = m, "dummy" = 1:nrow(dat), "del" = del, "binnedcoeff" = dat$binnedcoeff, "realcoeff" = dat$realcoeff)
#                 sim.selcoeffs <- rbind(sim.selcoeffs, temp)
#             }
#             write.csv(sim.selcoeffs, paste0(datadir, full_dataset, "_selection_coefficients.csv"), row.names=FALSE, quote=FALSE)
#         }
#     }
# }
# 
# print("Processing selection coefficients, deep mutational scanning")
# for (dataset in dms.datasets)
# {
#     for (bl in branch.lengths)
#     {
#         sim.selcoeffs <- data.frame("dataset" = character(), "bl" = numeric(), "method" = character(), "dummy" = numeric(), "binnedcoeff" = numeric(), "realcoeff" = numeric())
#         for (m in c("true", methods))
#         {
#             if (m == "true")
#             {
#                 directory <- truedir
#             }
#             else
#             {
#                 directory <- datadir
#             }
#             
#             dat <- read.csv(paste0(directory, dataset, "_bl", bl, "_", m, "_selcoeffs.csv"))
#             temp <- data.frame("dataset" = dataset, "bl" = bl, "method" = m, "dummy" = 1:nrow(dat), "binnedcoeff" = dat$binnedcoeff, "realcoeff" = dat$realcoeff)
#             sim.selcoeffs <- rbind(sim.selcoeffs, temp)
#         }
#     }
#     write.csv(sim.selcoeffs, paste0(datadir, full_dataset, "_selection_coefficients.csv"), row.names=FALSE, quote=FALSE)
# }



# 
# kl.breaks <- seq(from=-10, to=10, by=1)
# 
# for (dataset in yeast.datasets)
# {
#     print(dataset)
#     true <- read.csv(paste0("../../simulation/true_simulation_parameters/", dataset, "_delstrong_true_selcoeffs.csv"))
#     true.hist <- hist(true$binnedcoeff, breaks = kl.breaks)
#     
#     for (m in c("nopenal", "mvn10", "d0.01", "d0.01", "phylobayes"))
#     {
#         dat <- read.csv(paste0(dataset, "_delstrong_bl0.5_", m, "_selcoeffs.csv"))
#         dat.hist <- hist(dat$binnedcoeff, breaks = kl.breaks)
#         print(m)
#         kl <- kl.divergence(true.hist, dat.hist)
#         print(kl)
#         stop()
#     }
# }
# 
