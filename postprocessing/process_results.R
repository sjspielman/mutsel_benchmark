# SJS
# Script to merge JSD, dN/dS, and selection coefficient inferences from swMutSel and Phylobayes into concise single data frame(s).

require(dplyr)
require(tidyr)

resdir       <- "dataframes/"
truedir      <- "../simulation/true_simulation_parameters/"
datadir      <- "dataframes/"

all.datasets   <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A", "HA", "NP", "LAC", "GAL4")
branch.lengths <- c("0.01", "0.5")
dms.datasets   <- c("HA", "NP", "LAC", "GAL4")
methods        <- c("nopenal", "d0.01", "d0.1", "mvn10", "mvn100" , "phylobayes")

dat <- data.frame("dataset"      = character(),
                  "bl"           = numeric(),  # Branch length 
                  "site"         = numeric(),  # Site in protein
                  "jsd"          = numeric(),  # Jensen-Shannon distance from true fitnesses
                  "diffsum"      = numeric(),  # Sum of absolute differences between inferred and true fitnesses
                  "entropy"      = numeric(),  # Entropy of inferred AA frequency distribution
                  "neff.aa"      = numeric(),  # Effective number of amino acids in the inferred distribution ( = exp(entropy) )
                  "dnds"         = numeric(),  # Predicted dN/dS of the inferred parameters
                  "method"       = factor(),   # Inference method
                  "true.dnds"    = numeric(),  # True dN/dS, of stationary distribution
                  "true.entropy" = numeric())  # True AA entropy, of stationary distribution
for (name in all.datasets)
{
    print(name)
    true <- read.csv(paste0(truedir, name, "_true_dnds_entropy.csv"))
    l <- nrow(dat)
    for (bl in branch.lengths)
    {
        for (meth in methods)
        {   
            if (bl == 0.5 & name == "LAC" & meth == "phylobayes") next
            results <- read.csv(paste0(datadir, name, "_bl", bl, "_", meth, "_statistics.csv"))
            temp <- data.frame("dataset" = name, "bl" = bl, "site" = results$site, "jsd" = results$jsd, "diffsum" = results$abs_sum_differences, "entropy" = results$entropy, "neff.aa" = exp(results$entropy), "dnds" = results$dnds, "method" = meth, "true.dnds" = true$dnds , "true.entropy" = true$entropy)
            dat <- rbind(dat, temp)
        }
    }
    
}
write.csv(dat, paste0(datadir, "inference_results.csv"), row.names=FALSE, quote=FALSE)


# 
# 
# print("Processing selection coefficients for KL divergence calculations")
# jsd.breaks <- seq(from=-10, to=10, by=0.5)
# jsd.data <- data.frame("dataset"  = character(),
#                         "del"    = character(),
#                         "bl"     = numeric(),
#                         "method" = numeric(),
#                         "jsd"     = numeric())
# for (dataset in yeast.datasets)
# {
#     print(dataset)
#     for (del in deltypes)
#     {
#         true <- read.csv(paste0(truedir, dataset, "_del", del, "_true_selcoeffs.csv"))
#         true.hist <- hist(true$binnedcoeff, breaks = jsd.breaks)
#         true.density <- true.hist$density
#         true.density[true.density == 0] <- 1e-20
#         
#         for (bl in branch.lengths){
#             
#             for (m in methods)
#             {
#             
#                 dat <- read.csv(paste0("dataframes/", dataset, "_del", del, "_bl", bl, "_", m, "_selcoeffs.csv"))
#                 dat.hist <- hist(dat$binnedcoeff, breaks = jsd.breaks)
#                 dat.density <- dat.hist$density
#                 dat.density[dat.density == 0] <- 1e-20
#                 mean.density <- (true.density + dat.density) / 2
#                 kl1 <- 0.5 * sum(true.density * log(true.density / mean.density))
#                 kl2 <- 0.5 * sum(dat.density * log(dat.density / mean.density))
#                 jsd <- sqrt(kl1 + kl2)
#                 temp <- data.frame("dataset" = dataset, "del" = del, "bl" = bl, "method" = m, "jsd" = jsd)
#                 jsd.data <- rbind(jsd.data, temp)
#             }
#         }
#     }
# }
# write.csv(jsd.data, paste0(datadir, "jsd_selcoeffs.csv"), row.names=FALSE, quote=FALSE)





