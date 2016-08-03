# SJS
# Script to merge JSD, dN/dS, and selection coefficient inferences from swMutSel and Phylobayes into concise single data frame(s).

require(dplyr)
require(tidyr)
require(purrr)

resdir       <- "dataframes/"
truedir      <- "../simulation/true_simulation_parameters/"
datadir      <- "dataframes/"

# Merge true selection coefficients

files <- dir(path = truedir, pattern = "*selcoeffs.csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(truedir, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "dummy"), sep="_true_selcoeffs.csv") %>%
  dplyr::select(-dummy) -> true.selcoeffs
write_csv(true.selcoeffs, paste0(resdir, "true_selection_coefficients.csv"))


all.datasets   <- c("1B4T_A", "1RII_A", "1V9S_B", "1G58_B", "1W7W_B", "2BCG_Y", "2CFE_A", "1R6M_A", "2FLI_A", "1GV3_A", "1IBS_A", "HA", "NP", "LAC", "Gal4")
branch.lengths <- c("0.01", "0.5")
dms.datasets   <- c("HA", "NP", "LAC", "Gal4")
methods        <- c("nopenal", "d0.01", "d0.1", "mvn10", "mvn100" , "phylobayes")

dat <- data.frame("dataset"      = character(),
                  "bl"           = numeric(),  # Branch length
                  "site"         = numeric(),  # Site in protein
                  "jsd"          = numeric(),  # Jensen-Shannon distance from true fitnesses
                  "diffsum"      = numeric(),  # Sum of absolute differences between inferred and true fitnesses
                  "entropy"      = numeric(),  # Entropy of inferred AA frequency distribution
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
