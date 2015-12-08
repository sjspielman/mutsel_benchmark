
require(dplyr)
require(tidyr)
require(cowplot)
# results_dir    <- "../../results/"
# datadir <- paste0(results_dir,"raw_results/")
# datasets   <- c("PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")#, "pb2")
# methods    <- c("nopenal","d1.0")
# 
# simdat <- data.frame("dataset" = character(), 
#                      "site" = numeric(),
#                   "dnds"    = numeric(),
#                   "method" = character())
#                   
# empdat <- data.frame("dataset" = character(), 
#                      "site" = numeric(),
#                   "dnds"    = numeric(),
#                   "method" = character())
# 
# for (dataset in datasets){
#     
#     # swmutsel empirical, i.e. true
#     dat2 <- read.table(paste0(datadir, "empirical/derived_dnds/", dataset, "_nopenal_dnds.txt"))
#     temp <- data.frame("dataset" = dataset, "site"=1:nrow(dat2), "dnds" = dat2$V1, "method" = "true")
#     empdat <- rbind(empdat, temp)    
# 
#     
#     # swmutsel, simluation
#     for (m in methods){
#         dat2 <- read.table(paste0(datadir, "simulated_from_empirical/derived_dnds/", dataset, "_simulated_", m, "_dnds.txt"))
#         temp <- data.frame("dataset" = dataset, "site"=1:nrow(dat2), "dnds" = dat2$V1, "method" = m)
#         simdat <- rbind(simdat, temp) 
#     }   
# 
#     # SLAC, empirical 
#     dat2 <- read.table(paste0(datadir, "empirical/slac/", dataset, "_SLAC.txt"), header=T)
#     slac_emp <- dat2$dN/mean(dat2$dS)
#     temp <- data.frame("dataset" = dataset, "site"=  1:length(slac_emp), "dnds" = slac_emp, "method" = "slac")
#     empdat <- rbind(empdat, temp) 
#     
#     # SLAC, simulated
#     dat2 <- read.table(paste0(datadir, "simulated_from_empirical/slac/", dataset, "_simulated_SLAC.txt"), header=T)
#     slac_sim <- dat2$dN/mean(dat2$dS)
#     temp <- data.frame("dataset" = dataset, "site"=  1:length(slac_sim), "dnds" = slac_sim, "method" = "slac")
#     simdat <- rbind(simdat, temp) 
#    
# }
# 


for (d in datasets) {
  empdat %>% filter(dataset == d) %>% spread(method,dnds) -> subemp
  simdat %>% filter(dataset == d) %>% spread(method,dnds) -> subsim
  newdf <- data.frame("true" = subemp$true, "nopenal" = subsim$nopenal, "slac" = subsim$slac, "d1.0" = subsim$d1.0)
 
  fit <- lm(true~slac,data = subemp)
  m <- round(fit$coefficients[[2]],3)
  p1 <- ggplot(subemp, aes(x = slac,y=true)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("empirical nopenal-slac") + ylab("nopenal") + xlab("slac") + scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))

  fit <- lm(nopenal~slac,data = subsim)
  m <- round(fit$coefficients[[2]],3)
  p2 <- ggplot(subsim, aes(x = slac,y=nopenal)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("simulation nopenal-slac") + ylab("nopenal") + xlab("slac")+ scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))

  newnewdf <- data.frame("sim" = subsim$nopenal, "emp" = subemp$true)
  fit <- lm(emp~sim,data = newnewdf)
  m <- round(fit$coefficients[[2]],3)
  p3 <- ggplot(newnewdf, aes(x = sim,y=emp)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("empirical-simulation nopenal") + ylab("emp") + xlab("sim")+ scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))



 
#   fit <- lm(slac~true,data = newdf)
#   m <- round(fit$coefficients[[2]],3)
#   p3 <- ggplot(newdf, aes(x = true,y=slac)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("simulation slac-true") + ylab("slac") + xlab("true")+ scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))
#  
#   fit <- lm(nopenal~true,data = newdf)
#   m <- round(fit$coefficients[[2]],3)
#   p4 <- ggplot(newdf, aes(x = true,y=nopenal)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("simulation nopenal-true") + ylab("nopenal") + xlab("true")+ scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))
  
#   fit <- lm(d1.0~true,data = newdf)
#   m <- round(fit$coefficients[[2]],3)
#   p5 <- ggplot(newdf, aes(x = true,y=d1.0)) + geom_point(size=3) + geom_smooth(method = "lm") + geom_abline(slope=1, intercept=0) + annotate("text",x=0.1,y=0.9,label=m) + ggtitle("simulation d1.0-true") + ylab("d1.0") + xlab("true")+ scale_y_continuous(limits=c(0,1))+ scale_x_continuous(limits=c(0,1))

  pfull <- plot_grid(p1,p2,p3,nrow=1)
  print(pfull)
  readline() 
}