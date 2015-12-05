
require(dplyr)
require(tidyr)

results_dir    <- "../../results/"
datadir <- paste0(results_dir,"raw_results/")
datasets   <- c("PF00106", "PF00149", "PF00188", "PF00300", "PF00512", "PF00753", "PF01261", "PF01551", "PF01636", "PF03144", "PF00141", "PF00158", "PF00226", "PF00482", "PF00520", "PF01061", "PF01546", "PF01584", "PF02775", "PF04542", "PF00168", "PF00486", "PF00535", "PF01590", "PF03466", "PF00271", "PF00501", "PF00571", "PF00593", "PF00126", "PF01266", "PF01336", "PF01926", "PF02518", "PF04055", "PF07715")#, "pb2")
methods    <- c("nopenal") #, "d0.01", "d0.1", "d1.0", "mvn10", "mvn100", "mvn1000", "phylobayes")

dat <- data.frame("dataset" = factor(), 
                      "sw_emp"    = numeric(),
                      "sw_sim" = numeric(),
                      "slac_emp"    = numeric(),
                      "slac_sim"  = numeric())


for (dataset in datasets){
    
    # swmutsel, empirical
    dat2 <- read.table(paste0(datadir, "empirical/derived_dnds/", dataset, "_nopenal_dnds.txt"))
    sw_emp <- dat2$V1
    
    # swmutsel, simluation
    dat2 <- read.table(paste0(datadir, "simulated_from_empirical/derived_dnds/", dataset, "_simulated_dnds.txt"))
    sw_sim <- dat2$V1
        
    # SLAC, empirical 
    dat2 <- read.table(paste0(datadir, "empirical/slac/", dataset, "_SLAC.txt"), header=T)
    slac_emp <- dat2$dN/mean(dat2$dS)

    # SLAC, simulated
    dat2 <- read.table(paste0(datadir, "simulated_from_empirical/slac/", dataset, "_simulated_SLAC.txt"), header=T)
    slac_sim <- dat2$dN/mean(dat2$dS)

    temp <- data.frame("dataset" = dataset, "sw_emp" = sw_emp, "sw_sim" = sw_sim, "slac_emp" = slac_emp, "slac_sim" = slac_sim)
    dat <- rbind(dat, temp)    
}



for (d in dat3$dataset) {
  dat %>% filter(dataset == d) -> subdat
  fit <- lm(sw_emp~slac_emp,data=subdat)
  empslope <- round(empfit$coefficients[[2]],3)
  empcor
  
  simfit <- lm(sw_sim~slac_sim,data=subdat)
  simslope <- round(simfit$coefficients[[2]],3)
  
  
  p1 <- ggplot(subdat, aes(x = slac_emp, y=sw_emp))+geom_point(size=3)+geom_smooth(method="lm",color="red") + geom_abline(slope=1,intercept=0)+ggtitle("empirical") + annotate("text",x=0.1,y=0.9,label=empslope)+scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,1))
  p2 <- ggplot(subdat, aes(x = slac_sim, y=sw_sim))+geom_point(size=3)+geom_smooth(method="lm",color="red")+geom_abline(slope=1,intercept=0)+ggtitle("simulated from emp.") + annotate("text",x=0.1,y=0.9,label=simslope)+scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,1))
  pfull <- plot_grid(p1,p2,nrow=1)
  print(pfull)
  readline()
  
}
