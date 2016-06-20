# Functions used in plot_figures.R


################### THIS FUNCTION IS FROM THE PACKAGE "SMATR", which cannot be loaded because it breaks my linear models. ################
slope.test <- function( y, x, test.value=1, V=matrix(0,2,2))
{

    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)

    resDF <- n - 2
    fCrit <- qf( 1-alpha, 1, resDF )

    dat <- cbind(y[iref], x[iref])
    r.factor <- 1
	  vr <- ( cov(dat) - V )*(n-1)
		vr <- t(dat)%*%dat - V*n
    r <- vr[1,2]/sqrt( vr[1,1]*vr[2,2] )

    bCI     <- matrix( NA, 1, 2 )
    varTest <- matrix( 0, 2, 2 )

     # linear regression code only. SJS
     b            <- vr[1,2]/vr[2,2]
     varRes       <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] )/resDF
     varB         <- varRes/vr[2,2] * r.factor
     bCI[1,1]     <- b - sqrt(varB)*sqrt(fCrit)
     bCI[1,2]     <- b + sqrt(varB)*sqrt(fCrit)
     varTest[1,1] <- vr[1,1] - 2*test.value*vr[1,2] + test.value^2*vr[2,2]
     varTest[1,2] <- vr[1,2] - test.value*vr[2,2]
     varTest[2,2] <- vr[2,2]
     rTest  <- varTest[1,2] / sqrt( varTest[1,1] ) / sqrt( varTest[2,2] )
     F      <- rTest^2/(1 - rTest^2)/r.factor*(n-2)
     pValue <- 1 - pf( F, 1, resDF)
     list( pValue, b )

}



# function to return pvalue from an lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Linear models on dnds, entropy inferences
summarize_dnds_entropy <- function(dat, siglevel){
  dat %>% group_by(dataset, type, method, bl) %>%
    do(rraw.dnds = cor(.$dnds, .$true.dnds),
       braw.dnds = glm(dnds ~ offset(true.dnds), dat=.),
       rraw.h = cor(.$entropy, .$true.entropy),
       braw.h = glm(entropy ~ offset(true.entropy), dat=.))   %>%
    mutate(r2.dnds   = rraw.dnds[[1]]^2,
           b.dnds = summary(braw.dnds)$coeff[[1]],
           b.dnds.sig = summary(braw.dnds)$coeff[4]<siglevel,
           r2.entropy = rraw.h[[1]]^2,
           b.entropy = summary(braw.h)$coeff[1],
           b.entropy.sig = summary(braw.h)$coeff[4]<siglevel) %>%
    select(-rraw.dnds, -rraw.h, -braw.dnds, -braw.h) -> dat.sum1
  dat %>% group_by(dataset, method) %>%
    do(x = slope.test(.$dnds, .$true.dnds, test.value=1)) %>%
    mutate(slope.dnds.sig = x[[1]]<siglevel, slope.dnds = x[[2]]) %>%
    select(-x) -> dat.sum2
  dat %>% group_by(dataset, method) %>%
    do(x = slope.test(.$entropy, .$true.entropy, test.value=1)) %>%
    mutate(slope.entropy.sig = x[[1]]<siglevel, slope.entropy = x[[2]]) %>%
    select(-x) -> dat.sum3
  part <- left_join(dat.sum2, dat.sum3)
  dat.sum <- left_join(dat.sum1, part)
  dat.sum
}





### Histogram plot grid ###
plot_sc_grid <- function(datasets){


  
}
for (d in natural.datasets){
  sc <- true.selcoeffs %>% filter(dataset == d) %>% mutate(method = "True")
  for (m in methods_levels){
    sc.temp <- read.csv(paste0("dataframes/",d,"_bl0.5_", m, "_selcoeffs.csv"))
    temp <- data.frame(dataset = d, binnedcoeffs = sc.temp$binnedcoeff, method = m)
    sc <- rbind(sc, temp)
  }
  sc$method <- factor(sc$method, levels = c("True", "nopenal", "mvn100", "mvn10", "d0.01", "d0.1", "phylobayes"), labels = c("True", "Unpenalized", "mvn100", "mvn10", "d0.01", "d0.1", "pbMutSel"))
  x <- ggplot(sc, aes(x = binnedcoeffs)) + geom_histogram(fill = "white", color = "black") + facet_grid(~method) + xlab("Selection Coefficients") + ylab("Count")
  x2 <- ggdraw(x) + draw_label(d, x = 0.965, y = 0.5, size=14, fontface = "bold")
  grid_list[[i]] <- x2
  i <- i + 1
}
