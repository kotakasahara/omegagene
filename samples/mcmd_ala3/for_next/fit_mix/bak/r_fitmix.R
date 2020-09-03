args <- commandArgs(trailingOnly = T) 
nvst <- as.integer(args[1])
fn.out <- args[2]

library(ggplot2)
 vstates <- 1:nvst
read.vs.data <- function(vs){
  #fn.dat <- paste("s",vs,".fort.11", sep="")
  #dat <- read.table(fn.dat)
  fn.fit <- paste("s",vs,".fort.12", sep="")
  fit <- read.table(fn.fit)
  #dat.p <- cbind(dat, rep(as.factor(vs), nrow(dat)), rep("plot", nrow(dat)))
  fit.p <- cbind(fit, rep(as.factor(vs), nrow(fit)), rep("fit", nrow(fit)))
  #colnames(dat.p) <- c("ene","dens","vs","plot")
  colnames(fit.p) <- c("ene","dens","vs","plot")
  #rbind(dat.p, fit.p)
  fit.p
}
data <- read.vs.data(1)
for( i in vstates[2:length(vstates)] ) {
  data <- rbind(data, read.vs.data(i) ) 
}

p <- ggplot(data, aes(ene, dens))
(p <- p+ geom_line(data=subset(data, plot=="fit"), aes(color=vs)))
ggsave(file=fn.out, plot=p, dpi=600, width=6.4, height=4.8)
