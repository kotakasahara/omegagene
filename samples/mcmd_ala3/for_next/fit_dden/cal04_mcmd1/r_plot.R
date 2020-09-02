args <- commandArgs(trailingOnly = T) 
fn.out <- args[1]

library(ggplot2)
read.data <- function(){
fn.dat <- paste("e1_fort.11", sep="")
dat <- read.table(fn.dat)
fn.fit <- paste("e1_fort.12", sep="")
fit <- read.table(fn.fit)
dat.p <- cbind(dat, rep(as.factor(1), nrow(dat)), rep("plot", nrow(dat)))
fit.p <- cbind(fit, rep(as.factor(1), nrow(fit)), rep("fit", nrow(fit)))
colnames(dat.p) <- c("ene","dens","vs","plot")
colnames(fit.p) <- c("ene","dens","vs","plot")
rbind(dat.p, fit.p)
}
data <- read.data()
p <- ggplot(data, aes(ene, dens))
p <- p+ geom_line(data=subset(data, plot=="fit"), aes(color=plot))
p <- p+ geom_point(data=subset(data, plot=="plot"), aes(color=plot))

ggsave(file=fn.out, plot=p, dpi=600, width=6.4, height=4.8)