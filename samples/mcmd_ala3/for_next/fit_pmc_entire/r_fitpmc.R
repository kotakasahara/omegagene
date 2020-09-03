args <- commandArgs(trailingOnly = T) 
mdnum <- args[1]
nvst <- as.integer(args[2])
iordmin <- as.integer(args[3])
iordmax <- as.integer(args[4])
fn.out1 <- args[5]
fn.out2 <- args[6]
library(ggplot2)

vstates <- 1:nvst

read.vs.data <- function(vs, dir="."){
  fn.dat <- paste(dir, paste("s",vs,"_d.derv.dat", sep=""), sep="/")
  fn.pdf <- paste(dir, paste("s",vs,".pdf", sep=""), sep="/")
  fn.f12 <- paste(dir, paste("b",vs,".fort.12", sep=""), sep="/")
  dat <- read.table(fn.dat)
  pdf <- read.table(fn.pdf)
  f12 <- read.table(fn.f12)
  dat.p <- cbind(dat, rep(as.factor(vs), nrow(dat)), rep("plot", nrow(dat)))
  pdf.p <- cbind(pdf, rep(as.factor(vs), nrow(pdf)), rep("pdf", nrow(pdf)))
  f12.p <- cbind(f12, rep(as.factor(vs), nrow(f12)), rep("f12", nrow(f12)))
  colnames(dat.p) <- c("ene","dens","vs","plot")
  colnames(pdf.p) <- c("ene","dens","vs","plot")
  colnames(f12.p) <- c("ene","dens","vs","plot")
  rbind(dat.p, pdf.p, f12.p)
}

read.iord <- function(mdnum, iord){
  dir = paste("cal", mdnum, "_mcmd1_", as.character(iord), sep="")
  data <- read.vs.data(1, dir)
  if (length(vstates) >= 2){
    for( i in vstates[2:length(vstates)] ) {
       data <- rbind(data, read.vs.data(i, dir) ) 
    }
  }
  data <- cbind(data, rep(iord, nrow(data)))
  colnames(data)[5] <- "iord"
  data
}

read.all <- function(mdnum, iords){
  data <- read.iord(mdnum, iords[1])
  for( i in iords[2:length(iords)]){
    d <- read.iord(mdnum, i)
    data <- rbind(data, d)
  }
  data
}

d <- read.all(mdnum, iordmin:iordmax)

 p1 <- ggplot(d, aes(ene, dens)) + geom_line(data=subset(d, plot=="f12"), aes(color=vs)) + geom_point(data=subset(d, plot=="pdf"), aes(color=vs), size=0.8) + facet_wrap(~iord, scales="free") 
 p2 <- ggplot(d, aes(ene, dens)) + geom_line(data=subset(d, plot=="plot"), aes(color=vs)) +facet_wrap(~iord, scales="free") 

ggsave(file=fn.out1, plot=p1, dpi=600, width=6.4, height=4.8)

ggsave(file=fn.out2, plot=p2, dpi=600, width=6.4, height=4.8)
