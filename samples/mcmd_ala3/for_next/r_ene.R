args <- commandArgs(trailingOnly = T) 
n.stages <- as.integer(args[1])
n.series <- as.integer(args[2])
fn.out <- args[3]

dt <- 1.0

read.ene <- function(series, stage, skip=1){
  label.ser <- paste("n",series,sep="")
  path.ene <- paste(stage, label.ser, "mult.ene", sep="/")
  dat.ene.raw <- read.table(path.ene)
  dat.ene.raw <- dat.ene.raw[seq(1, nrow(dat.ene.raw), skip), ]
  n.data <- length(dat.ene.raw)
  dat.ene.lab <- data.frame(seq(1:n.data)*skip*dt, dat.ene.raw, rep(label.ser, n.data), rep(as.character(stage), n.data))
  colnames(dat.ene.lab) <- c("time", "value", "series", "stage")
  dat.ene.lab
}

library(ggplot2)
## functions for generating plot

read.ene.stage <- function(stage, series_set, skip=1){
  ene.ser <- read.ene(series_set[1], stage, skip)
  for ( i_series in 2:length(series_set) ){
    tmp.ene.ser <- read.ene(series_set[i_series], stage, skip)
    ene.ser <- rbind(ene.ser, tmp.ene.ser)
  }
  ene.ser
}
read.ene.all <- function(stage_set, series_set, skip=1){
  ene.stg <- read.ene.stage(stage_set[1], series_set, skip)
  last.time <- max(ene.stg$time)
  if ( length(stage_set) > 1){
    for (i_stg in 2:length(stage_set)){
      tmp.ene <- read.ene.stage(stage_set[i_stg], series_set, skip)
       tmp.ene$time <- tmp.ene$time + last.time
       last.time <- max(tmp.ene$time)
      ene.stg <- rbind(ene.stg, tmp.ene)
     }
  }
  ene.stg
}
plot.ene <- function(ene.ser){
  p <- ggplot(ene.ser)
  fn.line <- paste(fn.out, sep="")
  p.line <- p + geom_line(aes(time,value, color=series))
  ggsave(file=fn.line, plot=p.line, dpi=600, width=6.4, height=4.8)
  return(p.line)
}

ene <- read.ene.all(1:n.stages, 1:n.series, 1000)

p <- plot.ene(ene)

