args <- commandArgs(trailingOnly = T) 
n.stages <- as.integer(args[1])
n.series <- as.integer(args[2])
fn.out <- args[3]

dt <- 1.0

read.ene.bin <- function(series, stage, skip=1, endi="little"){
  label.ser <- paste("n",series,sep="")
  path.ene <- paste(stage, label.ser, "mult.ene", sep="/")
  f <- file(path.ene, "rb")
  magic <- readBin(f, integer(), n=1, size=4L, signed=T, endian=endi)
  size <- readBin(f, integer(), n=1, size=4L, signed=T, endian=endi)
  col <- readBin(f, integer(), n=1, size=4L, signed=T, endian=endi)
  dat.ene.raw <- raw(0)
  while ( length(buf <- readBin(f, double(), n=10000, size=8L, signed=T, endian=endi) ) >0){
    dat.ene.raw <- c(dat.ene.raw, buf)
  }

  dat.ene.raw <- dat.ene.raw[seq(1, length(dat.ene.raw), skip)]
  n.data <- length(dat.ene.raw)
  dat.ene.lab <- data.frame(seq(1:n.data)*skip*dt, dat.ene.raw, rep(label.ser, n.data), rep(as.character(stage), n.data))
  colnames(dat.ene.lab) <- c("time", "value", "series", "stage")
  close(f)
  dat.ene.lab
}

library(ggplot2)
## functions for generating plot

read.ene.stage <- function(stage, series_set, skip=1){
  ene.ser <- read.ene.bin(series_set[1], stage, skip)
  for ( i_series in 2:length(series_set) ){
    tmp.ene.ser <- read.ene.bin(series_set[i_series], stage, skip)
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

