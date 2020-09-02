args <- commandArgs(trailingOnly = T) 
n.stages <- as.integer(args[1])
fn.out <- args[2]

stages <- 1:n.stages
read_pdf <- function(stage){
  fn <- paste(stage,"/dden.dat", sep="")
  tbl.2 <- c()
  if(file.access(fn) == 0){
    tbl.1<- read.table(fn)
    tbl.2 <- cbind(tbl.1, as.factor(rep(stage, nrow(tbl.1))) )
    colnames(tbl.2) <- c("E","density","stage")
  }
  tbl.2
}
pdf <- c()
pdf <- read_pdf(stages[1])
for (i in stages[c(2:length(stages))]){
    pdf.tmp <- read_pdf(i)
    pdf <- rbind(pdf, pdf.tmp)
}
library(ggplot2)
(p <- ggplot(pdf, aes(E, density, color=stage, group=stage)) + geom_point())
ggsave(file=fn.out, plot=p, dpi=600, width=6.4, height=4.8)