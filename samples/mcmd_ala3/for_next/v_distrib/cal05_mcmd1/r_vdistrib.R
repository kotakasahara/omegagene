args <- commandArgs(trailingOnly = T) 
n.vs <- as.integer(args[1])
fn.out <- args[2]
vs <- 1:n.vs
read_pdf <- function(vsnum){
  fn <- paste("v_pdf/s",vsnum,".pdf", sep="")
  tbl.2 <- c()
  if(file.access(fn) == 0){
    tbl.1<- read.table(fn)
    tbl.2 <- cbind(tbl.1, as.factor(rep(vsnum, nrow(tbl.1))) )
    colnames(tbl.2) <- c("E","density","vs")
  }
  tbl.2
}
pdf <- c()
pdf <- read_pdf(vs[1])
for (i in vs[c(2:length(vs))]){
    pdf.tmp <- read_pdf(i)
    pdf <- rbind(pdf, pdf.tmp)
}
library(ggplot2)
p <- ggplot(pdf, aes(E, density, color=vs, group=vs)) + geom_line()
ggsave(file=fn.out, plot=p, dpi=600, width=6.4, height=4.8)
