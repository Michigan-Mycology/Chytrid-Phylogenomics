library(tidyverse)
library(scales)

plt_pane<-function(plt_list, rows, cols) {
  n = length(plt_list)
  nper_row = n/rows
  sqftage = rows*cols
  blank_cell<-ggplot()+geom_blank()+theme(plot.background = element_rect(colour="white"),
                                          panel.background = element_rect(colour="white"))
  grob_list = lapply(plt_list, ggplotGrob)
  if (length(grob_list) < sqftage) {
    last<-length(grob_list)
    print(last)
    for (idx in (last+1):sqftage) {
      grob_list[[idx]]<-ggplotGrob(blank_cell)
    }
  }
  binding_expr = "rbind("
  i=1
  for (p in 1:rows) {
    columns<-rep(sprintf("grob_list[[%d]]",seq(i,(i+cols-1),1)))
    binding_expr<-str_c(binding_expr,
                        "cbind(",
                        str_c(columns,collapse=","),
                        ")"
    )
    if (p != rows) {
      binding_expr<-str_c(binding_expr, ",")
    }
    i=i+cols
  }
  binding_expr<-str_c(binding_expr,")")
  print(binding_expr)
  grid::grid.newpage()
  grid::grid.draw(eval(parse(text=binding_expr)))
}

GatkAfHistGen <- function(df, isolate) {
  ggplot(data=df) +
    geom_histogram(aes(x=p), fill = muted("blue"), colour="black", alpha=0.5, binwidth=0.01) +
    geom_histogram(aes(x=q), fill = muted("lightblue"), colour="black", alpha=0.5, binwidth=0.01) +
    #geom_histogram(aes(x=combined), fill = "darkgreen", colour="black", alpha=0.25, binwidth=0.01) +
    #scale_fill_manual(name = "Coverage Bins", values=c("x < 60" = "blue", "60 < x < 120" = "red", "x > 120" = "green")) +
    annotate("text", x=Inf, y=Inf, label=paste(dim(df)[1], "SNPs", sep=" "), vjust=1.2, hjust=1.2, cex=5) +
    xlab("Allele Frequency") +
    ggtitle(paste(isolate, ", GATK SNP Allele Frequencies", sep=""))
}

KmerHistGen <- function(prefix, coverages = NULL) {
  if (!is.null(coverages)) {
    mean_coverage = round(coverages[coverages$strain == basename(prefix),]$mean_coverage, 2)
    median_coverage = round(coverages[coverages$strain == basename(prefix),]$median_coverage, 2)
  }
  hist = read.table(paste(prefix, ".khist", sep=""), col.names = c("depth", "rawcount", "count"), sep="\t")
  peaks = read.table(paste(prefix, ".peaks", sep=""), comment.char = "#", col.names=c("start","center","stop","max","volume"), sep="\t")
  ggplot(hist) +
    geom_point(aes(x=depth, y=count)) + 
    geom_line(aes(x=depth, y=count)) +
    ylim(0,max(hist$count)*1.25) +
    xlim(0,min(peaks$stop)*3.00) +
    geom_vline(xintercept = peaks[1,]$center, colour = "blue", alpha=0.5) +
    geom_vline(xintercept = peaks[2,]$center, colour = "blue", alpha=0.5) +
    #annotate("text", x=Inf, y=Inf, label=paste("mean coverage = ", mean_coverage, "x", sep=""), vjust=1.2, hjust=1.2, cex=5) +
    #annotate("text", x=Inf, y=Inf, label=paste("median coverage = ", median_coverage, "x", sep=""), vjust=3.0, hjust=1.2, cex=5) +
    ggtitle(paste(basename(prefix), ", 23-mer histogram", sep=""))
}


args = commandArgs(trailingOnly = TRUE)
df_path = args[1]
isolate = args[2]
print(args)
df = read_delim(args[1], delim="\t") %>%
  filter (p != 0.00 & p != 1.00)

afhist = GatkAfHistGen(df, isolate)

pdf(file = paste(isolate, ".pdf", sep=""), width = 11, height = 8.5, onefile=FALSE)
if (!is.na(args[3])) {
  kmerhist = KmerHistGen(args[3])
  plt_list = list()
  plt_list[[1]] = kmerhist
  plt_list[[2]] = afhist
  plt_pane(plt_list, 1, 2)
} else {
  afhist
}
dev.off()
#ggsave(plot=plt, filename = paste(isolate, ".pdf", sep=""), device="pdf", width = 22, height = 17, units = "in")
