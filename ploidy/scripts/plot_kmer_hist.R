library(tidyverse)

KmerHistGen <- function(prefix, dirpath, coverages = NULL) {
  if (!is.null(coverages)) {
    mean_coverage = round(coverages[coverages$strain == basename(prefix),]$mean_coverage, 2)
    median_coverage = round(coverages[coverages$strain == basename(prefix),]$median_coverage, 2)
  }
  hist = read.table(file.path(dirpath, paste(prefix, ".khist", sep="")), col.names = c("depth", "rawcount", "count"), sep="\t")
  peaks = read.table(file.path(dirpath, paste(prefix, ".peaks", sep="")), comment.char = "#", col.names=c("start","center","stop","max","volume"), sep="\t")
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

files = list.files(args[1])

# Get prefixes only
prefixes = c()
for (f in files) {
  if (endsWith(f, "khist")) {
    p = gsub("[.]khist", "", f)
    prefixes = c(prefixes, p)
  }
}

for (p in prefixes) {
  this_hist = KmerHistGen(p, args[1])
  ggsave(filename = file.path(getwd(), paste(p, ".khist", ".pdf", sep="")), plot = this_hist, width = 8.5, height = 11, units = "in", device="pdf")
}
