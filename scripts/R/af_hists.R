library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(reshape2)
library(magrittr)
library(stringr)
library(plotrix)
source("~/scripts/rbind.na.r")
source("~/scripts/cbind.na.r")

pull_strain <- function(fname) {
  strsplit(fname, "[.]")[[1]][1]
}

GatkAfHistGen <- function(path, isolate, dpfilt = FALSE) {
  if (class(path) == "data.frame") {
    #isolate = unlist(path[isolate,]$strain)
    print(isolate)
    path = path[path$strain == isolate,]$fname
    print(path)
  }
  print(df)
  df = read.table(path, header=FALSE)
  df %<>% 
    separate("V1", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
    filter(is.na(DELETE)) %>% 
    select(-DELETE,-NON) %>%
    mutate_all(as.numeric) %>%
    mutate(AFp = p/(p+q), AFq = q/(p+q), cov = p+q) %>%
    filter (AFp != 0.00 & AFp != 1.00) %>%
    mutate_all(as.numeric)

  if (dpfilt == TRUE) {
    cov_summ = summary(df$cov)
    df %<>%
      filter(cov >= cov_summ["1st Qu."] & cov <= cov_summ["3rd Qu."])
  }
  
  df = as.data.frame(cbind.na(as.matrix(df), rbind(as.matrix(df$AFp), as.matrix(df$AFq))))
  colnames(df)[6] = "combined"
  print(df)
  ggplot(data=df) +
    geom_histogram(aes(x=AFp), fill = "blue", colour="black", alpha=0.25, binwidth=0.01) +
    geom_histogram(aes(x=AFq), fill = "red", colour="black", alpha=0.25, binwidth=0.01) +
    geom_histogram(aes(x=combined), fill = "darkgreen", colour="black", alpha=0.25, binwidth=0.01) +
    #scale_fill_manual(name = "Coverage Bins", values=c("x < 60" = "blue", "60 < x < 120" = "red", "x > 120" = "green")) +
    annotate("text", x=Inf, y=Inf, label=paste(dim(df)[1]/2, "SNPs", sep=" "), vjust=1.2, hjust=1.2, cex=5) +
    xlab("Allele Frequency") +
    ggtitle(paste(isolate, ", GATK SNP Allele Frequencies", sep=""))
}
GatkAfHistGen("~/work/pursuit_paper/ploidy/", "JEL0388")
GatkAfHistGen("~/work/pursuit_paper/ploidy/JEL0388_stats.tsv", "JEL0388", dpfilt = TRUE)
af = GatkAfHistGen("~/work/pursuit_paper/ploidy/vcf_AD/JEL0728_consensus_contigs.fasta_clip_500.vcf.AD", "JEL0728", dpfilt = TRUE)

MeanStatPlot <- function(path, isolate, dpfilt = FALSE) {
  if (class(path) == "data.frame") {
    #isolate = unlist(path[isolate,]$strain)
    print(isolate)
    path = path[path$strain == isolate,]$fname
    print(path)
  }
  df = read.table(path, header=FALSE)
  df %<>% 
    separate("V2", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
    filter(is.na(DELETE)) %>% 
    select(-DELETE,-NON) %>%
    mutate_all(as.numeric) %>%
    mutate(AFp = p/(p+q), AFq = q/(p+q), cov = p+q) %>%
    filter (AFp != 0.00 & AFp != 1.00) %>%
    mutate_all(as.numeric) %>%
    rename(MQRankSum = V1) %>%
    mutate(AFp_round = round(AFp, 2))
  
  if (dpfilt == TRUE) {
    cov_summ = summary(df$cov)
    df %<>%
      filter(cov >= cov_summ["1st Qu."] & cov <= cov_summ["3rd Qu."])
  }
  
  df %<>%
    group_by(AFp_round) %>%
    summarise (
      MQRankSum_mean = mean(MQRankSum),
      MQRankSum_stderr = std.error(MQRankSum),
      count=n())
  
  mq = ggplot(data=df, aes(x=AFp_round, MQRankSum_mean)) + 
    geom_line() +
    geom_errorbar(aes(ymin=MQRankSum_mean-MQRankSum_stderr, ymax=MQRankSum_mean+MQRankSum_stderr),width=0.025, linetype="solid", col="black") +
    geom_point(aes(color = count), size=2) +
    scale_color_gradient(low="blue", high="red") +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept = 0.50, linetype="dashed") +
    guides(
      color = FALSE
    ) +
    xlim(0,1)
  
  mq
}


###### prep fname-strain DF ####
fs = as.data.frame(list.files("~/work/omics/pursuit/ploidy/MQRS_AD", pattern="*.MQRS_AD.tsv"))
colnames(fs) = c("fname")
fs %<>%
  mutate_all(as.character) %>%
  mutate(strain = unlist(lapply(fname, pull_strain))) %>%
  mutate(fname = paste("~/work/omics/pursuit/ploidy/MQRS_AD", fname, sep="/"))

kh = as.data.frame(list.files("~/work/omics/pursuit/ploidy/kmerhist/", pattern="*peaks"))
colnames(kh) = c("strain")
kh %<>%
  mutate_all(as.character) %>%
  mutate(strain = gsub("_23[.]peaks", "", strain)) %>%
  mutate(kmerhist_prefix = paste("~/work/omics/pursuit/ploidy/kmerhist/", strain, sep="/"))

fs %<>%
  left_join(kh, by="strain")
#####
coverages = read.table("~/work/omics/pursuit/ploidy/kmerhist/all.contigs.mean.tsv", sep="\t")
colnames(coverages) = c("contig", "depth", "strain")
coverages %<>% 
  select (strain, depth) %>%
  mutate (strain = as.character(strain)) %>%
  group_by (strain) %>%
  summarise(mean_coverage = mean(depth), median_coverage = median(depth))
####
AFpHistGen("~/work/omics/pursuit/ploidy/MQRS_AD/JEL0888.L50.MQRS_AD.tsv", isolate="JEL085", dpfilt = TRUE)

MQRS_plts = lapply(fs$strain, MeanStatPlot, path=fs, dpfilt=TRUE)
AFhists = lapply(fs$strain, AFpHistGen, path=fs, dpfilt=TRUE)
Khists = lapply(fs$kmerhist_prefix, KmerHistGen, coverages=coverages)

dev.off()
pdf(
  file=paste("~/work/omics/pursuit/ploidy/batch_ploidy_plots/batch_plots_47.pdf", sep=""),
  width=8.5,
  height=11
  )
for (i in 1:47) {
  plt_pane(
    list(AFhists[[i]], Khists[[i]], MQRS_plts[[i]]),
    rows = 2,
    cols =2
  )
}
dev.off()


AFpHistGen <- function(path, isolate, dpfilt = FALSE) {
  if (class(path) == "data.frame") {
    #isolate = unlist(path[isolate,]$strain)
    print(isolate)
    path = path[path$strain == isolate,]$fname
    print(path)
  }
  df = read.table(path, header=FALSE)
  df %<>% 
    filter(V1 == 0.00) %>%
    separate("V2", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
    filter(is.na(DELETE)) %>% 
    select(-DELETE,-NON) %>%
    mutate_all(as.numeric) %>%
    mutate(AFp = p/(p+q), AFq = q/(p+q), cov = p+q) %>%
    filter (AFp != 0.00 & AFp != 1.00) %>%
    mutate_all(as.numeric)
  
  if (dpfilt == TRUE) {
    cov_summ = summary(df$cov)
    df %<>%
      filter(cov >= cov_summ["1st Qu."] & cov <= cov_summ["3rd Qu."])
  }
  print(df)
  
  #df = as.data.frame(cbind.na(as.matrix(df), rbind(as.matrix(df$AFp), as.matrix(df$AFq))))
  #colnames(df)[6] = "combined"

  ggplot(data=df) +
    geom_histogram(aes(x=AFp), fill = "blue", colour="black", alpha=0.25, binwidth=0.01) +
    #geom_histogram(aes(x=AFq), fill = "red", colour="black", alpha=0.25, binwidth=0.01) +
    #geom_histogram(aes(x=combined), fill = "darkgreen", colour="black", alpha=0.25, binwidth=0.01) +
    #scale_fill_manual(name = "Coverage Bins", values=c("x < 60" = "blue", "60 < x < 120" = "red", "x > 120" = "green")) +
    annotate("text", x=Inf, y=Inf, label=paste(dim(df)[1], "SNPs", sep=" "), vjust=1.2, hjust=1.2, cex=5) +
    xlim(0,1) +
    xlab("Allele Frequency") +
    ggtitle(isolate)
}
AFpHistGen("~/work/pursuit_paper/ploidy/JEL0388_stats.tsv", "JEL0388", dpfilt = TRUE)
###########################
###Batch GATK histograms###
###########################
abspath="/home/aimzez/work/pursuit_paper/ploidy/vcf_AD_L50/"
gatk_todo = as.data.frame(list.files(abspath, pattern="*.vcf.AD")) %>%
  mutate_all(as.character)
colnames(gatk_todo) = c("fname")
gatk_todo$strain = lapply(gatk_todo$fname, pull_strain)
gatk_todo$fname = paste(abspath, gatk_todo$fname, sep="")
gatk_todo %<>% 
  mutate(fname = as.character(fname))

gatk_hists = lapply(unlist(gatk_todo$strain), GatkAfHistGen, path = gatk_todo, dpfilt=TRUE)
plt_pane(plt_list = gatk_hists[1:9], rows=3, cols=3)

##########################

PilonAfHistGen <- function(path, isolate) {
  if (class(path) == "data.frame") {
    #isolate = unlist(path[isolate,]$strain)
    print(isolate)
    path = path[path$strain == isolate,]$fname
    print(path)
  }
  df = read.table(path, sep="\t", header=FALSE, col.names = c("p"))
  df %<>% 
    filter (p != 0.00 & p != 1.00) %>%
    mutate (q = 1.00-p)
  print(df)
    
  df = as.data.frame(cbind.na(as.matrix(df), rbind(as.matrix(df$p), as.matrix(df$q))))
  colnames(df)[3] = "combined"
  ggplot(df) +
    geom_histogram(aes(x=p), binwidth=0.01, color="black", fill="blue", alpha=0.3) +
    geom_histogram(aes(x=q), binwidth=0.01, color="black", fill="red", alpha=0.3) +
    geom_histogram(aes(x=combined), binwidth=0.01, color="black", fill="darkgreen", alpha=0.3) +
    xlab("Allele Frequency") +
    ylab("count") +
    xlim(0,1) +
    annotate("text", x=Inf, y=Inf, label=paste(dim(df)[1]/2, "SNPs", sep=" "), vjust=1.2, hjust=1.2, cex=5) +
    ggtitle(paste(isolate, ", Pilon-filtered SNP Allele Frequencies", sep=""))
}
PilonAfHistGen("~/work/pursuit_paper/ploidy/pilon_vcf_AF/ARG085_pilon.vcf.AF", "ARG085")

###########################
###Batch Pilon histograms###
###########################
abspath="/home/aimzez/work/pursuit_paper/ploidy/pilon_vcf_AF/"
pilon_todo = as.data.frame(list.files(abspath, pattern="*.vcf.AF")) %>%
  mutate_all(as.character)
colnames(pilon_todo) = c("fname")
pilon_todo$strain = lapply(pilon_todo$fname, pull_strain)
pilon_todo$fname = paste(abspath, pilon_todo$fname, sep="")
pilon_todo %<>% 
  mutate(fname = as.character(fname))

pilon_hists = lapply(unlist(pilon_todo$strain), PilonAfHistGen, path = pilon_todo)
plt_pane(plt_list = pilon_hists[1:9], rows=3, cols=3)

#########################

KmerHistGen <- function(prefix, coverages = NULL) {
  if (!is.null(coverages)) {
    mean_coverage = round(coverages[coverages$strain == basename(prefix),]$mean_coverage, 2)
    median_coverage = round(coverages[coverages$strain == basename(prefix),]$median_coverage, 2)
  }
  hist = read.table(paste(prefix, "_23.khist", sep=""), col.names = c("depth", "rawcount", "count"), sep="\t")
  peaks = read.table(paste(prefix, "_23.peaks", sep=""), comment.char = "#", col.names=c("start","center","stop","max","volume"), sep="\t")
  ggplot(hist) +
    geom_point(aes(x=depth, y=count)) + 
    geom_line(aes(x=depth, y=count)) +
    ylim(0,max(hist$count)*1.25) +
    #ylim(0,3e7)+
    xlim(0,min(peaks$stop)*3.00) +
    geom_vline(xintercept = peaks[1,]$center, colour = "blue", alpha=0.5) +
    geom_vline(xintercept = peaks[2,]$center, colour = "blue", alpha=0.5) +
    #annotate("text", x=Inf, y=Inf, label=paste("mean coverage = ", mean_coverage, "x", sep=""), vjust=1.2, hjust=1.2, cex=5) +
    #annotate("text", x=Inf, y=Inf, label=paste("median coverage = ", median_coverage, "x", sep=""), vjust=3.0, hjust=1.2, cex=5) +
    ggtitle(paste(basename(prefix), ", 23-mer histogram", sep=""))
}
KmerHistGen("~/work/omics/neozygites/khist/ARSEF_5376")
KmerHistGen("~/work/omics/neozygites/khist/NeoP_MDA")
KmerHistGen("~/work/omics/neozygites/khist/Neo_30")

KmerHistGen("~/DATA/pursuit/ploidy/pacbio_kmercount/Chyhya1")

##################################
#### Batch kmer hists ############
##################################
coverages = read.table("~/work/pursuit_paper/ploidy/depths/all.contig_means.tsv", sep="\t")
colnames(coverages) = c("contig", "depth", "strain")
coverages %<>% 
  select (strain, depth) %>%
  mutate (strain = as.character(strain)) %>%
  group_by (strain) %>%
  summarise(mean_coverage = mean(depth), median_coverage = median(depth))

abspath="/home/aimzez/work/pursuit_paper/ploidy/kmerhist/"
khist_todo = list.files(abspath, pattern="*.khist")
khist_todo = unlist(lapply(khist_todo, sub, pattern="_23.khist", replacement=""))
khist_todo = paste(abspath, khist_todo, sep="")
khists = lapply(khist_todo, KmerHistGen, coverages=coverages)
plt_pane(plt_list = khists[1:9], rows=3, cols=3)
plt_pane(plt_list = khists[10:19], rows=4, cols=3)

##################################

for (i in seq(1,21,3)) {
  plt_pane(list(khists[i][[1]], gatk_hists[i][[1]], #pilon_hists[i][[1]],
    khists[i+1][[1]], gatk_hists[i+1][[1]], #pilon_hists[i+1][[1]], 
    khists[i+2][[1]], gatk_hists[i+2][[1]]), #pilon_hists[i+2][[1]]),
    rows=3, cols=2)
}
plt_pane(list(khists[19][[1]], gatk_hists[19][[1]]), rows=3, cols=2)
         
data = read.table("~/work/pursuit_paper/ploidy/stats.JEL0388", header=1)
hist(data$QUAL, breaks=1000, xlim=c(0,2000), ylim=c(0,500))
abline(v=mean(data[,4]), col="red")
abline(v=85.6, col="blue")
abline(v=663.6, col = "blue")

af=read.table("~/work/pursuit_paper/ploidy/vcf_AD_pilon/JEL0388_pilon_rawsort.AF", header=FALSE)
af$V2 = 1 - af$V1
afb = cbind.na(af, rbind(as.matrix(af$V1),as.matrix(af$V2)))
colnames(afb)[3] = "V3"
ggplot(afb) +
  geom_histogram(aes(x=V1), binwidth = 0.01, colour = "black", fill="blue", alpha=0.4) +
  geom_histogram(aes(x=V2), binwidth = 0.01, colour="black", fill="red", alpha=0.4) +
  geom_histogram(aes(x=V3), binwidth = 0.01, colour="black", fill="darkgreen", alpha=0.4) +
  xlim(0.01,0.99) +
  xlab("Allele Frequency") +
  ggtitle("JEL0388 PILON-filtered Allele Frequences")

gt60 = read.table("~/work/pursuit_paper/ploidy/JEL0388_consensus_contigs.fasta_clip_500.DPFilter_gt60lt120.vcf.AD", header=FALSE)
lt60 = read.table("~/work/pursuit_paper/ploidy/JEL0388_consensus_contigs.fasta_clip_500.DPFilter_gt19lt60.vcf.AD", header=FALSE)
lt19 = read.table("~/work/pursuit_paper/ploidy/JEL0388_consensus_contigs.fasta_clip_500.DPFilter_lt19.vcf.AD", header=FALSE)
gt120 = read.table("~/work/pursuit_paper/ploidy/JEL0388_consensus_contigs.fasta_clip_500.DPFilter_gt120.vcf.AD", header=FALSE)

gt60 %<>% 
  separate("V1", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
  filter(is.na(DELETE)) %>% 
  select(-DELETE,-NON) %>%
  mutate_all(as.numeric) %>%
  mutate(AFp = p/(p+q), AFq = q/(p+q)) %>%
  mutate(group = "GT60")
gt60 %<>%
  select(-p,-q) %>% 
  mutate(idx = seq.int(nrow(gt60)))

gt60_melt = melt(gt60, id.vars = c("idx","group"))

lt60 %<>% 
  separate("V1", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
  filter(is.na(DELETE)) %>% 
  select(-DELETE,-NON) %>%
  mutate_all(as.numeric) %>%
  mutate(AFp = p/(p+q), AFq = q/(p+q)) %>%
  mutate(group = "LT60")
lt60 %<>%
  select(-p,-q) %>%
  mutate(idx = seq.int(nrow(lt60)))

ggplot(lt60, aes(x=AFp, y=p+q)) +
  geom_point(alpha=0.1) +
  xlim(0.01,0.99)

lt60_melt = melt(lt60, id.vars = c("idx","group"))

gt120 %<>% 
  separate("V1", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
  filter(is.na(DELETE)) %>% 
  select(-DELETE,-NON) %>%
  mutate_all(as.numeric) %>%
  mutate(AFp = p/(p+q), AFq = q/(p+q)) %>%
  mutate(group = "GT120")
gt120 %<>%
  select(-p,-q) %>%
  mutate(idx = seq.int(nrow(gt120)))
gt120_melt = melt(gt120, id.vars = c("idx","group"))

lt19 %<>% 
  separate("V1", into=c("p","q","NON", "DELETE"), sep = ",") %>% 
  filter(is.na(DELETE)) %>% 
  select(-DELETE,-NON) %>%
  mutate_all(as.numeric) %>%
  mutate(AFp = p/(p+q), AFq = q/(p+q)) %>%
  mutate(group = "LT19")
lt19 %<>%
  select(-p,-q) %>%
  mutate(idx = seq.int(nrow(lt19)))
lt19_melt = melt(lt19, id.vars = c("idx","group"))

plt = rbind(lt19_melt, gt60_melt, lt60_melt, gt120_melt)

ggplot(plt, aes(x=value)) +
  geom_histogram(data = subset(plt, group=="LT19"), aes(fill="x < 19"), colour="black", alpha=0.25, binwidth=0.01) +
  geom_histogram(data = subset(plt, group=="LT60"), aes(fill="x < 60"), colour="black", alpha=0.25, binwidth=0.01) +
  geom_histogram(data = subset(plt, group=="GT60"), aes(fill="60 < x < 120"), colour="black", alpha=0.25, binwidth=0.01) +
  geom_histogram(data = subset(plt, group=="GT120"), aes(fill="x > 120"), colour="black", alpha=0.25, binwidth=0.01) +
  scale_fill_manual(name = "Coverage Bins", values=c("x < 19" = "yellow", "x < 60" = "blue", "60 < x < 120" = "red", "x > 120" = "green")) +
  xlim(0.01,0.99) +
  annotate("text", x=0.80, y=1100, label="mean coverage = 43.38x") +
  ggtitle("JEL0388 Allele Frequency Histogram")

