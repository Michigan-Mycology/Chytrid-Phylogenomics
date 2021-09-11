library(tidyverse)
library(scales)
library(patchwork)
library(readxl)

GatkAfHistGen <- function(df, isolate, approx_norm_binom_stdev = NULL, stats = NULL, title = isolate) {

  out_plt = ggplot(data=df) +
    geom_histogram(aes(x=p), fill = "black", colour="black", alpha=1.0, binwidth=0.01) +
    geom_histogram(aes(x=q), fill = "black", colour="black", alpha=1.0, binwidth=0.01) +
    #geom_histogram(aes(x=combined), fill = "darkgreen", colour="black", alpha=0.25, binwidth=0.01) +
    #scale_fill_manual(name = "Coverage Bins", values=c("x < 60" = "blue", "60 < x < 120" = "red", "x > 120" = "green")) +
    xlab("Allele Frequency") +
    ylab("Count") +
    ggtitle(title) +
    theme_bw() +
    theme (
      panel.grid = element_blank(),
      plot.title = element_text(size=8)
    )

  if (!is.null(approx_norm_binom_stdev)) {
    # Put the approximated normal distribution (of binomial) on the AF plot
    x = seq(0, 1, 0.01)
    y = dnorm(x = x, mean = 0.5, sd = approx_norm_binom_stdev)
    dist_points = as_tibble(cbind(x,y))
    dist_hist = hist(df$p, breaks = x)
    max_dist_hist = max(dist_hist$counts)
    scaling_factor = max_dist_hist/max(y)
    print(dist_points)
    dist_points = dist_points %>%
      mutate(y_scaled = y * scaling_factor)
    expect_range = c(0.5-approx_norm_binom_stdev, 0.5+approx_norm_binom_stdev)

    out_plt = out_plt +
      geom_line(data = dist_points, aes(x,y_scaled), color = "red") +
      geom_area(data = subset(dist_points, x>=expect_range[1] & x<=expect_range[2]), aes(x=x, y=y_scaled), fill = "red", alpha = 0.15)
  }

  yrange = ggplot_build(out_plt)$layout$panel_scales_y[[1]]$range$range[2]
  out_plt = out_plt +
    ylim(0, yrange*1.1)

  if (!is.null(stats)) {
    p_in_binom_expect = stats %>%
      filter(isolate == ploidy_file_prefix) %>%
      pull(p_in_binom_expect)

    snp_density = stats %>%
      filter(isolate == ploidy_file_prefix) %>%
      pull(snp_density)

    mean_coverage = stats %>%
      filter(isolate == ploidy_file_prefix) %>%
      pull(mean_coverage)

    out_plt = out_plt +
      annotate("text", x=Inf, y=Inf, label=paste(round(mean_coverage,1), "x coverage", sep =""), sep=" ", vjust=1.3, hjust=1.2, cex=2.5) +
      annotate("text", x=Inf, y=Inf, label=paste(dim(df)[1], "SNPs", sep=" "), vjust=2.7, hjust=1.2, cex=2.5) +
      annotate("text", x=Inf, y=Inf, label=paste(round(p_in_binom_expect,2), "in expect"), sep=" ", vjust=4.1, hjust=1.2, cex=2.5) +
      annotate("text", x=Inf, y=Inf, label=paste(round(snp_density,8), "SNP/bp"), sep=" ", vjust=5.5, hjust=1.1, cex=2.5)
  }

  out_plt
}

KmerHistGen <- function(prefix, coverages = NULL, title = basename(prefix)) {
  if (!is.null(coverages)) {
    mean_coverage = round(coverages[coverages$strain == basename(prefix),]$mean_coverage, 2)
    median_coverage = round(coverages[coverages$strain == basename(prefix),]$median_coverage, 2)
  }
  hist = read.table(paste(prefix, ".khist", sep=""), col.names = c("depth", "rawcount", "count"), sep="\t")

  #hist_quartiles = quantile(hist$count, prob = c(0.25, 0.5, 0.75))
  #hist_iqr = hist_quartiles[["75%"]] - hist_quartiles[["25%"]]

  #y_max = hist_quartiles[["75%"]] + ((1.5)*(hist_iqr))
  peaks = read.table(paste(prefix, ".peaks", sep=""), comment.char = "#", col.names=c("start","center","stop","max","volume"), sep="\t")
  lower_peak_center = min(peaks$center)
  print(lower_peak_center)
  print(peaks)
  peaks = peaks %>%
    arrange(desc(volume)) %>%
    top_n(2, wt = volume) %>%
    filter(center <= 10*lower_peak_center)
  print(peaks)
  #if ( (y_max < max(peaks$max)) | (y_max > max(peaks$max)*10 )) {
  #  y_max = max(peaks$max)*1.25
  #}
  ggplot(hist) +
    geom_point(aes(x=depth, y=count)) +
    geom_line(aes(x=depth, y=count)) +
    ylim(0,hist[peaks[1,2],3]*1.5) +
    xlim(0,max(peaks$stop)*2.0) +
    geom_vline(xintercept = peaks[1,]$center, colour = "blue", alpha=0.5) +
    geom_vline(xintercept = peaks[2,]$center, colour = "blue", alpha=0.5) +
    xlab("Depth") +
    ylab("Count") +
    #annotate("text", x=Inf, y=Inf, label=paste("mean coverage = ", mean_coverage, "x", sep=""), vjust=1.2, hjust=1.2, cex=5) +
    #annotate("text", x=Inf, y=Inf, label=paste("median coverage = ", median_coverage, "x", sep=""), vjust=3.0, hjust=1.2, cex=5) +
    ggtitle(title) +
    theme_bw()+
    theme (
      panel.grid = element_blank(),
      plot.title = element_text(size=8)
    )
}

#### To add stats that Tim wants on plots ####
CHYTRID_PHYLO = "~/dev/Chytrid-Phylogenomics"

l50_genome_sizes = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, l50_assembly_length)
isolates = read_xlsx(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Isolates.xlsx"))

PATH_PREFIX = "/scratch/amsesk/pursuit/ploidy"
binom_expect = read_delim(file.path(PATH_PREFIX, "all_AF_ranges_from_binom.tsv"), delim="\t", col_names=F) %>%
  rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  select(-ploidy_file_prefix) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

snp_densities = read_delim(file.path(CHYTRID_PHYLO, "figures/1_tree", "all.snp_contig_counts_strainified.tsv"), delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(num_snps = sum(num_snps)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  left_join(l50_genome_sizes) %>%
  left_join(binom_expect) %>%
  select(SPECIES.TREE.LABEL, ploidy_file_prefix, num_snps, l50_assembly_length, mean_coverage, p_in_binom_expect) %>%
  mutate(snp_density = num_snps/l50_assembly_length) %>%
  rename(label = SPECIES.TREE.LABEL)


#### Plot here ####
snp_stat_dir = "/scratch/amsesk/pursuit/ploidy/snp_stats/"
approx_norm_binom_stdev_tbl = "/scratch/amsesk/pursuit/ploidy/all_AF_ranges_from_binom_withApproxStDev.tsv"

approx_norm_binom_stdev_tbl = read_delim(approx_norm_binom_stdev_tbl, delim = "\t", col_names = F)
approx_norm_binom_stdev_tbl = approx_norm_binom_stdev_tbl %>%
  select(X1,X4) %>%
  rename(iso=X1, approx_norm_binom_stdev=X4)

snp_stat_files = list.files(snp_stat_dir)

ploidyfix2stl = read_xlsx("/scratch/amsesk/pursuit/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, ploidy_file_prefix) %>%
  mutate(SPECIES.TREE.LABEL = gsub("[_.]v[0-9][.]*[0-9]*", "", SPECIES.TREE.LABEL)) %>%
  mutate(SPECIES.TREE.LABEL = gsub("[.]LCG", "", SPECIES.TREE.LABEL)) %>%
  mutate(SPECIES.TREE.LABEL = gsub("[_]", " ", SPECIES.TREE.LABEL)) %>%
  mutate(SPECIES.TREE.LABEL = gsub("LCG$", "", SPECIES.TREE.LABEL)) %>%
  separate(SPECIES.TREE.LABEL, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = ploidyfix2stl  %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

ploidyfix2stl = ploidyfix2stl %>%
  mutate(SPECIES.TREE.LABEL = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(SPECIES.TREE.LABEL = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), SPECIES.TREE.LABEL))

kmer_hists = list()
af_hists = list()
i=1
for (ss in snp_stat_files) {
  isolate = gsub("[.]snp[_]stats[.]tsv$", "", ss)
  title = ploidyfix2stl %>%
    filter(ploidy_file_prefix == isolate) %>%
    pull(SPECIES.TREE.LABEL)
  this_stddev = approx_norm_binom_stdev_tbl %>%
    filter(iso == isolate) %>%
    pull(approx_norm_binom_stdev)

  df = read_delim(ss, delim="\t") %>%
    filter(p != 0.00 & p != 1.00)

  afhist = GatkAfHistGen(df, isolate, this_stddev, stats = snp_densities, title = title)

  kmer_pref = paste("/scratch/amsesk/pursuit/ploidy/kmerhist/", isolate, "_23", sep = "")
  kmerhist = KmerHistGen(kmer_pref, title = title)

  af_hists[[i]] = afhist
  kmer_hists[[i]] = kmerhist
  i = i + 1
}

pageid = 1
for (page_idx in seq(1,i,4)) {
  plt_indices = seq(page_idx,page_idx+3,1)
  rm = which(plt_indices > 112)
  if (length(rm) != 0 ) {
    plt_indices = plt_indices[-rm]
  }
  for (p in plt_indices) {
    if (p == min(plt_indices)) {
      page = kmer_hists[[p]] + af_hists[[p]]
    } else {
      page = page + kmer_hists[[p]] + af_hists[[p]]
    }
  }

  page = page + plot_layout(ncol=4)

  page_name = paste("/home/amsesk/dev/Chytrid-Phylogenomics/figures/suppl_hists/page", pageid, ".pdf", sep="")
  ggsave(page_name,
         plot = page,
         height = 8.5,
         width =11,
         units = "in",
         device = cairo_pdf)

  pageid = pageid +1
}


afhist = GatkAfHistGen(df, isolate, this_value, stats = snp_densities)

pdf(file = paste(isolate, ".pdf", sep=""), width = 11, height = 8.5, onefile=FALSE)
if (!is.na(args[3])) {
  kmerhist = KmerHistGen(args[3])
  kmerhist + afhist
  #plt_list = list()
  #plt_list[[1]] = kmerhist
  #plt_list[[2]] = afhist
  #plt_pane(plt_list, 1, 2)
} else {
  afhist
}
dev.off()

#### Plot with Rscript from shell ####
args = commandArgs(trailingOnly = TRUE)
df_path = args[1]
isolate = args[2]
approx_norm_binom_stdev_tbl = args[4]

if (!is.null(approx_norm_binom_stdev_tbl)) {
  approx_norm_binom_stdev_tbl = read_delim(approx_norm_binom_stdev_tbl, delim = "\t", col_names = F)
  approx_norm_binom_stdev_tbl = approx_norm_binom_stdev_tbl %>%
    select(X1,X4) %>%
    rename(iso=X1, approx_norm_binom_stdev=X4)

  this_value = approx_norm_binom_stdev_tbl %>%
    filter(iso == isolate) %>%
    pull(approx_norm_binom_stdev)
} else {
  this_value = NULL
}

print(args)
df = read_delim(args[1], delim="\t") %>%
  filter (p != 0.00 & p != 1.00)

afhist = GatkAfHistGen(df, isolate, this_value, stats = snp_densities)

pdf(file = paste(isolate, ".pdf", sep=""), width = 11, height = 8.5, onefile=FALSE)
if (!is.na(args[3])) {
  kmerhist = KmerHistGen(args[3])
  kmerhist + afhist
  #plt_list = list()
  #plt_list[[1]] = kmerhist
  #plt_list[[2]] = afhist
  #plt_pane(plt_list, 1, 2)
} else {
  afhist
}
dev.off()
#ggsave(plot=plt, filename = paste(isolate, ".pdf", sep=""), device="pdf", width = 22, height = 17, units = "in")
