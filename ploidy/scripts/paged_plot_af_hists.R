library(tidyverse)
library(scales)
library(readxl)
library(tublerone)
library(argparse)
library(readxl)
library(patchwork)

produce_single_plot = function(args) {
  isolate = args$isolate
  isolates_sheet = read_xlsx(args$isolates_sheet)
  snps = read_delim(args$snp_stats, delim = "\t") %>%
    filter(p != 0.00 & p != 1.00)
  extra_snp_stats = read_delim(args$snp_contig_counts, delim = "\t", col_names = F) %>%
    rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
    group_by(ploidy_file_prefix) %>%
    summarise(num_snps = sum(num_snps)) %>%
    left_join(isolates_sheet %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix))

  if (!is.null(args$expected_dist) && !is.null(args$traits_sheet)) {

    print(args$isolate)
    expected_dist = read_delim(args$expected_dist, delim="\t", col_names = F) %>%
      rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3, approx_stdev = X4) %>%
      left_join(isolates_sheet %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
      filter(ploidy_file_prefix == isolate) %>%
      select(-ploidy_file_prefix) %>%
      select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect, approx_stdev)
    approx_stdev = expected_dist %>%
      pull(approx_stdev)
    print(expected_dist)

    l50_genome_sizes = read_delim(args$traits_sheet, delim = "\t") %>%
      select(SPECIES.TREE.LABEL, l50_assembly_length, coding)

    extra_snp_stats = extra_snp_stats %>%
      left_join(expected_dist) %>%
      left_join(l50_genome_sizes) %>%
      mutate(snp_density = num_snps/l50_assembly_length)

    title = expected_dist %>% pull(SPECIES.TREE.LABEL)
    title = gsub("[_]", " ", title)
    afhist = GatkAfHistGen(snps, isolate, approx_stdev, stats = extra_snp_stats, title = "")
  }  else if (!is.null(args$expected_dist) || !is.null(args$traits_sheet)) {
    stop("-t|--traits_sheet and -e|--expected_dist must either be set together or not at all")
  } else {
    afhist = GatkAfHistGen(snps, isolate)
  }

  if (!is.null(args$kmerhist)) {
    library(patchwork)
    kmerhist = KmerHistGen(args$kmerhist, xmax = args$kmer_xmax, ymax = args$kmer_ymax, title = title)
    kmerhist + afhist
  } else {
    afhist
  }
}

#### Plot with Rscript from shell ####
parser = argparse::ArgumentParser(description = "Plot an allele frequency histogram.")
parser$add_argument('-s', '--snp_stats_dir', action = "store", help = "Path to *.snp_stats.tsv file output from `vcf_to_af.py`: A 5-column (with column names) tab-separated text file with contig, p, q, depth, and mqrs.", required = TRUE)
parser$add_argument('-c', '--snp_contig_counts', action = "store", help = "Path to strainified output from `snp_contig_counts.py`: A 5-column (no column names) tab-separated text file with isolate name, contig, number of SNPs, SNP density, and contig length.")
parser$add_argument('-t', '--isolates_sheet', action = "store", help = "Path to isolate metadata sheet (in xlsx format) for ploidy_file_prefix and actual isolate name resolution. Must have `SPECIES.TREE.LABEL` and `ploidy_file_prefix` columns.", required = TRUE)
parser$add_argument('-r', '--traits_sheet', action = "store", help = "Path to isolate traits sheet (in tsv format) for accessing l50 assembly length. Must have `SPECIES.TREE.LABEL` and `l50_assembly_length`` columns.")
parser$add_argument('-e', '--expected_dist', action = "store", help = "Path to output from `binomial_dists.R`: A 4 column (no column names) tab-separated text file with isolate name, mean depth, proportion of SNPs in the expected region of the normal-approximated null binomial distribution, and the approximate standard deviation of that distribution.")
parser$add_argument('-k', '--kmerhist_dir', action = "store", help = "Prefix of *.peaks and *.khist files from kmercountexact. Do not pass to exclude.")
parser$add_argument('--kmer_maxes_tsv', action = "store", help = "3-column tab-separated sheet of isolate name, xmax, and ymax for plots where auto-limited didn't work very well.")
parser$add_argument('--plots_per_page', action = "store", help = "Number of plots to print per page.", default = "4")
args = parser$parse_args()

plts = list()
i=1
for (f in list.files(args$snp_stats_dir)) {
  args = parser$parse_args()

  if (!is.null(args$kmer_maxes_tsv)) {
    kmer_maxes = read_delim(args$kmer_maxes_tsv, col_names = F, delim = "\t")
  }

  args$kmer_xmax = NULL
  args$kmer_ymax = NULL

  args$isolate = gsub("[.]snp[_]stats[.]tsv", "", f)
  args$snp_stats = file.path(args$snp_stats_dir, f)
  args$kmerhist = file.path(args$kmerhist_dir, paste(args$isolate, "_23", sep = ""))

  if (!is.null(args$kmer_maxes_tsv)) {
    these_kmer_maxes = kmer_maxes %>%
      filter(X1 == args$isolate)
    if (dim(these_kmer_maxes)[1] == 1) {
      args$kmer_xmax = these_kmer_maxes %>% pull(X2)
      args$kmer_ymax = these_kmer_maxes %>% pull(X3)
      if (is.na(args$kmer_xmax)) {
        args$kmer_xmax = NULL
      }
      if (is.na(args$kmer_ymax)) {
        args$kmer_ymax = NULL
      }
    }
  }
  plts[[i]] = produce_single_plot(args)

  i=i+1
}

plots_per_page = as.numeric(args$plots_per_page)
page_number = 1
for (i in seq(1, length(plts), plots_per_page)) {

  d = "1233
  4455"
  page = plts[[i]] + plts[[i+1]] + plts[[i+2]] + plts[[i+3]] + plot_layout(design=d)
  ggsave(filename = paste("AF_hists_", page_number, ".pdf", sep=""), plot=page, width = 8.5, height = 11, onefile=FALSE)

  page_number = page_number + 1
}

stop()
setwd("/home/aimzez/DATA/pnas_rev/ploidy")

args = NULL
args$expected_dist = "all_AF_ranges_from_binom_withApproxStDev.tsv"
args$isolates_sheet = "../Pursuit_Isolates.xlsx"
args$kmer_maxes_tsv = "manual_maxes_for_plots_nolookit.tsv"
args$kmerhist_dir = "kmerhist"
args$plots_per_page = "4"
args$snp_contig_counts = "all.snp_contig_counts_strainified.tsv"
args$snp_stats_dir = "snp_stats"
args$traits_sheet = "../Pursuit_Phylo_Traits.tsv"
args$isolate = "Backusella_circina_FSU_941.Bacci1.v1"
args$snp_stats = "snp_stats/Backusella_circina_FSU_941.Bacci1.v1.snp_stats.tsv"
args$kmerhist = "kmerhist/Backusella_circina_FSU_941.Bacci1.v1_23"

plts=list()
produce_single_plot(args)

#pdf(file = paste("weird", ".pdf", sep=""), width = 11, height = 8.5, onefile=TRUE)
#produce_single_plot(args)
#plts[[1]]
#dev.off()
