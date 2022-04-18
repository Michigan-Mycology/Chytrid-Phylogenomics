library(tidyverse)
library(scales)
library(readxl)
library(tublerone)
library(argparse)
library(readxl)

#### Plot with Rscript from shell ####
parser = argparse::ArgumentParser(description = "Plot an allele frequency histogram.")
parser$add_argument('-s', '--snp_stats', action = "store", help = "Path to *.snp_stats.tsv file output from `vcf_to_af.py`: A 5-column (with column names) tab-separated text file with contig, p, q, depth, and mqrs.", required = TRUE)
parser$add_argument('-c', '--snp_contig_counts', action = "store", help = "Path to strainified output from `snp_contig_counts.py`: A 5-column (no column names) tab-separated text file with isolate name, contig, number of SNPs, SNP density, and contig length.")
parser$add_argument('-t', '--isolates_sheet', action = "store", help = "Path to isolate metadata sheet (in xlsx format) for ploidy_file_prefix and actual isolate name resolution. Must have `SPECIES.TREE.LABEL` and `ploidy_file_prefix` columns.", required = TRUE)
parser$add_argument('-r', '--traits_sheet', action = "store", help = "Path to isolate traits sheet (in tsv format) for accessing l50 assembly length. Must have `SPECIES.TREE.LABEL` and `l50_assembly_length`` columns.")
parser$add_argument('-i', '--isolate', action = "store", help = "Name of the isolate.", required = TRUE)
parser$add_argument('-e', '--expected_dist', action = "store", help = "Path to output from `binomial_dists.R`: A 4 column (no column names) tab-separated text file with isolate name, mean depth, proportion of SNPs in the expected region of the normal-approximated null binomial distribution, and the approximate standard deviation of that distribution.")
parser$add_argument('-k', '--kmerhist', action = "store", help = "Prefix of *.peaks and *.khist files from kmercountexact. Do not pass to exclude.")
parser$add_argument('--kmer_xmax', action = "store", help = "Upper limit of x-axis for plots where auto-limited didn't work very well.")
parser$add_argument('--kmer_ymax', action = "store", help = "Upper limit of y-axis for plots where auto-limited didn't work very well.")
args = parser$parse_args()

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

  expected_dist = read_delim(args$expected_dist, delim="\t", col_names = F) %>%
    rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3, approx_stdev = X4) %>%
    left_join(isolates_sheet %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
    filter(ploidy_file_prefix == isolate) %>%
    select(-ploidy_file_prefix) %>%
    select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect, approx_stdev)
  approx_stdev = expected_dist %>%
    pull(approx_stdev)

  l50_genome_sizes = read_delim(args$traits_sheet, delim = "\t") %>%
    select(SPECIES.TREE.LABEL, l50_assembly_length)

  extra_snp_stats = extra_snp_stats %>%
    left_join(expected_dist) %>%
    left_join(l50_genome_sizes) %>%
    mutate(snp_density = num_snps/l50_assembly_length)

  afhist = GatkAfHistGen(snps, isolate, approx_stdev, stats = extra_snp_stats)
}  else if (!is.null(args$expected_dist) || !is.null(args$traits_sheet)) {
  stop("-t|--traits_sheet and -e|--expected_dist must either be set together or not at all")
} else {
  afhist = GatkAfHistGen(snps, isolate)
}

pdf(file = paste(isolate, ".pdf", sep=""), width = 11, height = 8.5, onefile=FALSE)
if (!is.null(args$kmerhist)) {
  library(patchwork)
  kmerhist = KmerHistGen(args$kmerhist, xmax = args$kmer_xmax, ymax = args$kmer_ymax)
  kmerhist + afhist
} else {
  afhist
}
dev.off()
