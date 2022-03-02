library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

LENGTH_PATH="~/DATA/pursuit/ploidy_final/assembly_l50_lengths/"
get_contig_lengths_from_file = function(asm) {
  d = read_delim(file.path(LENGTH_PATH, paste(asm, ".l50.contig_lengths", sep="")), delim="\t", col_names=F) %>%
    rename(contig = X1, length = X2)
  d
}

#setwd("/home/aimzez/DATA/pursuit/ploidy")
#args=c("contig.depth/Basidiobolus_heterosporus_B8920.N168.v1.sorted.bam.depth.ContigMean", "snps_stats/Basidiobolus_heterosporus_B8920.N168.v1.snp_stats.tsv")

### Read in depth by contig from a BAM file
contig_depths = read_delim(args[1], delim="\t", col_names=T)

ploidy_data_table = read_delim("~/DATA/pursuit/ploidy_final/ploidy_data_table.csv", delim=",", col_names=F) %>%
  rename(Strain =X1, Asm =X2) %>%
  select(-X3, -X4)

strain = gsub("[.]sorted[.]bam[.]depth[.]ContigMean[.]L50", "", basename(args[1]))
strain_asm = ploidy_data_table %>%
  filter (Strain == strain) %>%
  .$Asm
strain_asm = basename(strain_asm)

stopifnot(length(strain_asm) == 1)

contig_lengths = get_contig_lengths_from_file(strain_asm[1])

### Set missing contig depths to 0
contig_lengths_and_depths = contig_lengths %>%
  left_join(contig_depths) %>%
  mutate(depth = ifelse(is.na(depth), 0.0, depth))

mean_coverage = round(weighted.mean(contig_lengths_and_depths$depth, contig_lengths_and_depths$length))

#### Generate binomial distributions to determine x-axis cutoffs ####

#min = 0
#max = 10000
#coverage = mean_coverage
#step = 1

#x = seq(min, max, step)

#y = dbinom(x, coverage, 0.5)

#df = as_tibble(cbind(x,y)) %>%
#  mutate(x = (x*(max/coverage))/max)

#stdev = (sqrt(coverage*0.5*(1-0.5))*(max/coverage))/max

#snp_fit_range = c(0.5-stdev, 0.5+stdev)

#### Approximate binomial distributions with normal distributions ONLY when coverage*0.5 >= 5, or coverage >= 10 ####
min = 0
max = 1000000
coverage = mean_coverage
step = 1

# For a normal distrubtion approximation of a binomial distribution...

# mean = n*p
approx_mean = coverage*0.5

# stdev = sqrt(n*p*q)
approx_stdev = sqrt(coverage*0.5*0.5)

x = seq(min, max, step)
y = dnorm(x = x, mean = approx_mean, sd = approx_stdev)

df = as_tibble(cbind(x,y)) %>%
  mutate(x = (x*(max/coverage))/max)

approx_stdev_scaled = approx_stdev*((max/coverage)/max)

snp_fit_range = c(0.5-approx_stdev_scaled, 0.5+approx_stdev_scaled)

expected_snps_in_range = 0.683

### Prepare outputs to another script
snp_stats = read_delim(args[2], delim="\t", col_names=T)

total_snps = dim(snp_stats)[1]

snps_in_range = snp_stats %>%
  filter(p >= snp_fit_range[1] & p <= snp_fit_range[2])

total_snps_in_range = dim(snps_in_range)[1]

### For binomial

cat(paste(strain), mean_coverage, total_snps_in_range/total_snps, approx_stdev_scaled,  sep = "\t")
cat("\n")


### For binomials APPROXIMATED with normal distributions
#if (mean_coverage > 10) {
#  cat(paste(strain), mean_coverage, total_snps_in_range/total_snps, expected_snps_in_range, sep = "\t")
#  cat("\n")
#}
