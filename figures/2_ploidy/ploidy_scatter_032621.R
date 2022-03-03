library(tidyverse)
library(readxl)
library(plotrix)
library(scales)

weighted.sd = function(values, weights) {
  n_values = length(values)
  wm = weighted.mean(values, weights)
  val_min_wm_sq = (values - wm)^2
  the_sum = sum(val_min_wm_sq)
  v = the_sum / n_values
  sd = sqrt(v)
  return(sd)
}

l50_genome_sizes = read_delim("~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv", delim="\t") %>%
  select(SPECIES.TREE.LABEL, l50_assembly_length)
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx")

ploidy_df = read_delim("~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv", delim="\t") %>%
  select(SPECIES.TREE.LABEL, coding) %>%
  rename(ploidy = coding)

snp_counts = read_delim("~/DATA/pursuit/ploidy_final/all.snp_contig_counts_strainified.tsv", delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  filter(ploidy_file_prefix == "Anasp1")

wm = weighted.mean(snp_counts$snp_density, snp_counts$contig_length)
anasp_snp_densities = snp_counts$snp_density

#anasp_snp_densities -

snp_densities = read_delim("~/DATA/pursuit/ploidy_final/all.snp_contig_counts_strainified.tsv", delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(snp_density_mean = weighted.mean(snp_density, contig_length), stdev = weighted.sd(snp_density, contig_length)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  select(SPECIES.TREE.LABEL, snp_density_mean, stdev)

#snp_densities = read_delim("~/DATA/pursuit/ploidy_final/all.snp_contig_counts_strainified.tsv", delim="\t", col_names = F) %>%
#  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
#  group_by(ploidy_file_prefix) %>%
#  summarise(num_snps = sum(num_snps)) %>%
#  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
#  left_join(l50_genome_sizes) %>%
#  select(SPECIES.TREE.LABEL, num_snps, l50_assembly_length) %>%
#  mutate(snp_density = num_snps/l50_assembly_length) %>%
#  group_by(SPECIES.TREE.LABEL) %>%
#  summarise(snp_density_mean = mean(snp_density), stddev = sd(snp_density))

binom_expect = read_delim("~/DATA/pursuit/ploidy_final/all_AF_ranges_from_binom.tsv", delim="\t", col_names=F) %>%
  rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  select(-ploidy_file_prefix) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

combined = binom_expect %>%
  left_join(snp_densities) %>%
  left_join(ploidy_df)

colpal = c(
  "lightgrey",
  "red",
  "blue"
)

ploidy_scatter = ggplot(combined, aes(x = p_in_binom_expect, y = snp_density_mean, fill = ploidy)) +
  geom_point(pch = 21, cex = 2.5) +
  geom_errorbar(aes(ymax = snp_density_mean + stdev, ymin = snp_density_mean - stdev)) +
  scale_fill_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1),
    plot.margin = unit(c(0,0,0,2), "in")
  ) +
  stat_ellipse(data=subset(combined, ploidy %in% c("2", "1")), geom="polygon", aes(fill=ploidy), alpha=0.4, type="norm", color="black")
ploidy_scatter
