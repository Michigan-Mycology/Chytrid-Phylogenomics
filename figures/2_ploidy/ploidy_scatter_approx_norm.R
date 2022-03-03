library(tidyverse)
library(readxl)
library(plotrix)
library(scales)


### Merge all_snp_contig_counts with Pursuit_Isolates to standardize names between file names and Pursuit_Isolates
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP, STRAIN, ploidy_file_prefix, "ploidy (12-11-20)")

counts = read_delim("~/DATA/pursuit/ploidy/all_snp_contig_counts.strainified.tsv", delim="\t", col_names = FALSE)

counts_isolates_joined_by_SPECIES.TREE.LABEL = counts %>%
  rename(SPECIES.TREE.LABEL = X1, contig = X2, n_snps = X3, snp_density = X4, contig_length = X5) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, contig, n_snps, snp_density, contig_length)

counts_isolates_joined_by_STRAIN = counts %>%
  rename(STRAIN = X1, contig = X2, n_snps = X3, snp_density = X4, contig_length = X5) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, contig, n_snps, snp_density, contig_length)

counts_isolates_joined_by_ploidy_file_prefix = counts %>%
  rename(ploidy_file_prefix = X1, contig = X2, n_snps = X3, snp_density = X4, contig_length = X5) %>%
  left_join(isolates) %>%
  select(-STRAIN, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-ploidy_file_prefix, -LTP) %>%
  select(SPECIES.TREE.LABEL, contig, n_snps, snp_density, contig_length)

merged_counts = rbind(counts_isolates_joined_by_SPECIES.TREE.LABEL, counts_isolates_joined_by_STRAIN, counts_isolates_joined_by_ploidy_file_prefix)

######

approx_norm_expect = read_delim("~/DATA/pursuit/ploidy/all_AF_ranges_from_binom_APPROX_norm.tsv", delim="\t", col_names=F)

approx_norm_expect_isolates_joined_by_SPECIES.TREE.LABEL = approx_norm_expect %>%
  rename(SPECIES.TREE.LABEL = X1, mean_coverage = X2, observed = X3, expected = X4) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, observed, expected)

approx_norm_expect_isolates_joined_by_STRAIN = approx_norm_expect %>%
  rename(STRAIN = X1, mean_coverage = X2, observed = X3, expected = X4) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, observed, expected)

approx_norm_expect_isolates_joined_by_ploidy_file_prefix = approx_norm_expect %>%
  rename(ploidy_file_prefix = X1, mean_coverage = X2, observed = X3, expected = X4) %>%
  left_join(isolates) %>%
  select(-STRAIN, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-ploidy_file_prefix, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, observed, expected)

merged_approx_norm_expect = rbind(approx_norm_expect_isolates_joined_by_SPECIES.TREE.LABEL, approx_norm_expect_isolates_joined_by_STRAIN, approx_norm_expect_isolates_joined_by_ploidy_file_prefix)

merged_counts_summary = merged_counts %>%
  group_by(SPECIES.TREE.LABEL) %>%
  summarise(snp_density_mean = mean(snp_density), stderror = std.error(snp_density)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, "ploidy (12-11-20)")) %>%
  rename(inferred_ploidy = `ploidy (12-11-20)`) %>%
  mutate(inferred_ploidy = ifelse(inferred_ploidy == "FOUND" | inferred_ploidy == "FIND" | is.na(inferred_ploidy), "Not_Yet_Inferred", inferred_ploidy))
merged_counts_summary = merged_counts_summary %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary$SPECIES.TREE.LABEL[order(merged_counts_summary$snp_density_mean, decreasing = TRUE)])) %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "Not_Yet_Inferred")))

merged_counts_summary$inferred_ploidy
colpal = c(
  "blue",
  "red",
  muted("blue"),
  muted("red"),
  "yellow",
  "grey"
)

# Observed-Expected with normal-approximated binomial distributions
merged_counts_summary2 = merged_counts %>%
  group_by(SPECIES.TREE.LABEL) %>%
  summarise(snp_density_mean = mean(snp_density), snp_density_stdev = sd(snp_density)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, "ploidy (12-11-20)")) %>%
  rename(inferred_ploidy = `ploidy (12-11-20)`) %>%
  mutate(inferred_ploidy = ifelse(inferred_ploidy == "FOUND" | inferred_ploidy == "FIND" | is.na(inferred_ploidy), "Not_Yet_Inferred", inferred_ploidy))
merged_counts_summary2 = merged_counts_summary2 %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary$SPECIES.TREE.LABEL[order(merged_counts_summary$snp_density_mean, decreasing = TRUE)])) %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "Not_Yet_Inferred"))) %>%
  left_join(merged_approx_norm_expect)


ggplot(merged_counts_summary2, aes(x = abs(observed-expected), y = snp_density_mean, fill = inferred_ploidy)) +
  geom_point(pch = 21, cex = 5) +
  geom_errorbar(aes(ymax = snp_density_mean + snp_density_stdev, ymin = snp_density_mean - snp_density_stdev)) +
  scale_fill_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1),
    plot.margin = unit(c(0,0,0,2), "in")
  )
