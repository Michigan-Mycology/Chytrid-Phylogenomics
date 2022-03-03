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

binom_expect = read_delim("~/DATA/pursuit/ploidy/all_AF_ranges_from_binom.tsv", delim="\t", col_names=F)

binom_expect_isolates_joined_by_SPECIES.TREE.LABEL = binom_expect %>%
  rename(SPECIES.TREE.LABEL = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

binom_expect_isolates_joined_by_STRAIN = binom_expect %>%
  rename(STRAIN = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates) %>%
  select(-ploidy_file_prefix, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-STRAIN, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

binom_expect_isolates_joined_by_ploidy_file_prefix = binom_expect %>%
  rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates) %>%
  select(-STRAIN, -`ploidy (12-11-20)`) %>%
  drop_na %>%
  select(-ploidy_file_prefix, -LTP) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

merged_binom_expect = rbind(binom_expect_isolates_joined_by_SPECIES.TREE.LABEL, binom_expect_isolates_joined_by_STRAIN, binom_expect_isolates_joined_by_ploidy_file_prefix)


#ploidy_strains = list.files("~/DATA/pursuit/ploidy/snp_contig_counts/strainified/") %>%
#  as_tibble %>%
#  mutate(value = gsub("[.]counts[.]strainified[.]tsv", "", value)) %>%
#  rename(ploidy_file_prefix = value)


merged_counts_summary = merged_counts %>%
  group_by(SPECIES.TREE.LABEL) %>%
  summarise(snp_density_mean = mean(snp_density), stderror = std.error(snp_density)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, "ploidy (12-11-20)")) %>%
  rename(inferred_ploidy = `ploidy (12-11-20)`) %>%
  mutate(inferred_ploidy = ifelse(inferred_ploidy == "FOUND" | inferred_ploidy == "FIND" | is.na(inferred_ploidy), "Not_Yet_Inferred", inferred_ploidy))
merged_counts_summary = merged_counts_summary %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary$SPECIES.TREE.LABEL[order(merged_counts_summary$snp_density_mean, decreasing = TRUE)])) %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "3N", "Not_Yet_Inferred")))

merged_counts_summary$inferred_ploidy
colpal = c(
  "blue",
  "red",
  muted("blue"),
  muted("red"),
  "yellow",
  "darkgreen",
  "grey"
)
### Make a point plot

ggplot(merged_counts_summary, aes(x = SPECIES.TREE.LABEL, y = snp_density_mean, fill = inferred_ploidy)) +
  geom_point(pch = 21, cex = 2) +
  geom_errorbar(aes(ymax = snp_density_mean + stderror, ymin = snp_density_mean - stderror)) +
  scale_fill_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1, size = 5),
    plot.margin = unit(c(0,0,0,2), "in")
  )

# Make a boxplot
merged_counts = merged_counts %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary$SPECIES.TREE.LABEL[order(merged_counts_summary$snp_density_mean, decreasing = TRUE)])) %>%
  left_join(merged_counts_summary, by="SPECIES.TREE.LABEL")
#box_color_map = merged_counts_summary %>% select(SPECIES.TREE.LABEL, inferred_ploidy)
ggplot(merged_counts, aes(x = SPECIES.TREE.LABEL, y = snp_density, fill=inferred_ploidy, color=inferred_ploidy)) +
  geom_boxplot() +
  geom_point(data=merged_counts_summary, aes(x=SPECIES.TREE.LABEL, y=snp_density_mean), fill="black", cex=2, pch=23) +
  scale_fill_manual(values = colpal) +
  scale_color_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1, size = 5),
    plot.margin = unit(c(0,0,0,2), "in")
  )

#Binom expect proportion by mean snp density
merged_counts_summary2 = merged_counts %>%
  group_by(SPECIES.TREE.LABEL) %>%
  summarise(snp_density_mean = mean(snp_density), snp_density_stdev = sd(snp_density)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, "ploidy (12-11-20)")) %>%
  rename(inferred_ploidy = `ploidy (12-11-20)`) %>%
  mutate(inferred_ploidy = ifelse(inferred_ploidy == "FOUND" | inferred_ploidy == "FIND" | is.na(inferred_ploidy), "Not_Yet_Inferred", inferred_ploidy))
merged_counts_summary2 = merged_counts_summary2 %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary2$SPECIES.TREE.LABEL[order(merged_counts_summary2$snp_density_mean, decreasing = TRUE)])) %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "3N", "Not_Yet_Inferred"))) %>%
  left_join(merged_binom_expect)

write.table(x=merged_counts_summary2, file="~/DATA/pursuit/ploidy/ploidy_scatter_data.tsv", sep="\t", quote=F, col.names = T, row.names = F)
ploidy_scatter = ggplot(merged_counts_summary2, aes(x = p_in_binom_expect, y = snp_density_mean, fill = inferred_ploidy)) +
  geom_point(pch = 21, cex = 2.5) +
  geom_errorbar(aes(ymax = snp_density_mean + snp_density_stdev, ymin = snp_density_mean - snp_density_stdev)) +
  scale_fill_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1),
    plot.margin = unit(c(0,0,0,2), "in")
  )

#Mean contig length by mean snp density
merged_counts_summary3 = merged_counts %>%
  group_by(SPECIES.TREE.LABEL) %>%
  summarise(snp_density_mean = mean(snp_density), snp_density_stderror = std.error(snp_density),
            contig_length_mean = mean(contig_length), contig_length_stderror = std.error(contig_length)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, "ploidy (12-11-20)")) %>%
  rename(inferred_ploidy = `ploidy (12-11-20)`) %>%
  mutate(inferred_ploidy = ifelse(inferred_ploidy == "FOUND" | inferred_ploidy == "FIND" | is.na(inferred_ploidy), "Not_Yet_Inferred", inferred_ploidy))
merged_counts_summary2 = merged_counts_summary2 %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=merged_counts_summary$SPECIES.TREE.LABEL[order(merged_counts_summary$snp_density_mean, decreasing = TRUE)])) %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "Not_Yet_Inferred")))


ggplot(merged_counts_summary2, aes(x = log10(contig_length_mean), y = snp_density_mean, fill = inferred_ploidy)) +
  geom_point(pch = 21, cex = 5) +
  geom_errorbar(aes(ymax = snp_density_mean + snp_density_stderror, ymin = snp_density_mean - snp_density_stderror)) +
  geom_errorbar(aes(xmax = log10(contig_length_mean + contig_length_stderror), xmin = log10(contig_length_mean - contig_length_stderror))) +
  scale_fill_manual(values = colpal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1),
    plot.margin = unit(c(0,0,0,2), "in")
  )

### Maybe haploid assemblies to run BLAT for
haps = merged_counts_summary %>%
  filter(inferred_ploidy == "1N" | inferred_ploidy == "1N?") %>%
  left_join(read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>% select(SPECIES.TREE.LABEL, `teamchytrid:Assembly/Final/genomes`))
View(haps)
