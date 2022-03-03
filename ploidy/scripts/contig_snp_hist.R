library(tidyverse)
library(scales)

args = commandArgs(trailingOnly = TRUE)

#args = c()
#args[1] = "~/DATA/pursuit/ploidy/contig_counts/ARG085.counts.tsv"


counts = read_delim(args[1], delim="\t", col_names = F) %>%
  rename(contig = X1, snp_count = X2, snp_density = X3, contig_length = X4)

total_asm_length = sum(counts$contig_length)

counts = counts %>%
  mutate(contig_length_pTotal = contig_length / total_asm_length) %>%
  mutate(entries_into_hist = contig_length_pTotal*10000) %>%
  group_by(snp_density) %>%
  summarise(entries_into_hist = sum(entries_into_hist)) %>%
  mutate(entries_into_hist = round(entries_into_hist, 0))

hist = c()
for (i in seq(1, dim(counts)[1], 1)) {
  hist = c(hist, rep(counts[i,1] %>% .$snp_density, counts[i,2] %>% .$entries_into_hist))
}
hist = tibble(hist) %>%
  rename(snp_density = hist)

strain = gsub("[.]counts[.]tsv", "", args[1])

density_range = range(hist)
max_density = density_range[2]
min_density = density_range[1]
max_min_diff = max_density - min_density
bin_width_for_plot = max_min_diff/50

plt = ggplot(hist, aes(x=snp_density)) +
  geom_histogram(binwidth = bin_width_for_plot, fill = muted("red"), alpha=0.5, color="black") +
  geom_density(aes(y = bin_width_for_plot * ..count..), fill=muted("red"), alpha=0.7) +
  ggtitle(strain) +

  #Got these numbers from looking at range of snp density across all contigs from all assemblies - to get single, shared x-axis
  #xlim(0,0.06) +

  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggsave(filename = paste(strain, "contig_snp_density", "pdf", sep="."), plot = plt, width = 8.5, height = 11, device = "pdf")
