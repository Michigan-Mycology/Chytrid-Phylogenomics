library(readxl)
library(tidyverse)

# 

sheet = read_xlsx("~/work/pursuit/sheets/non_chytrids.xlsx") %>%
  select(Got_it, SPECIES.TREE.LABEL, Assembly, SRX) %>%
  filter(Got_it == 1) %>%
  filter(!is.na(SRX))

write_delim(sheet, col_names = FALSE, delim="\t", path = "non_chytrds_srx.tsv", quote_escape=FALSE)


asmlen = read_delim("~/DATA/pursuit/chytrids_asm_length.tsv", delim="\t", col_names = FALSE)
colnames(asmlen) = c("assembly", "assembly_length")

sheet = read_delim("~/DATA/pursuit/ploidy/data_table_revised.csv", delim=",", col_names=FALSE) %>%
  select(X1, X2) %>%
  mutate(X2 = basename(X2))
colnames(sheet) = c("strain", "assembly")

snpraw = read_delim("~/DATA/pursuit/chytrid_snp_raw.tsv", delim="\t", col_names=FALSE)
colnames(snpraw) = c("strain", "snps") 

combined = sheet %>%
  left_join(asmlen) %>%
  left_join(snpraw) %>%
  mutate (snp_density = snps/assembly_length)

write_delim(combined, delim=',', file = "~/DATA/pursuit/chytrid_snp_data_combined.csv", quote=FALSE)

#######################

newdata = read_delim('~/DATA/pursuit/new_data/vcf/newdata_snp_data_combined.csv', delim=',', col_names=TRUE)
chytrids = read_delim('~/DATA/pursuit/ploidy/vcf/chytrids_snp_data_combined.csv', delim=',', col_names=TRUE)


combined = bind_rows(newdata, chytrids)
combined = combined %>%
  mutate(strain = factor(strain, levels=combined$strain[order(combined$snp_density, decreasing = TRUE)]))

ggplot(combined, aes(strain, snp_density)) +
  geom_text(aes(label=snps), angle=65, nudge_y = 0.0005) +
  geom_bar(stat="identity") +
  scale_x_discrete(expand=c(0.025,0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust=1, angle=45)
  )
