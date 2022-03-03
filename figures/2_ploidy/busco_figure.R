library(readxl)
library(tidyverse)

genome_sizes = read_delim("~/work/pursuit/sheets/Phylogenomics_Ultrastructure_renamed.tsv", delim="\t") %>%
  select(SPECIES.TREE.LABEL, assembly_length)
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx")

snp_densities = read_delim("~/DATA/pursuit/ploidy_final/all.snp_contig_counts_strainified.tsv", delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(num_snps = sum(num_snps)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  left_join(genome_sizes) %>%
  select(SPECIES.TREE.LABEL, num_snps, assembly_length) %>%
  mutate(snp_density = num_snps/assembly_length)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, ploidy_file_prefix, `teamchytrid:Assembly/Final/genomes`) %>%
  rename(assembly = `teamchytrid:Assembly/Final/genomes`) %>%
  mutate(assembly = gsub("[.]gz", "", assembly))

busco_sheet = read_xlsx("~/work/pursuit/sheets/genome_buscos.xlsx") %>%
  left_join(isolates) %>%
  left_join(snp_densities) %>%
  mutate(perc_dup = complete_and_duplicated_buscos/complete_buscos)
busco_sheet = busco_sheet %>%
  mutate(SPECIES.TREE.LABEL = factor(SPECIES.TREE.LABEL, levels=busco_sheet$SPECIES.TREE.LABEL[order(busco_sheet$perc_dup)])) %>%
  drop_na
busco_sheet

ggplot(busco_sheet, aes(x=SPECIES.TREE.LABEL, y=perc_dup)) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=6)
  )

