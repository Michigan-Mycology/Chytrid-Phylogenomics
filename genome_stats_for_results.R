library(tidyverse)
library(readxl)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx")

genome_stats = read_delim("~/work/Chytrid-Phylogenomics/figures/Supplmental/Supplementary_Table_1_genomes_summary.tsv", delim="\t")

genome_stats = genome_stats %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, `this project`)) %>%
  filter(`this project` == "Y")

summary(genome_stats$`Assembly Length`)
