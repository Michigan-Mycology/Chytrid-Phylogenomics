library(tidyverse)
library(phytools)

setwd("/home/amsesk/data/pursuit/more_hmm/final_gene_trees/")

concat_tree = read.newick(file.path(CHYTRID_PHYLO, "figures/1_tree", "combined_tree_filtered_support_RENAMED.tre"))

cellwall_occ = as_tibble(concat_tree$tip.label) %>%
  rename(SPECIES.TREE.LABEL = value)
for (f in list.files(getwd())) {
  if (endsWith(f, ".tre.renamed")) {
    marker = gsub("[.]tre[.]renamed", "", f)
    tree = read.newick(f)
    these_tips = as_tibble(tree$tip.label) %>%
      separate(col = value, into=c("SPECIES.TREE.LABEL"), sep = "[&]") %>%
      mutate(!!marker := 1)
    cellwall_occ = cellwall_occ %>%
      left_join(these_tips)
  }
}

cellwall_occ = cellwall_occ %>%
  mutate_at(c("AGS", "DHCR7", "ERG11", "ERG5", "FKS1", "MYSc_Myo17"),
            .funs = list(~ ifelse(is.na(.), 0, 1)) ) %>%
  gather(key = "gene", value = "present", -SPECIES.TREE.LABEL) %>%
  rename(ID = SPECIES.TREE.LABEL) %>%
  mutate(present = factor(present, levels=c(0,1)))
