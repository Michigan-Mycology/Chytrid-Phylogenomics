library(phytools)
library(tidyverse)
library(readxl)
library(ggtree)

tree = read.newick("~/DATA/pursuit/phylogeny_final/concat/treefiles/unpart_concat_NPB_from-consensus.treefile")
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet=1) %>%
  select (LTP, ASTRAL_TaxID)

new = c()
for (t in tree$tip.label) {
  new = c(new, isolates %>% filter(LTP == t) %>% .$ASTRAL_TaxID)
}

tree$tip.label = unlist(lapply(X = tree$tip.label, FUN = gsub, pattern = "[&].*$", replacement = ""))

tree$tip.label = new
write.tree(tree, file = "~/DATA/pursuit/phylogeny_final/concat/treefiles/unpart_concat_NFB.treefile.renamed_with_ASTRAL_names")
