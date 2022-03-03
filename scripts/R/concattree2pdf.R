library(tidyverse)
library(phytools)
library(ggtree)
library(tublerone)
library(readxl)

treepath = "~/DATA/pursuit/phylogeny_final/concat/treefiles/unpart_concat_NPB_from-consensus.treefile"

tree = read.newick(treepath)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet=1)

old = as.list(isolates$LTP)
names(old) = isolates$LTP
old_to_new = lapply(old, FUN = function(x) isolates %>% select(SPECIES.TREE.LABEL, LTP) %>% filter (LTP == x) %>% .$SPECIES.TREE.LABEL)

tree = tublerone::phylo_rename_tips(tree, old_to_new)
ggtree(tree) +
  geom_tiplab()
