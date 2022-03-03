library(phytools)
library(tidyverse)
library(readxl)

name_astral_tips = function(lab, isolates) {
  new = isolates %>%
    filter(lab == ASTRAL_TaxID) %>%
    select (LTP) %>%
    .$LTP
  print(new)
  new
}

args = commandArgs(trailingOnly = TRUE)
tree_path = args[1]
isolates_path = args[2]

tree_path="~/DATA/neozygites/phylogeny/ASTRAL/iqtree_NPB_astral.tre"
isolates_path="~/work/neozygites/Neozygites_Isolates.xlsx"

isolates = read_xlsx(isolates_path) %>%
  select(ASTRAL_TaxID, LTP, SPECIES.TREE.LABEL)

tree = read.newick(tree_path)
tree$tip.label = unlist(lapply(tree$tip.label, name_astral_tips, isolates = isolates))

tree$edge.length <- compute.brlen( tree, 1 )$edge.length
#tree$edge.length = NULL

ggtree(tree) + geom_tiplab() + geom_nodelab()
write.tree(tree, file = gsub("[.]tre", ".tre.renamed", tree_path))
write.table(tree$tip.label, file = paste(tree_path, '.tiplabels', sep=''), quote=FALSE, sep="\n", row.names = FALSE, col.names = FALSE)
