library(tidyverse)
library(ggtree)
library(readxl)
library(phytools)

astralify_tips = function(phylo, isolates) {
  new_tips = c()
  for (tip in phylo$tip.label) {
    ltp = stringr::str_split(tip, "[|]")[[1]][1]
    new = isolates %>%
      filter(LTP == ltp) %>%
      .$ASTRAL_TaxID

    new_tips = c(new_tips, new)
  }
  phylo$tip.label = new_tips
  return (phylo)
}

SUFFIX = "contree"

args = commandArgs(trailingOnly = T)
fs = list.files(args[1])

isolates = read_xlsx("~/work/neozygites/Neozygites_Isolates.xlsx")

for (f in fs) {
  if (!endsWith(f, SUFFIX)) {
    next
  }
  f = file.path(args[1], f)
  tree = read.newick(f)
  renamed_tree = astralify_tips(tree, isolates)
  write.tree(renamed_tree,
             file = gsub(SUFFIX, paste("ASTRAL.renamed", SUFFIX, sep = "."), f))
}
