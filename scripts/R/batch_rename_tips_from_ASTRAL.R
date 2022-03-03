library(phytools)
library(tidyverse)
library(readxl)

nametips = function(lab, isolates) {
  new = isolates %>% 
    filter(ASTRAL_TaxID == lab) %>%
    select (SPECIES.TREE.LABEL) %>%
    .$SPECIES.TREE.LABEL
  new
}

args = commandArgs(trailingOnly = TRUE)
tree_path = args[1]
isolates_path = args[2]

setwd(tree_path)
isolates = read_xlsx(isolates_path) %>%
  select(ASTRAL_TaxID, SPECIES.TREE.LABEL)

files = list.files(path = tree_path)
#markers = lapply(files, (function(f) strsplit(f, split="[.]")[[1]][1]) )

flist = as.list(files)
#names(flist) = unlist(markers)

for (f in flist) {
  t = read.newick(f)
  t$tip.label = sapply(as.list(t$tip.label), nametips, isolates = isolates)
  write.tree(t, file = paste(f, ".fullnames"))
}