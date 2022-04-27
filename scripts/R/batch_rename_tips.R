library(phytools)
library(tidyverse)
library(readxl)

nametips = function(lab, isolates) {
  spl = strsplit(lab, split="[|]")
  new = isolates %>%
    filter(LTP == spl[[1]][1]) %>%
    select (SPECIES.TREE.LABEL) %>%
    .$SPECIES.TREE.LABEL
  paste(new, lab, sep="&")
}

args = commandArgs(trailingOnly = TRUE)
tree_path = args[1]
isolates_path = args[2]

setwd(tree_path)
isolates = read_xlsx(isolates_path) %>%
  select(SPECIES.TREE.LABEL, LTP)

files = list.files(path = tree_path)
markers = lapply(files, (function(f) strsplit(f, split="[.]")[[1]][1]) )

flist = as.list(files)
names(flist) = unlist(markers)

for (f in flist) {
  t = read.newick(f)
  t$tip.label = sapply(as.list(t$tip.label), nametips, isolates = isolates)
  write.tree(t, file = gsub("[.]tre", ".tre.renamed", f))
}
