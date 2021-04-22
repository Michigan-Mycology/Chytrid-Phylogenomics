library(phytools)
library(tidyverse)
library(readxl)

nametips = function(lab) {
  spl = strsplit(lab, split="[|]")
  spl[[1]][1]
}

args = commandArgs(trailingOnly = TRUE)
tree_path = args[1]

setwd(tree_path)

files = list.files(path = tree_path)
markers = lapply(files, (function(f) strsplit(f, split="[.]")[[1]][1]) )

flist = as.list(files)
names(flist) = unlist(markers)

for (f in flist) {
  t = read.newick(f)
  t$tip.label = sapply(as.list(t$tip.label), nametips)
  write.tree(t, file = gsub("[.]contree", ".contree.strip_renamed", f))
}
