library(phytools)
library(tidyverse)
library(readxl)
library(argparse)

nametips = function(lab, isolates) {
  spl = strsplit(lab, split="[|]")
  new = isolates %>%
    filter(LTP == spl[[1]][1]) %>%
    select (SPECIES.TREE.LABEL) %>%
    .$SPECIES.TREE.LABEL
  paste(new, lab, sep="&")
}

parser = argparse::ArgumentParser()
parser$add_argument("-g", "--genetrees", action = "store", required = T, help = "Path to direction containing the gene trees `to rename.")
parser$add_argument("-i", "--isolates", action = "store", required = T, help = "Path to isolates sheet with `SPECIES.TREE.LABEL` and `LTP` columns.")
parser$add_argument("--suffix", action = "store", required = T, help = "Suffix of gene tree files.")

args = parser$parse_args()
tree_path = args$genetrees
isolates_path = args$isolates
suffix = paste(args$suffix, sep ="")

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
  write.tree(t, file = gsub(suffix, paste(suffix, ".renamed", sep=""), f, fixed = T))
}
