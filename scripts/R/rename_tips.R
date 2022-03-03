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

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP)

f = tree_path

t = read.newick(f)
t$tip.label = sapply(as.list(t$tip.label), nametips, isolates = isolates)
write.tree(t, file = gsub(".treefile", ".treefile.renamed", f))
