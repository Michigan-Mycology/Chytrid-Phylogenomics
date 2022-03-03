library(tidyverse)
library(ggtree)
library(phytools)

astral_tree = read.newick("~/DATA/pursuit/astral_tree/ASTRAL_20201027.noBrlen.tre.renamed")
astral_tree = ape::root(astral_tree, outgroup=c("Capsaspora_owczarzaki_ATCC_30864.v2",
                                                "Monosiga_brevicolis_MX1.v1",
                                                "Drosophila_melanogaster.v6"
))

astral = ggtree(astral_tree, branch.length = "none") %>%
  ggtree::rotate(node = 186) %>%
  ggtree::rotate(node = 145) %>%
  ggtree::rotate(node = 156) %>%
  ggtree::rotate(node = 169)
astral = astral +
  geom_tiplab(cex=2.5) +
  geom_text(data=subset(astral$data, !isTip), aes(x=x-0.5, y=y+0.5, label=label), cex=2.0) +
  scale_x_continuous(limits=c(0,35))


astral_cleaned_tiplabs = astral
astral_cleaned_tiplabs$data = astral_cleaned_tiplabs$data %>%
  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
  mutate(label = gsub("LCG$", "", label)) %>%
  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = astral_cleaned_tiplabs$data %>%
  filter(isTip) %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

astral_cleaned_tiplabs$data = astral_cleaned_tiplabs$data %>%
  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(label = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), label))

ggsave(filename = "~/work/Chytrid-Phylogenomics/figures/Supplmental/ASTRAL_tree.pdf",
       plot = astral_cleaned_tiplabs, device = "pdf", width = 8.5, height = 11, units = "in")
