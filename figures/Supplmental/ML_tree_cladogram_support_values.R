library(tidyverse)
library(ggtree)

tree = read.newick("/home/aimzez/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_all_support.tre")

plt_tree = ggtree(tree, branch.length = "none")

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.double) %>%
  mutate(gcf = gcf/100) %>%
  mutate(gcf = ifelse(is.na(gcf), 0, gcf)) %>%
  mutate(gcf = round(gcf,2)) %>%
  mutate(qpic = round(qpic,2)) %>%
  unite(col = "support_above", bootstrap, qpic, sep = "/", remove = F) %>%
  unite(col = "support_below", gcf, astral, sep = "/", remove = F)

clado = ggtree(tree, branch.length = "none")  %>%
  ggtree::rotate(node = 146) %>%
  ggtree::rotate(node=225) %>%
  ggtree::rotate(node = 228)
clado = clado %<+% tree_data

clado = clado +
  geom_tiplab(cex=2.5) +
  geom_nodelab(aes(x = x-0.5, y=y+0.5, label=support_above), size=1.5) +
  geom_nodelab(aes(x = x-0.5, y=y-0.5, label=support_below), size=1.5) +
  scale_x_continuous(limits = c(0,28))

clado_cleaned_tiplabs = clado
clado_cleaned_tiplabs$data = clado_cleaned_tiplabs$data %>%
  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
  mutate(label = gsub("LCG$", "", label)) %>%
  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = clado_cleaned_tiplabs$data %>%
  filter(isTip) %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

clado_cleaned_tiplabs$data = clado_cleaned_tiplabs$data %>%
  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(label = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), label))

ggsave(filename='~/work/Chytrid-Phylogenomics/figures/Supplmental/ML_tree_cladogram_support.pdf',
       plot = clado_cleaned_tiplabs, device = "pdf", width = 8.5, height = 11, units = "in")

