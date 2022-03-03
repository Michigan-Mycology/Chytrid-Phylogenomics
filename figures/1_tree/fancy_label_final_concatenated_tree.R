library(tidyverse)
library(ggtree)
library(phytools)
library(scales)

tree = read.newick("~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_all_support.tre")
tree$node.label
n_tips = length(tree$tip.label)

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

plt_tree = ggtree(tree)

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.numeric) %>%
  mutate(gcf = gcf/100)

final = ggtree(tree) %<+% tree_data
final = final +
  theme_tree2() +
  geom_tiplab(cex=2.5) +
  geom_point(data = subset(final$data, !isTip & !is.na(qpic)), aes(x=x-max(x), y=y, fill=qpic), color="black", pch=21, cex=2.0) +
  geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=x-max(x)-20, y=y, fill=gcf), color="black", pch=22, cex=2.0) +
  geom_point(data = subset(final$data, !isTip & !is.na(astral)), aes(x=x-max(x)-40, y=y, fill=astral), color="black", pch=23, cex=2.0) +

  scale_fill_continuous(type="viridis")

  #geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=x-max(x)-20, y=y), fill="black", color="black", pch=24)

final = revts(final) +
  scale_x_continuous(limits = c(-1100,500), labels = abs)
  #geom_nodepoint(aes(color = qpic)) +
final

ggsave(plot = final, file = "~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_nodepoint.pdf", device = "pdf", width = 8.5, height = 11, units = "in")

