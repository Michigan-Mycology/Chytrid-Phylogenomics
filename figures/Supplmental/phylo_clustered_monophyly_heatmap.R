library(phytools)
library(tidyverse)
library(readxl)
library(magrittr)
library(scales)
library(ggtree)
library(patchwork)

### Adapted from https://thackl.github.io/ggtree-composite-plots, thackl on GitHub
tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}

ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip=FALSE,
                       expand_limits=expansion(0,.6), ...){

  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")

  # match the tree limits
  print(expand_limits)
  limits <- range(ggtree$data$y, na.rm = T)
  limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
  limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]

  print(limits)

  if(flip){
    mapping <- modifyList(aes_(x=~x), mapping)
    data <- mutate(data, x=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_x_continuous(limits=limits, expand=c(0,0))
  }else{
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  gg
}

#### Tree ####
tree = read.newick("/home/aimzez/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_filtered_support.tre")

plt_tree = ggtree(tree)

#plt_tree %>%
#  ggtree::rotate(node = 146) %>%
#  ggtree::rotate(node=225) %>%
#  ggtree::rotate(node = 228) +
#  geom_text(aes(label=node)) +
#  geom_tiplab(cex=2.5)

final = ggtree(tree)  %>%
  ggtree::rotate(node = 146) %>%
  ggtree::rotate(node=225) %>%
  ggtree::rotate(node = 228)
  #collapse(node=229) %>%
  #collapse(node=241) %>%
  #collapse(node=272)
final

##########################
### Prefilter ####
##########################

matrix_path = "~/DATA/pursuit/phylogeny_trimfix/round1_raw/monophyly_matrix.csv"
isolates_xlsx = "~/work/pursuit/sheets/Pursuit_Isolates.xlsx"

mat = read_delim(matrix_path, delim=",")

colnames(mat)[1] = "LTP"
isolates = read_xlsx(isolates_xlsx) %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>%
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP)

named_mat %<>%  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = as.character(value)) %>%
  mutate(value = as.factor(value))

named_mat$value = factor(named_mat$value, levels=c("-1","0","1"))
named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

named_mat = named_mat %>%
  rename(label = SPECIES.TREE.LABEL)


colors = c("grey", muted("red"), muted("blue"))
names(colors) = levels(named_mat$value)
hm_orig = ggtreeplot(final, named_mat, aes(x=marker, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=colors) +
  guides(fill =F) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

##########################
### Postfilter ####
##########################

matrix_path = "~/DATA/pursuit/phylogeny_trimfix/round6_rmpoly/mafft_monophyly_matrix.csv"
isolates_xlsx = "~/work/pursuit/sheets/Pursuit_Isolates.xlsx"

mat = read_delim(matrix_path, delim=",")

colnames(mat)[1] = "LTP"
isolates = read_xlsx(isolates_xlsx) %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>%
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP)

named_mat %<>%  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = as.character(value)) %>%
  mutate(value = as.factor(value))

named_mat$value = factor(named_mat$value, levels=c("-1","0","1"))
named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

named_mat = named_mat %>%
  rename(label = SPECIES.TREE.LABEL)


colors = c("grey", muted("red"), muted("blue"))
names(colors) = levels(named_mat$value)
hm_post = ggtreeplot(final, named_mat, aes(x=marker, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=colors) +
  guides(fill =F) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

combined = final + hm_orig + hm_post + plot_layout(widths=c(2,3,3))
ggsave(filename = "~/Dropbox (University of Michigan)/dissertation/dissertation_talk/clustered_phyly_hm.pdf", plot = combined, device=cairo_pdf, width = 13.2, height = 6.5, units = "in")
