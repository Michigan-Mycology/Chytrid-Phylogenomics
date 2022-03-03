library(tidyverse)
library(phytools)
library(ggtree)
#library(tublerone)
library(readxl)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(scales)
library(patchwork)

#### Functions ####
rename_ultra = function(gs_combined, genus, ultra) {
  is_gs = gs_combined %in% ultra$Taxon
  if (is_gs) {
    ultraTax = ultra %>%
      filter(Taxon == gs_combined) %>%
      .$Taxon

    return(ultraTax)

  } else {
    is_g = genus %in% ultra$Taxon
    if(is_g) {
      ultraTax = ultra %>%
        filter(Taxon == genus) %>%
        .$Taxon

      return(ultraTax)
    } else {
      return("NO")
    }
  }
}

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
  limits <- range(ggtree$data$y, na.rm = T)
  print(limits)
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

#### Some important variables ####
raberns_prettiest = c("N", "2", "7", "15", "19", "20")

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet=1)

#### Prepare Ultrasctructural Character Table for Tree ####
ultra = read_xlsx("~/work/pursuit/sheets/Phylogenomics_Ultrastructure.xlsx") %>%
  mutate(Taxon = gsub("[ ]","_", Taxon))
ultra_merger = isolates %>%
  select(SPECIES.TREE.LABEL) %>%
  mutate(hold = SPECIES.TREE.LABEL) %>%
  separate(SPECIES.TREE.LABEL, into=c("genus", "species"), sep="[_]") %>%
  mutate(genus_merge = genus) %>%
  mutate(species_merge = species) %>%
  unite(gs_combined, c(genus_merge,species_merge), sep="_")

ultra_merged = ultra_merger %>%
  rowwise() %>%
  mutate(Taxon = rename_ultra(gs_combined, genus, ultra)) %>%
  select(hold, Taxon) %>%
  rename(SPECIES.TREE.LABEL = hold) %>%
  left_join(ultra, by="Taxon") %>%
  select(-Taxon) %>%
  #select(SPECIES.TREE.LABEL, all_of(raberns_prettiest)) %>%
  rename(
    "Electron-opaque plugs in flagellum" = "2",
    "Kinetosome-associated plates" = "7",
    "Non-flagellated centriole positioning" = "15",
    "Paracrystalline inclusion" = "19",
    "Ribosomal organization" = "20"
  )

renamed = final$data %>%
  filter(isTip) %>%
  select(label) %>%
  rename(SPECIES.TREE.LABEL = label) %>%
  left_join(ultra_merged)

write.table(renamed, file = "~/work/pursuit/sheets/Phylogenomics_Ultrastructure_renamed.tsv", quote=F, row.names = F, sep="\t")

ultra_merged_long = ultra_merged %>%
  gather(key="char", value="state", -SPECIES.TREE.LABEL)

#### Read and plot Tree ####
tree = read.newick("/home/aimzez/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_filtered_support.tre")

plt_tree = ggtree(tree)

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.numeric) %>%
  mutate(gcf = gcf/100)

final = ggtree(tree) %>%
  collapse(node=229) %>%
  collapse(node=241) %>%
  collapse(node=272)
final = final %<+% tree_data

final = final +
  geom_text2(aes(subset=(node == 229)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 229)), label = "Dikarya", cex=3.0, vjust=0.4, hjust = -0.5) +
  geom_text2(aes(subset=(node == 241)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 241)), label = "Mucoromycota", cex=3.0, vjust=0.4, hjust = -0.25) +
  geom_text2(aes(subset=(node == 272)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 272)), label = "Non-Fungal Eukaryota", cex=3.0, hjust =-0.2,vjust=.45)

x_max = max(final$data$x[!is.na(final$data$x)])
final = final +
  theme_tree2() +
  geom_tiplab(cex=2.7)

horiz_just = 35
vertical_just = 1
nodelab_pos = as_tibble(final$data) %>%
  select(node, label, branch.length, x, y, isTip) %>%
  filter(!isTip) %>%
  mutate(repos_node_label = ifelse(branch.length < 35, TRUE, FALSE)) %>%
  mutate(node_x = x - x_max) %>%
  mutate(node_y = y) %>%
  mutate(nodelab_x = ifelse(repos_node_label, x - x_max - horiz_just, x - x_max)) %>%
  mutate(nodelab_y = ifelse(repos_node_label, y + vertical_just, y)) %>%
  select(-x, -y, -isTip, -label, -branch.length)
final$data = final$data %>%
  left_join(nodelab_pos, by="node")

final = final +
  geom_point(data = subset(final$data, !isTip & !is.na(qpic)), aes(x=nodelab_x, y=nodelab_y, fill=qpic), color="black", pch=21, cex=2.0) +
  geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=nodelab_x-20, y=nodelab_y, fill=gcf), color="black", pch=22, cex=2.0) +

  # Shown only where clade was not represented in astral tree
  geom_point(data = subset(final$data, !isTip & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW" & is.na(astral)), aes(x=nodelab_x-40, y=nodelab_y, fill=astral), color="black", pch=23, cex=2.0) +
  scale_fill_continuous(type="viridis", na.value="red")

final = revts(final) +
  geom_segment(data = subset(final$data, !isTip & repos_node_label & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW"), aes(x=node_x, y=node_y, xend=nodelab_x, yend=nodelab_y-0.35), alpha=1.0) +
  #geom_text(data = subset(final$data, !isTip), aes(x=nodelab_x, y=nodelab_y, label=label), size=1.8, color=muted("blue")) +
  scale_x_continuous(limits = c(-1100,800), labels = abs, breaks = c(-1000, -750, -500, -250, 0)) +
  scale_y_continuous(expand=expansion(0, 0.6)) +
  xlab("Millions of Years Ago (Mya)") +
  guides(
    fill = FALSE
  )

final

#### Make ploidy track
#### Map ploidy on to tips ####
ploidy = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP, STRAIN, ploidy_file_prefix, "ploidy (12-11-20)") %>%
  rename(inferred_ploidy = "ploidy (12-11-20)") %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "Not_Yet_Inferred"))) %>%
  select(SPECIES.TREE.LABEL, inferred_ploidy) %>%
  rename(label = SPECIES.TREE.LABEL)

ploidy_colors = c(
  "blue",
  "red",
  muted("blue"),
  muted("red"),
  "yellow",
  "grey"
)

ploidy_track = ggtreeplot(final, ploidy, aes(x=1, y=y, fill=inferred_ploidy), flip=F, expand_limits = expansion(0,0.6)) +
  geom_point(pch=21, color="black", size=3) +
  #geom_tile() +
  scale_fill_manual(values=ploidy_colors) +
  theme_minimal() +
  guides(fill=F) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

#### Add simpler characters like zoosporic? nutritional mode? Development? Ef1a-like protein?

#### Make list of ultrastructure tracks ####
final_pts = final
columns = list()
i=1
pals = c("YlOrRd", "BrBG", "BuGn", "Set1", "Accent")
pch = c(21,22,23,24,25)
for (c in colnames(ultra_merged)[-1]) {
  sub = ultra_merged_long %>%
    filter(char == c) %>%
    rename(label = SPECIES.TREE.LABEL)
  plt_tbl = final$data %>%
    filter(isTip == T) %>%
    select(label, y) %>%
    left_join(sub)

  colscale = brewer.pal(length(unique(plt_tbl %>% .$state))-1, pals[i])
  col = ggtreeplot(final, plt_tbl, aes(x=1, y=y, fill=state), flip=FALSE, expand_limits = expansion(0,0.6)) +
    geom_tile(pch = pch[i], color="black") +
    coord_fixed()+
    theme_minimal() +
    guides() +
    scale_fill_manual(values=colscale)+
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      panel.grid = element_blank(),
      plot.margin = unit(c(0,0,0,0),"pt"),
      panel.spacing = unit(0, "cm"),
      panel.background = element_blank(),
      #plot.background = element_rect(fill="lightblue")
    )
  columns[[i]] = col
  i = i +1
}

columns[[1]]+columns[[2]]+columns[[3]]+columns[[4]]+columns[[5]] + plot_layout(nrow=1)

final +
  ploidy_track +
  columns[[1]] +
  columns[[2]] +
  columns[[3]] +
  columns[[4]] +
  columns[[5]] +
  plot_layout(nrow=1, widths=c(15,1,0.5,0.5,0.5,0.5,0.5))

pts = final
for (c in colnames(ultra_merged)[-1]) {
  print(c)
  pts = pts +
    geom_facet(panel = c, data = subset(ultra_merged_long, char == c), geom = geom_point, mapping=aes(x = 1, color = state))
}

gheatmap(plt_tree, ultra_merged_hm, offset=0.5, legend_title = "Character State", color=NULL) +
  scale_fill_manual(values=pal)
