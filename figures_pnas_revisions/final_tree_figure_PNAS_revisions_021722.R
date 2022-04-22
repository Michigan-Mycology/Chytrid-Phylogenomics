library(tidyverse)
library(readxl)
library(ggtree)
library(phytools)
library(scales)
library(RColorBrewer)
library(patchwork)
library(ggstance)
library(tublerone)
library(ggtreeExtra)

CHYTRID_PHYLO="/home/aimzez/work/Chytrid-Phylogenomics"

#### Tree ####
tree = read.newick(file.path(CHYTRID_PHYLO, "figures/1_tree", "combined_tree_filtered_support_RENAMED.tre"))

plt_tree = ggtree(tree)

#plt_tree %>%
#  ggtree::rotate(node = 146) %>%
#  ggtree::rotate(node=225) %>%
#  ggtree::rotate(node = 228) +
#  geom_text(aes(label=node)) +
#  geom_tiplab(cex=2.5)

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.numeric) %>%
  mutate(gcf = gcf/100) %>%
  mutate(gcf = ifelse(is.na(gcf), 0, gcf))

final = ggtree(tree)  %>%
  ggtree::rotate(node = 146) %>%
  ggtree::rotate(node=225) %>%
  ggtree::rotate(node = 228) %>%
  collapse(node=229) %>%
  collapse(node=241) %>%
  collapse(node=272) %>%
  collapse(node=256)
final = final %<+% tree_data

final = final +
  geom_text2(aes(subset=(node == 229)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
  geom_text2(aes(subset=(node == 229)), label = "Dikarya", cex=3.0, vjust=0.4, hjust = -0.5, color = "black") +
  geom_text2(aes(subset=(node == 241)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
  geom_text2(aes(subset=(node == 241)), label = "Mucoromycota", cex=3.0, vjust=0.4, hjust = -0.5, color = "black") +
  geom_text2(aes(subset=(node == 272)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
  geom_text2(aes(subset=(node == 272)), label = "Non-Fungal Eukaryota", cex=3.0, hjust =-0.5,vjust=.45, color = "black") +
  geom_text2(aes(subset=(node == 256)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
  geom_text2(aes(subset=(node == 256)), label = "Zoopagomycota", cex=3.0, hjust =-0.5,vjust=.45, color = "black")

x_max = max(final$data$x[!is.na(final$data$x)])
final = final +
  theme_tree2() +
  geom_tiplab(size=2.3, color = 'black')

#### The offset (sometimes) node point stuff ####
#horiz_just = 35
#vertical_just = 1
#nodelab_pos = as_tibble(final$data) %>%
#  select(node, label, branch.length, x, y, isTip) %>%
#  filter(!isTip) %>%
#  mutate(repos_node_label = ifelse(branch.length < 35, TRUE, FALSE)) %>%
#  mutate(node_x = x - x_max) %>%
#  mutate(node_y = y) %>%
#  mutate(nodelab_x = ifelse(repos_node_label, x - x_max - horiz_just, x - x_max)) %>%
#  mutate(nodelab_y = ifelse(repos_node_label, y + vertical_just, y)) %>%
#  select(-x, -y, -isTip, -label, -branch.length)
#final$data = final$data %>%
#  left_join(nodelab_pos, by="node")

#final = final +
  #geom_point(data = subset(final$data, !isTip & !is.na(qpic)), aes(x=nodelab_x, y=nodelab_y, fill=gcf), color="black", pch=21, cex=2.0) +
  #geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=nodelab_x-20, y=nodelab_y, fill=gcf), color="black", pch=22, cex=2.0) +

  #Shown only where clade was not represented in astral tree
  #geom_point(data = subset(final$data, !isTip & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW" & is.na(astral)), aes(x=nodelab_x-20, y=nodelab_y, fill=astral), color="black", pch=23, cex=2.0) +
  #scale_fill_continuous(type="viridis", na.value="red")

######### Node pos for just astral points #########
nodelab_pos = as_tibble(final$data) %>%
  select(node, label, branch.length, x, y, isTip) %>%
  filter(!isTip) %>%
  mutate(node_x = x - x_max) %>%
  mutate(node_y = y) %>%
  select(-x, -y, -isTip, -label, -branch.length)
final$data = final$data %>%
  left_join(nodelab_pos, by="node")
###################################################

### Add ASTRAL points
final = final +
  geom_point(data = subset(final$data, !isTip & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW" & is.na(astral)), aes(x=node_x, y=node_y, fill=astral), color="black", pch=23, cex=2.0) +
  scale_fill_continuous(type="viridis", na.value="red")

#### Color edge by GCF ####
#final = final +
#  aes(color=gcf) +
#  scale_color_gradient(low = "blue", high = "orange", na.value="black")
###########################

#### Annotate line width by GCF ####
final = final +
  aes(lwd=gcf) +
  #geom_text(aes(label=gcf)) +
  scale_size_continuous(range = c(0.5,1.25))

####################################

final = revts(final) +
  #geom_segment(data = subset(final$data, !isTip & repos_node_label & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW"), aes(x=node_x, y=node_y, xend=nodelab_x, yend=nodelab_y-0.35), alpha=1.0) +
  #geom_text(data = subset(final$data, !isTip), aes(x=nodelab_x, y=nodelab_y, label=label), size=1.8, color=muted("blue")) +
  scale_x_continuous(limits = c(-1100,800),
                     labels = c("1000", rep("",4), "750", rep("",4), "500", rep("",4), "250", rep("",4),"0"),
                     breaks = seq(-1000,0,50),
                     expand=expansion(0,0)) +
  scale_y_continuous(expand=expansion(0, 0.6)) +
  xlab("Millions of Years Ago (Mya)") +
  guides(
    fill = FALSE,
    color = F,
    size = F
  ) +
  theme(
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.length.x = unit(0, "pt"),
    axis.title.x = element_blank(),
    axis.line.x = element_line(),
    plot.margin = unit(c(0,0,0,0),"cm"),
  )

final

#### Geological Time Scale ####
geodates = read_delim(file.path(CHYTRID_PHYLO, "figures/1_tree", "geologic_time_scale_dates.csv"), delim=',')

geoplot = geodates %>%
  select(eon, period, start_mya) %>%
  group_by(eon, period) %>%
  summarise(start_mya = max(start_mya)) %>%
  ungroup %>%
  arrange(start_mya) %>%
  mutate(start = 0-start_mya) %>%
  mutate(end = lag(start)) %>%
  mutate(period_ymin = 1, period_ymax = 2) %>%
  mutate(period_disp = str_sub(period, 1, 1) ) %>%
  mutate(period_disp = ifelse(period == "Quaternary", "", period_disp)) %>%
  mutate(period_disp = ifelse(period == "Neogene", "", period_disp)) %>%
  mutate(end = ifelse(period == "Quaternary", 0, end)) %>%
  mutate(period_label_pos_y = (period_ymin+period_ymax)/2, period_label_pos_x = (start+end)/2) %>%

  # Filter out time before tree root
  filter(end > -1000) %>%

  # Change start dates of stuff on edge of past cutoff (i.e 1000 MYA)
  mutate(start = ifelse(start <= -1000, -1000, start))

eons = geoplot %>%
  group_by(eon) %>%
  summarise(start = min(start), end = max(end)) %>%
  mutate(eon_ymin = 0, eon_ymax = 1) %>%
  mutate(eon_label_pos_y = (eon_ymin+eon_ymax)/2, eon_label_pos_x = (start+end)/2)

period_colors = rep(c("white", "grey"),8)[-16]
names(period_colors) = geoplot$period
period_colors = as.list(period_colors)
period_colors[["Proterozoic"]] = "grey20"
period_colors[["Phanerozoic"]] = "grey20"

timescale = ggplot(geoplot) +
  geom_rect(aes(xmin = start, xmax = end, ymin = period_ymin, ymax = period_ymax, fill = period), color = "black") +
  geom_rect(data = eons, aes(xmin = start, xmax = end, ymin = eon_ymin, ymax = eon_ymax, fill = eon), color = "black", inherit.aes = F) +
  geom_text(aes(label = period_disp, x = period_label_pos_x, y = period_label_pos_y), size=2)+
  geom_text(data = eons, aes(label=eon, x=eon_label_pos_x, y=eon_label_pos_y), color="white", size =2) +
  scale_x_continuous(limits = c(-1100,800), expand = c(0,0)) +
  scale_fill_manual(values=period_colors) +
  guides (fill=F) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0),"cm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt")
  )

tree_ylim = c(
  min(final$data$y, na.rm = T),
  max(final$data$y, na.rm = T)
)
geoplot_sub = geoplot %>%
  filter(period %in% c("Cryogenian", "Cambrian", "Silurian", "Carboniferous", "Triassic", "Cretaceous", "Neogene"))
final = final +
  geom_rect(data = geoplot_sub, aes(xmin=start, xmax=end, ymin=tree_ylim[1], ymax=tree_ylim[2]), inherit.aes = F, alpha = 0.2)
final+ timescale + plot_layout(nrow=2, heights=c(15,1))

final_cleaned_tiplabs = final
final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
  mutate(label = gsub("LCG$", "", label)) %>%
  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = final_cleaned_tiplabs$data %>%
  filter(isTip) %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(label = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), label))

### Just print the tree and timescale, sans right-side data
just_the_tree = final_cleaned_tiplabs + timescale + plot_layout(nrow=2, heights=c(15,1))
ggsave(
  filename = file.path(CHYTRID_PHYLO, "figures_pnas_revisions/1_tree", "Tree_Figure_2_042222.pdf"),
  plot = just_the_tree,
  device=cairo_pdf,
  width = 8.5,
  height = 11,
  units = "in")


################################################
################## FIGURE 2 ####################
################################################

#### Tree ####
isolates = read_xlsx("~/data/pursuit/Pursuit_Isolates.xlsx")

tree = read.newick(file.path(CHYTRID_PHYLO, "figures/1_tree", "combined_tree_filtered_support_RENAMED.tre"))

old_to_new = isolates$SPECIES.TREE.LABEL
names(old_to_new) = isolates$LTP

#tree = tublerone::phylo_rename_tips(tree, old_to_new)
tree = tublerone::phylo_root(tree, outgroup = c("Capsaspora_owczarzaki_ATCC_30864.v2", "Drosophila_melanogaster.v6", "Monosiga_brevicolis_MX1.v1"))
tree$tip.label
plt_tree = ggtree(tree, open.angle = 0.15) %>%
  ggtree::rotate(node = 145) %>%
  ggtree::rotate(node = 144) %>%
  ggtree::rotate(node = 234)
plt_tree = plt_tree +
  geom_tiplab(size =2, color = 'black') +
  geom_text(aes(x=x, y=y, label=node), color='red')
plt_tree
#plt_tree %>%
#  ggtree::rotate(node = 146) %>%
#  ggtree::rotate(node=225) %>%
#  ggtree::rotate(node = 228) +
#  geom_text(aes(label=node)) +
#  geom_tiplab(cex=2.5)

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.numeric) %>%
  mutate(gcf = gcf/100) %>%
  mutate(gcf = ifelse(is.na(gcf), 0, gcf))


#%%
######
final = ggtree(tree, layout = "fan", size = 0.50, open.angle = 80)  %>%
  ggtree::rotate(node = 145) %>%
  ggtree::rotate(node = 144) %>%
  ggtree::rotate(node = 234)
  #collapse(node=235) %>%
  #collapse(node=247) %>%
  #collapse(node=273) #%>%
  #collapse(node=256)
final = final %<+% tree_data

final = final +
  geom_highlight(node = 155, extend =1000, alpha =0.1, fill = "darkgoldenrod1", color = "black") + #chytrid
  geom_highlight(node = 148, extend =1000, alpha =0.1, fill = "red", color = "black") + #neos
  geom_highlight(node = 266, extend =1000, alpha =0.1, fill = "green", color = "black") + #blastos
  geom_highlight(node = 127, extend =1000, alpha =0.1, fill = "darkturquoise", color = "black") + #olpidium
  geom_highlight(node = 256, extend =1000, alpha =0.1, fill = "pink", color = "black") + #zoopags
  geom_highlight(node = 229, extend =1000, alpha =0.1, fill = "darkblue", color = "black") + #dikarya
  geom_highlight(node = 241, extend =1000, alpha =0.1, fill = "cadetblue1", color = "black") + #mucoro
  geom_highlight(node = 7, extend =1000, alpha =0.1, fill = "coral2", color = "black") #paraphelidium

#final = final +
#  geom_text2(aes(subset=(node == 229)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
#  geom_text2(aes(subset=(node == 145)), label = "Dikarya", cex=3.0, vjust=0.4, hjust = -0.5, color = "black") +
#  geom_text2(aes(subset=(node == 241)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
#  geom_text2(aes(subset=(node == 144)), label = "Mucoromycota", cex=3.0, vjust=0.4, hjust = -0.5, color = "black") +
#  geom_text2(aes(subset=(node == 272)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
#  geom_text2(aes(subset=(node == 273)), label = "Non-Fungal Eukaryota", cex=3.0, hjust =-0.5,vjust=.45, color = "black") #+
#  geom_text2(aes(subset=(node == 256)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45, color = "black") +
#  geom_text2(aes(subset=(node == 256)), label = "Zoopagomycota", cex=3.0, hjust =-0.5,vjust=.45, color = "black")

x_max = max(final$data$x[!is.na(final$data$x)])
final = final +
  theme_tree2() #+
  #geom_tiplab(size=2.3, color = 'black')

#### Character heatmap ####
character_sheet = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, Flagellar.State, `Electron-opaque_plug_in_axoneme_core_and_between_axoneme_and_flagellar_membrane (2)`,
         `Mode (M)`, `Nutritional_Mode (N)`, EFL, EF1a, Cbl_group, MMCoA_group, MetH, Selenocystein,
         Ancestral_combined:Whi5_Fungal) %>%
  rename(opaque_flag_plug = `Electron-opaque_plug_in_axoneme_core_and_between_axoneme_and_flagellar_membrane (2)`) %>%
  mutate_at(vars(Ancestral_combined:Whi5_Fungal), ~ifelse(is.na(.), 0, .)) %>%
  mutate(EFL = ifelse(EFL == 0, "EFL-", "EFL+")) %>%
  mutate(EF1a = ifelse(EF1a == 0, "EF1a-", "EF1a+")) %>%
  mutate(E2F_Ancestral = ifelse(E2F_Ancestral == 0, "EF2-", "EF2+")) %>%
  mutate(Rb_Ancestral = ifelse(Rb_Ancestral == 0, "Rb-", "Rb+")) %>%
  mutate(SBF_Fungal = ifelse(SBF_Fungal == 0, "SBF-", "SBF+")) %>%
  mutate(Whi5_Fungal = ifelse(Whi5_Fungal == 0, "Whi5-", "Whi5+")) %>%
  mutate(Flagellar.State = ifelse(Flagellar.State == 0, "not_zoosporic", "zoosporic")) %>%
  mutate(opaque_flag_plug = ifelse(opaque_flag_plug == 0, "OFP-", "OFP+")) %>%
  mutate(Selenocystein = ifelse(Selenocystein == 0, "Sel-", "Sel+")) %>%
  select(-Ancestral_combined, -Fungal_combined)

character_sheet_long = character_sheet %>%
  gather(key="char", value="state", -SPECIES.TREE.LABEL) %>%
  mutate(state = as.factor(state)) %>%
  rename(ID = SPECIES.TREE.LABEL)
character_sheet_long = character_sheet_long %>%
  mutate(char = factor(character_sheet_long$char, levels = names(pals)))

unipal = c("white", "mediumpurple1", "purple", "magenta4", "wheat1", "turquoise2",
           "white", "blue", "white", "darkorange3", "white", "yellow", "green", "springgreen1",
           "white", "white", "seagreen4", "tomato1", "white", "darkorange3", "violetred1",
           "sienna4", "white", "darkorange1", "white", "darkolivegreen", "white", "darkorange1", "red")

unipal_reduced = c("white", "mediumpurple1", "purple", "magenta4",
           "white", "blue", "white", "darkorange3", "white", "yellow",
           "white", "white", "darkorange3", "white", "darkorange1",
           "white", "darkolivegreen", "white", "darkorange1", "red", "red")

character_sheet_long = character_sheet_long %>%
  filter(char %in% c("Flagellar.State", "Mode (M)", "EFL", "EF1a",
                     "Cbl_group", "MMCoA_group", "MetH", "Selenocystein",
                     "E2F_Ancestral", "Rb_Ancestral", "SBF_Fungal", "Whi5_Fungal"))
character_sheet_long = character_sheet_long %>%
  mutate(char = factor(character_sheet_long$char, levels = c("Flagellar.State", "Mode (M)", "EFL", "EF1a",
                                                             "Cbl_group", "MMCoA_group", "MetH", "Selenocystein",
                                                             "E2F_Ancestral", "Rb_Ancestral", "SBF_Fungal", "Whi5_Fungal")))


fruitdata1 = character_sheet_long %>% filter( char %in% c("Flagellar.State", "EFL", "EF1a"))
fruitdata1 = fruitdata1 %>%
  mutate(char = factor(fruitdata1$char, levels = c("Flagellar.State", "EFL", "EF1a")))

fig2 = final +
  geom_fruit(data = fruitdata1,
             geom = geom_tile,
             mapping = aes(y= ID, x = char, fill = state),
             stat = "identity",
             color = "black",
             offset = 0.05,
             size = 0.5,
             pwidth = 0.1,
             alpha = 1.0) +
  geom_fruit(data = character_sheet_long %>% filter(char %in% c("Cbl_group", "MMCoA_group", "MetH", "Selenocystein")),
             geom = geom_tile,
             mapping = aes(y= ID, x = char, fill = state),
             stat = "identity",
             color = "black",
             offset = 0.05,
             size = 0.5,
             pwidth = 0.2,
             alpha = 1.0) +
  geom_fruit(data = character_sheet_long %>% filter(char %in% c("E2F_Ancestral", "Rb_Ancestral", "SBF_Fungal", "Whi5_Fungal")),
             geom = geom_tile,
             mapping = aes(y= ID, x = char, fill = state),
             stat = "identity",
             color = "black",
             offset = 0.05,
             size = 0.5,
             pwidth = 0.2,
             alpha = 1.0) +
  scale_fill_manual(values = unipal_reduced) +
  guides(fill = "none") +
  theme(
    axis.text = element_blank(),
    axis.line.x = element_blank()
  )
fig2


geom_fruit(data = character_sheet_long %>% filter(!char == "Mode (M)"),
           geom = geom_tile,
           mapping = aes(y= ID, x = char, fill = state),
           stat = "identity",
           color = "black",
           offset = 0.02,
           size = 0.5,
           pwidth = 0.6,
           alpha = 1.0)


ggsave(
  filename = file.path(CHYTRID_PHYLO, "figures_pnas_revisions/2_character_states", "CharStates_Fig_2_021722.pdf"),
  plot = fig2,
  device=cairo_pdf,
  width = 8.5,
  height = 11,
  units = "in")


#### SNP density bar charts ####

ploidy_df = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, coding, known_diploid_mitosis, UM_ploidy) %>%
  rename(ploidy = coding) %>%
  mutate(ploidy= ifelse(ploidy == 1, "haploid", "diploid"))

l50_genome_sizes = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, l50_assembly_length)
isolates = read_xlsx(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Isolates.xlsx"))

snp_densities = read_delim(file.path(CHYTRID_PHYLO, "figures/1_tree", "all.snp_contig_counts_strainified.tsv"), delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(num_snps = sum(num_snps)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  left_join(l50_genome_sizes) %>%
  select(SPECIES.TREE.LABEL, num_snps, l50_assembly_length) %>%
  mutate(snp_density = num_snps/l50_assembly_length) %>%
  left_join(ploidy_df %>% select(SPECIES.TREE.LABEL, ploidy)) %>%
  rename(ID = SPECIES.TREE.LABEL)

fig2_ploidy = fig2 +
  geom_fruit(data = snp_densities, geom = geom_bar,
             mapping = aes(y = ID, x = num_snps, color = ploidy),
             fill = "lightskyblue",
             offset = 0.10,
             stat = "identity",
             pwidth = 0.65,
             orientation = "y") +
  scale_color_manual(values =c('blue', 'red')) +
  guides(color="none")
ggsave(
  filename = file.path(CHYTRID_PHYLO, "figures_pnas_revisions/2_character_states", "CharStates_Fig_2_021722_withPloidy_snpDense.pdf"),
  plot = fig2_ploidy,
  device=cairo_pdf,
  width = 8.5,
  height = 11,
  units = "in")



#### Genome size bubbles ####
genome_sizes = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, assembly_length) %>%
  rename(label = SPECIES.TREE.LABEL)

plt_tbl = final$data %>%
  filter(isTip == T) %>%
  select(label, y) %>%
  left_join(genome_sizes) %>%
  mutate(assembly_length = assembly_length/1e6)

genome_size_column = tublerone::ggtreeplot(final, plt_tbl, aes(x=1, y=y, size=assembly_length), flip=FALSE, expand_limits = expansion(0.00,0.6)) +
  geom_point(color="black", fill="burlywood1", pch=21) +
  theme_minimal() +
  guides(size=F) +
  xlim(0.95,1.05) +
  scale_size(range=c(2,8)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    panel.grid = element_blank(),
    plot.margin = unit(c(0,0,0,0),"pt"),
    panel.spacing = unit(0, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
#### Ploidy track ####
ploidy_df = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, coding, known_diploid_mitosis, UM_ploidy) %>%
  rename(label = SPECIES.TREE.LABEL, ploidy = coding) %>%
  #mutate(ploidy= ifelse(UM_ploidy == 0, "?", ploidy)) %>%
  mutate(ploidy = as.factor(ploidy), known_diploid_mitosis = as.factor(known_diploid_mitosis))

ploidy_track = tublerone::ggtreeplot(final, ploidy_df, aes(x=1, y=y, fill=ploidy), flip=FALSE, expand_limits = expansion(0.00,0.6)) +
  #geom_point(pch=22, aes(x=1, y=y, fill=known_diploid_mitosis), color="black", size=3.0, inherit.aes = F) +
  geom_point(pch=21, color="black", size=2.5) +
  scale_fill_manual(values = c("lightgrey", "red", "blue")) +
  guides(fill=F)+
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    panel.grid = element_blank(),
    plot.margin = unit(c(0,0,0,0),"pt"),
    panel.spacing = unit(0, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
ploidy_track

#### Character heatmap ####
character_sheet = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, Flagellar.State, `Electron-opaque_plug_in_axoneme_core_and_between_axoneme_and_flagellar_membrane (2)`,
         `Mode (M)`, `Nutritional_Mode (N)`, EFL, EF1a, Cbl_group, MMCoA_group, MetH, Selenocystein,
         Ancestral_combined:Whi5_Fungal) %>%
  rename(opaque_flag_plug = `Electron-opaque_plug_in_axoneme_core_and_between_axoneme_and_flagellar_membrane (2)`) %>%
  mutate_at(vars(Ancestral_combined:Whi5_Fungal), ~ifelse(is.na(.), 0, .))
character_sheet_long = character_sheet %>%
  gather(key="char", value="state", -SPECIES.TREE.LABEL) %>%
  mutate(state = as.factor(state))

pals = list(
  "Flagellar.State" = c("white", "red"),
  "Mode (M)" = c("turquoise2", "green"),
  "Nutritional_Mode (N)" = c("wheat1", "sprintgreen1", "tomato1", "violetred1", "sienna4"),
  "opaque_flag_plug" = c("white", "seagreen4"),
  "EFL" = c("white", "yellow"),
  "EF1a" = c("white", "blue"), #All 1's for now, so no need for white
  "Cbl_group" = c("white", "mediumpurple1", "purple", "magenta4"),
  "MMCoA_group" = c("white", "mediumpurple1", "purple", "magenta4"),
  "MetH" = c("white", "magenta4"),
  "Selenocystein" = c("white", "darkolivegreen"),
  "Ancestral_combined" = c("white", "darkorange", "darkorange3"),
  "Fungal_combined" = c("white", "darkseagreen1", "darkseagreen3"),
  "E2F_Ancestral" = c("white", "darkorange3"),
  "Rb_Ancestral" = c("white", "darkorange3"),
  "SBF_Fungal" = c("white", "darkorange1"),
  "Whi5_Fungal" = c("white", "darkorange1")
)

columns = list()
i=1
for (c in colnames(character_sheet)[-1]) {
  sub = character_sheet_long %>%
    filter(char == c) %>%
    rename(label = SPECIES.TREE.LABEL)
  plt_tbl = final$data %>%
    filter(isTip == T) %>%
    select(label, y) %>%
    left_join(sub)
  print(c)
  print(pals[[c]])
  col = tublerone::ggtreeplot(final, plt_tbl, aes(x=1, y=y, fill=state), flip=FALSE, expand_limits = expansion(0.00,0.6)) +
    geom_tile(color="black") +
    coord_fixed()+
    theme_minimal() +
    guides(fill=F) +
    scale_fill_manual(values=pals[[c]], na.value="lightgrey")+
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

#### SNP density bar charts ####
l50_genome_sizes = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, l50_assembly_length)
isolates = read_xlsx(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Isolates.xlsx"))

snp_densities = read_delim(file.path(CHYTRID_PHYLO, "figures/1_tree", "all.snp_contig_counts_strainified.tsv"), delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(num_snps = sum(num_snps)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  left_join(l50_genome_sizes) %>%
  select(SPECIES.TREE.LABEL, num_snps, l50_assembly_length) %>%
  mutate(snp_density = num_snps/l50_assembly_length) %>%
  rename(label = SPECIES.TREE.LABEL)

#Cut for now: We didn't do ploidy for it
#snp_density_max_x = max(snp_densities$snp_density)
#grey_bars = isolates  %>%
#  filter(!SPECIES.TREE.LABEL %in% snp_densities$label) %>%
#  filter(!SPECIES.TREE.LABEL %in% c("Coelomomyces_lativittatus_CIRM-AVA-1-Amber.LCG",
#                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Orange.LCG",
#                                    "Rozella_rhizoclosmatii")) %>%
#  select(SPECIES.TREE.LABEL) %>%
#  mutate(num_snps =NA, l50_assembly_length=NA, snp_density=NA) %>%
#  rename(label=SPECIES.TREE.LABEL)#
#
#snp_densities = as_tibble(
#  rbind(snp_densities, grey_bars)
#)


snp_density_bar = tublerone::ggtreeplot(final, snp_densities, aes(x=snp_density, y=y), flip=FALSE, expand_limits = expansion(0.00,0.6)) +
  geom_barh(stat = "identity", color="black", fill="darkgreen") +
  scale_x_continuous(expand=c(0,0)) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    panel.border = element_blank()
  )

snp_density_bar
#####
#### After everything is merged and aligned, remove extraneous text from tip labels
#final_cleaned_tiplabs = final
#final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
#  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
#  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
#  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
#  mutate(label = gsub("LCG$", "", label))

#final_cleaned_tiplabs$data[123,4] = "Martensiomyces pterosporus CBS 209.56"
#final_cleaned_tiplabs$data[124,4] = "Ramicandelaber brevisporus CBS 109374"
#final_cleaned_tiplabs$data[119,4] = "Conidiobolus thromboides FSU 785"
#final_cleaned_tiplabs$data[129,4] = "Coelomomyces lativittatus CIRM-AVA-1"
#final_cleaned_tiplabs$data[126,4] = "Syncephalis fuscata S228"
#final_cleaned_tiplabs$data[127,4] = "Olpidium bornovanus UCB F19785"
#final_cleaned_tiplabs$data[32,4] = "Blyttiomyces helicus"
#final_cleaned_tiplabs$data[7,4] = "Paraphelidium tribonemae X-108"



#### Remove strain ids from final tip labels
#final_cleaned_tiplabs = final
#final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
#  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
#  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
#  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
#  mutate(label = gsub("LCG$", "", label)) %>%
#  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
#  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species)))

final_cleaned_tiplabs = final
final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
  mutate(label = gsub("LCG$", "", label)) %>%
  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = final_cleaned_tiplabs$data %>%
  filter(isTip) %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

final_cleaned_tiplabs$data = final_cleaned_tiplabs$data %>%
  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(label = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), label))
#### Put it together
character_heatmap = columns[[1]] + columns[[2]] + columns[[3]] + columns[[4]] + columns[[5]] +
  columns[[6]] + columns[[7]] + columns[[8]] + columns[[9]] + columns[[10]] + columns[[13]] + columns[[14]] + columns[[15]] + columns[[16]] +
  plot_layout(nrow=1)

pwdes = "
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
111111111111111111345555556
222222222222222222#########
"
combined_tree = final_cleaned_tiplabs +
  timescale +
  genome_size_column +
  #ploidy_track +
  character_heatmap +
  #snp_density_bar +
  plot_layout(design=pwdes)
  #plot_layout(nrow=1, widths = c(5.5,0.3,0.3,2,1))

#widths = c(5.5,0.3,0.3,2,1)
combined_tree
ggsave(filename = file.path(CHYTRID_PHYLO, "figures/1_tree/", "Tree-Figure_2_083021.pdf"), plot = combined_tree, device=cairo_pdf, width = 8.5, height = 11, units = "in")



