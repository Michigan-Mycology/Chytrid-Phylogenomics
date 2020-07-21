library(phytools)
library(tidyverse)
library(readxl)
library(ggtree)

mancur = read_excel("~/DATA/phylogeny/gene_tree_curation_round2.xlsx", sheet = "combined") %>%
  mutate(delete_internal_nodes = as.character(delete_internal_nodes))

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, TAX.GROUP, LTP)

scores = read_delim("~/DATA/phylogeny/hit_report_all.csv", delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  mutate(tiplab_display = tiplab)

filt_path = "/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/fast_gene_trees_renamed"
path = filt_path
files = list.files(path = path)

setwd("~/DATA/phylogeny/round8_mancur_of_round7")      
for (f in files) {
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)
  
  this_marker_cat = mancur %>%
    dplyr::filter (marker == this_marker) %>%
    .$category
  
  if (toupper(this_marker_cat) == "GOOD") {
    
    # This marker is good - print tiplabels as is.
    tree = read.newick( paste(path, f, sep="/") )
    outname = paste(this_marker, "_tiplabs", sep = "")
    write.table(tree$tip.label, outname, quote=FALSE, row.names = FALSE)
    
  } 
  
  else if (toupper(this_marker_cat) == "BAD") {
    
    #This marker is bad - do nothing."
    next
  
  } 
  
  else if ( toupper(this_marker_cat) == "SALVAGEABLE" ) {
    # This tip is salvageable, but tips OR nodes need to be deleted.
    
    tree = read.newick( paste(path, f, sep="/") )
    tree = midpoint.root(tree)
    
    this_marker_row = mancur %>% 
      dplyr::filter(marker == this_marker)
    
    # Tips and Nodes both have entries
    if ( !is.na(this_marker_row$delete_internal_nodes) & !is.na(this_marker_row$delete_tips) ) {
      print(paste(this_marker, "There are entries in both columns."))
    }
    
    # Only Nodes have an entry
    else if ( !is.na(this_marker_row$delete_internal_nodes) & is.na(this_marker_row$delete_tips) ) {
      nodes_to_delete = this_marker_row %>%
        .$delete_internal_nodes
      spl = unlist(strsplit(nodes_to_delete, split=";"))
      clean_nodes = c()
      for (s in spl) {
        s = gsub("^[ ]", "", s)
        clean_nodes = c(clean_nodes,s)
      }
      clean_nodes = as.numeric(clean_nodes)
      
      filtered_tree = tree
      for (c in clean_nodes) {
        
        # By node number
        tips_to_drop = as_tibble(tree$tip.label[getDescendants(tree, c)]) %>%
          drop_na() %>%
          .$value
        filtered_tree = drop.tip(filtered_tree, tips_to_drop)
        
        # By edge number
        #tips_to_drop = as_tibble(tree$tip.label[getDescendants(tree, tree$edge[c,2])]) %>%
          #drop_na() %>%
          #.$value
        #filtered_tree = drop.tip(filtered_tree, tips_to_drop)
      }

      stopifnot(length(filtered_tree$tip.label[filtered_tree$tip.label %in% tips_to_drop]) == 0)
      
      outname = paste(this_marker, "_tiplabs", sep = "")
      write.table(filtered_tree$tip.label, outname, quote=FALSE, row.names = FALSE)
      
      d = scores %>%
        dplyr::filter(marker == this_marker)
      ntips = length(filtered_tree$tip.label)
        
      annotated_filtered_tree = ggtree(filtered_tree, branch.length = 0.1) %<+% d +
        #geom_label(aes(x=branch, label=edge_num), cex=2, label.padding=unit(0.05, "cm")) + 
        geom_tiplab(aes(color=score), cex=3) +
        scale_color_continuous(type="viridis") +
        ggtitle(label=this_marker) +
        theme_tree()
        
        x_max = layer_scales(annotated_filtered_tree)$x$range$range[2]
        annotated_filtered_tree = annotated_filtered_tree +
          ggplot2::xlim(0, x_max*1.5)
        
        ggsave(filename = paste(f,".tree_with_clade_rm", ".pdf", sep=""), plot = annotated_filtered_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
      }
    
    # Only Tips have an entry
    else if ( is.na(this_marker_row$delete_internal_nodes) & !is.na(this_marker_row$delete_tips) ) {
      tips_to_delete = this_marker_row %>%
        .$delete_tips
      spl = unlist(strsplit(tips_to_delete, split=";"))
      clean_tips = c()
      for (s in spl) {
        s = gsub("−","-", s)
        s = gsub("[ ][-][ ][0-9.]*", "", s)
        s = gsub("^[ ]", "", s)
        clean_tips = c(clean_tips, s)
      }
      filtered_tips = tree$tip.label[!tree$tip.label %in% clean_tips]
      stopifnot ( length(tree$tip.label) - length(spl) == length(filtered_tips) )
      
      outname = paste(this_marker, "_tiplabs", sep = "")
      write.table(filtered_tips, outname, quote=FALSE, row.names = FALSE)
      
    }
    
    # Neither Nodes or Tips have entries
    else if ( is.na(this_marker_row$delete_internal_nodes) & is.na(this_marker_row$delete_tips) ) {
      print(paste(this_marker, "No entries"))
    }
    
    # This shouldn't happen as long as all conditions are covered.
    else {
      print(paste(this_marker, "Missed a condition."))
    }
  }
  else {
    
    # This shouldn't happen. Print something out if it does.
    print(paste(this_marker, "has a bad category name", this_marker_cat))
    
  }
}
