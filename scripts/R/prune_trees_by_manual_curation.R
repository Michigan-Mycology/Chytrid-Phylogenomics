library(phytools)
library(tidyverse)
library(readxl)
library(ggtree)
library(optparse)

option_list = list(
  make_option(c("-m", "--mancur"), type="character", help = "Path to manual curation sheet. Format is important. Use the example sheet on the Chytrid-Phylogenomics repo."),
  make_option(c("-i", "--isolates"), type="character", help = "Path to isolates sheet with at least two columns named `SPECIES.TREE.LABEL` and `LTP`."),
  make_option(c("-r", "--hitreport"), type="character", help = "Path to `hit_report_all.tsv` from `gtfilter.py`"),
  make_option(c("-g", "--genetrees"), type="character", help = "Path to directory containing all the renamed gene trees."),
  make_option(c("-o", "--outdir"), type="character", help = "Path to directory to write output `*toget` files.")
)
opt_parser = OptionParser(option_list = option_list, add_help_option = T)
opt = parse_args(opt_parser)

mancur = read_excel(opt$mancur) %>%
  rename(marker = `tree number`, delete_tips = tips_names) %>%
  gather(key = "category", value = "value", -delete_tips, -marker) %>%
  filter(!is.na(value)) %>%
  mutate(category = toupper(category)) %>%
  mutate(delete_internal_nodes = NA)

isolates = read_xlsx(opt$isolates) %>%
  select(SPECIES.TREE.LABEL, LTP)

scores = read_delim(opt$hitreport, delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  mutate(tiplab_display = tiplab)


path = opt$genetrees
files = list.files(path = path)

dropped_tips = matrix(nrow=0, ncol=2)
for (f in files) {
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)
  outname = file.path(opt$outdir, paste(this_marker, "_tiplabs", sep = ""))

  print(this_marker)

  this_marker_cat = mancur %>%
    dplyr::filter (marker == this_marker) %>%
    .$category

  print(this_marker_cat)

  if (this_marker_cat == "GOOD") {

    # This marker is good - print tiplabels as is.
    tree = read.newick( paste(path, f, sep="/") )
    write.table(tree$tip.label, outname, quote=FALSE, row.names = FALSE)

  }

  else if (this_marker_cat == "BAD") {

    #This marker is bad - do nothing."
    next

  }

  else if ( this_marker_cat == "SALVAGEABLE" ) {
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

      all_tips_to_drop =
      filtered_tree = tree
      for (c in clean_nodes) {

        # By node number
        tips_to_drop = as_tibble(tree$tip.label[getDescendants(tree, c)]) %>%
          drop_na() %>%
          .$value
        filtered_tree = drop.tip(filtered_tree, tips_to_drop)

        #By edge number
        #tips_to_drop = as_tibble(tree$tip.label[getDescendants(tree, tree$edge[c,2])]) %>%
        #  drop_na() %>%
        #  .$value
        #filtered_tree = drop.tip(filtered_tree, tips_to_drop)

        dropped_tips = rbind(dropped_tips, c(this_marker, toString(tips_to_drop)))
      }

      stopifnot(length(filtered_tree$tip.label[filtered_tree$tip.label %in% tips_to_drop]) == 0)

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

        if (s %in% tree$tip.label) {
          clean_tips = c(clean_tips, s)
        }
      }

      filtered_tips = tree$tip.label[!tree$tip.label %in% clean_tips]
      stopifnot ( length(tree$tip.label) - length(clean_tips) == length(filtered_tips) )

      write.table(filtered_tips, outname, quote=FALSE, row.names = FALSE)
      dropped_tips = rbind(dropped_tips, c(this_marker, toString(clean_tips)))

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
