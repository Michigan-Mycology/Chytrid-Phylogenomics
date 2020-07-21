#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:28:51 2020

@author: aimzez
"""
#%%
import sys
import pandas as pd
import numpy as np
sys.path.append("/home/aimzez/work/Chytrid-Phylogenomics/scripts/python")
from look_for_paralogs import MultiMarkerGeneTrees, final_score_filter

#%%
trees = MultiMarkerGeneTrees("/home/aimzez/DATA/phylogeny/round4_scorefilt_monogenus/fast_gene_trees_renamed")

#%%
metadata = pd.read_csv("~/DATA/phylogeny/hit_report_all.csv", sep="\t", header=None)
metadata.columns = ["gene", "marker", "evalue", "score"]
new = metadata.gene.str.split("|", expand=True)
metadata["isolate"] = new[0]
metadata["gene"] = new[1]
    
trees.annotate_scores(metadata)
#%%
final_score_filter(trees, "/home/aimzez/DATA/phylogeny/round4_scorefilt_monogenus/final_score_filt_unaln/")
