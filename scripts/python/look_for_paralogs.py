#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 18:10:43 2020

@author: aimzez
"""

#%%
import sys
import os
from ete3 import Tree
import pandas as pd
import numpy as np

#%% Object to handle calculations and matricies
class MultiMarkerGeneTrees(object):
    def __init__ (self, tree_paths, suffix = ".aa.tre.renamed"):
        self.trees = {}
        tree_files = tree_paths
        for t in tree_files:
            this_marker = os.path.basename(t).replace(suffix,"")
            tree = Tree(t)
            for tip in [n for n in tree.get_leaves()]:
                spl = tip.name.split("&")
                spl_2 = spl[1].split("|")
                tree_name = spl[0]
                isolate = spl_2[0]
                phead = spl_2[1]
                tip.add_feature("isolate", isolate)
                tip.add_feature("gene", phead)
                tip.add_feature("genus", tip.name.split("_")[0])
            self.trees[this_marker] = tree
        
    
    def annotate_scores (self, metadata):
        for marker,tree in self.trees.items():
            sub = metadata[(metadata.marker == marker)].groupby("isolate")
            for isolate,group in sub:
                for tup in group.itertuples():
                    leaf = [l for l in tree.get_leaves() if l.gene == tup.gene and l.isolate == tup.isolate]
                    if(len(leaf) == 0): 
                        print(f"Missing leafs: {tup}")
                        continue
                    elif(len(leaf) > 1):
                        print(f"Multiple leafs: {tup}")
                        sys.exit()
                        continue
                    else:
                        leaf = leaf[0]
                        leaf.add_feature("score", tup.score)
                        leaf.add_feature("evalue", tup.evalue)
    
    def filter_pipeline (self):
        pass
            
    
    def max_score_difference(self):
        ddict = {}
        for marker,tree in self.trees.items():
            u_isolates = set([x.isolate for x in tree.get_leaves()])
            d = {k: v for k,v in zip(u_isolates, [0]*len(u_isolates))}
            
            for i in u_isolates:
                i_scores = [l.score for l in tree.get_leaves() if l.isolate == i]
                d[i] = max(i_scores) - min(i_scores)
                
            ddict[marker] = d
        
        return ddict

    def monophyly_matrix(self):
        ddict = {}
        for marker,tree in self.trees.items():
            u_isolates = set([x.isolate for x in tree.get_leaves()])
            d = {k: v for k,v in zip(u_isolates, [0]*len(u_isolates))}
            
            for i in u_isolates:
                b, r, _ = tree.check_monophyly([i], "isolate", unrooted=True)
                if b:
                    d[i] = 1
            
            ddict[marker] = d
        
        return ddict
    
#%%
def filter_pipeline(trees, outdir, remove_remaining_polyphyletic_taxa = False):
    for marker,tree in trees.trees.items():
        #print(f"------{marker}")
        u_isolates = u_isolates = set([x.isolate for x in tree.get_leaves()])
        filtered_d = {}
        for i in u_isolates:
            high_name = ""
            low_name = ""
            high_names = []
            low_names = []
            high_idx = None
            low_idx = None
            
            
            i_names = [x.name for x in tree.get_leaves() if x.isolate == i]
            i_scores = [x.score for x in tree.get_leaves() if x.isolate == i]
            i_genus = [x.genus for x in tree.get_leaves() if x.isolate == i][0]
            
            high_idx = i_scores.index(max(i_scores))
            high_name = i_names[high_idx]
            
            #print(i, i_names)
            if len(i_scores) == 1:
                
                filtered_d[i] = [high_name]
                continue #There is only one hit. Take it.
            
            else: #There is more than one hit. Compare scores.
                bool_within_30p = [x >= max(i_scores)*0.70 for x in i_scores]
                if all(bool_within_30p):
                    
                    b,r,_ = tree.check_monophyly([i], "isolate", unrooted=True)
                    if b:
                        filtered_d[i] = [high_name]
                        continue
                    else:
                        if genus_is_monophyletic(i_genus, tree):
                            filtered_d[i] = [high_name]
                            continue
                        else:
                            # Either remove this isolate from the tree ...
                            if remove_remaining_polyphyletic_taxa:
                                continue
                            
                            # ... Or put all the remaining hits back in the tree
                            else:
                                filtered_d[i] = i_names
                                continue #All the hit scores are within 30% of the maximum score. Ask monophyly. If yes, take highest. If no, look @ gene tree.
                else:
                    not_within_30p = [x for x in bool_within_30p if not x]
                    if len(not_within_30p) == 1:
                        low_idx = i_scores.index(min(i_scores))
                        low_name = i_names[low_idx]
                        tree_copy = tree.copy()
                        tree_copy.get_leaves_by_name(low_name)[0].delete()
                        #print(f"One tip with low score, {low_name}. Deleting it and reasking monophyly.")
                        b,r,_ = tree_copy.check_monophyly([i], "isolate", unrooted=True)
                            
                        if b:
                            #print(f"{i} is monophyletic now. Let's take the best hit: {high_name}\n")
                            filtered_d[i] = [high_name]
                            continue #There is only one low hit, and removing it DOES make the isolate monophyletic. Take the highest hit.
                        else:
                            if genus_is_monophyletic(i_genus, tree_copy):
                                filtered_d[i] = [high_name]
                                continue
                            else:
                                #print(f"{i} is still not monophyletic. Remove the low one, put all others in the tree again.")
                                i_names.pop(low_idx)
                                # Either remove this isolate from the tree ...
                                if remove_remaining_polyphyletic_taxa:
                                    continue
                                
                                # ... Or put all the remaining hits back in the tree
                                else:
                                    filtered_d[i] = i_names
                                    continue #There is only one low hit, but removing it DOESN'T make the isolate monophyletic. Keep all high hits in tree.
                    else:
                        is_within_30p = [x for x in bool_within_30p if x]
                        if len(is_within_30p) == 1:
                            #print(f"{i}: singular high hit, {high_name}")
                            filtered_d[i] = [high_name]
                            continue #There is only one high hit. Drop lows and take it. 
                        
                        else: #There are multiple high hits. Drop all lows and reask monophyly. If yes, take highest. If no, look @ gene tree.
                            high_indices = [i for i in range(0,len(i_scores)) if bool_within_30p[i]]
                            high_names = [i_names[i] for i in high_indices]
                            
                            low_indices = [i for i in range(0,len(i_scores)) if not bool_within_30p[i]]
                            low_names = [i_names[i] for i in low_indices]
                            tree_copy = tree.copy()
                            for leaf in low_names:
                                tree_copy.get_leaves_by_name(leaf)[0].delete()
                            
                            b,r,_ = tree_copy.check_monophyly([i], "isolate", unrooted=True)
                            
                            if b:
                                #print(f"{i} is monophyletic now. Let's take the best hit: {high_name}\n")
                                filtered_d[i] = [high_name]
                                continue # There are multiple high hits, but removing all the low ones makes the isolate monophyletic. Take it!
                            else:
                                if genus_is_monophyletic(i_genus, tree_copy):
                                    filtered_d[i] = [high_name]
                                    continue
                                
                                # Isolate is STILL polyphyletic
                                else:
                                    # Either remove this isolate from the tree ...
                                    if remove_remaining_polyphyletic_taxa:
                                        continue
                                    
                                    # ... Or put all the remaining hits back in the tree
                                    else:
                                        filtered_d[i] = high_names
                                        continue # There are multiple high hits, and removing all the low ones DOESN'T make the isolate monophytletic. Keep all the high hits in the tree.
            
            print (f"ERROR: Iteration for {i} in {marker} did nothing.")
            sys.exit(1)
          
        with open(os.path.join(outdir, f"{marker}.toget"), 'w') as out:
            for key,values in filtered_d.items():
                 #print(key, values)
                out.write("\n".join([x.split("&")[1] for x in values]))
                out.write("\n")
    return filtered_d
#%%
def final_score_filter (trees, outpath):
    for marker, tree in trees.trees.items():
        scores = [x.score for x in tree.get_leaves()]
        mean = np.mean(scores)
        std = np.std(scores)
    
        tree_copy = tree.copy()
        for leaf in [x.name for x in tree.get_leaves() if x.score <= mean-(std*1.5)]:
            print(leaf)
            tree_copy.get_leaves_by_name(leaf)[0].delete()
        
        with open(os.path.join(outpath, f"{marker}.toget"), 'w') as out:
            for leaf in tree_copy.get_leaves():
                out.write(leaf.name.split("&")[1])
                out.write("\n")
            
        
#%%
def genus_is_monophyletic(genus, tree):
    b,r,_ = tree.check_monophyly([genus], "genus", unrooted=True)
    return b
                      
#%%
if __name__ == "__main__":
    #Objectified so far
    #directory = "/home/aimzez/DATA/phylogeny/round1_raw/fast_gene_trees_renamed" ### THIS POINTED TO FILTERED TREES - REDO TO CONFIRM?
    directory = "/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/fast_gene_trees_renamed"
    files = [os.path.join(directory,t) for t in os.listdir(directory) if t.endswith("renamed")]
    markers = [os.path.basename(x.replace(".aa.tre.renamed", "")) for x in files]
    
    metadata = pd.read_csv("~/DATA/phylogeny/hit_report_all.csv", sep="\t", header=None)
    metadata.columns = ["gene", "marker", "evalue", "score"]
    new = metadata.gene.str.split("|", expand=True)
    metadata["isolate"] = new[0]
    metadata["gene"] = new[1]
    
    all_trees = MultiMarkerGeneTrees(files)
    #%% Post-Trimal, Pre-Filter Monophyly Matrix
    mono_ddict = all_trees.monophyly_matrix()
    frame = pd.DataFrame(mono_ddict).fillna(-1)
    frame.to_csv("/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/monophyly_matrix.csv", sep=",")
    #%%
    all_trees.annotate_scores(metadata)
    
    #%%
    fd = filter_pipeline(
        all_trees, 
        outdir = "/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/toget", 
        remove_remaining_polyphyletic_taxa = True
        )
    
    #%% Post-Trimal, Post-Filter Monophyly Matrix
    directory = "/home/aimzez/DATA/phylogeny/round4_scorefilt_monogenus/fast_gene_trees_renamed"
    files = [os.path.join(directory,t) for t in os.listdir(directory) if t.endswith("renamed")]
    markers = [os.path.basename(x.replace(".aa.tre.renamed", "")) for x in files]
    all_trees = MultiMarkerGeneTrees(files)
    mono_ddict = all_trees.monophyly_matrix()
    frame = pd.DataFrame(mono_ddict).fillna(-1)
    frame.to_csv("/home/aimzez/DATA/phylogeny/post_trimal_post_filter_monophyly_matrix.csv", sep=",")

 #%% Post-Trimal, Post-Filter, Post-SCOREFILT Monophyly Matrix
    directory = "/home/aimzez/DATA/phylogeny/round5_final_score_filt/fast_gene_trees_renamed"
    files = [os.path.join(directory,t) for t in os.listdir(directory) if t.endswith("renamed")]
    markers = [os.path.basename(x.replace(".aa.tre.renamed", "")) for x in files]
    all_trees = MultiMarkerGeneTrees(files)
    mono_ddict = all_trees.monophyly_matrix()
    frame = pd.DataFrame(mono_ddict).fillna(-1)
    frame.to_csv("/home/aimzez/DATA/phylogeny/post_trimal_post_filter_final_score_filt_monophyly_matrix.csv", sep=",")
    
    #%%
    ddict = all_trees.max_score_difference()
    frame = pd.DataFrame(ddict)
    frame.to_csv("/home/aimzez/DATA/phylogeny/max_score_difference.csv", sep=",")

'''
#%% Generate monophyly matrix 
directory = "/home/aimzez/DATA/phylogeny/fast_gene_trees/"
files = [os.path.join(directory,t) for t in os.listdir(directory) if t.endswith("renamed")]
markers = [os.path.basename(x.replace(".aa.tre.renamed", "")) for x in files]

ddict = {}

for t in files:
    this_marker = os.path.basename(t).replace(".aa.tre.renamed","")
    tree = Tree(t)
    for tip in [n for n in tree.get_leaves()]:
        spl = tip.name.split("&")
        spl_2 = spl[1].split("|")
        tree_name = spl[0]
        isolate = spl_2[0]
        phead = spl_2[1]
        tip.add_feature("isolate", isolate)
        tip.add_feature("gene", phead)
    
    u_isolates = set([x.isolate for x in tree.get_leaves()])
    d = {k: v for k,v in zip(u_isolates, [0]*len(u_isolates))}

    for i in u_isolates:
        b, r, _ = tree.check_monophyly([i], "isolate", unrooted=True)
        if b:
            d[i] = 1
    
    ddict[this_marker] = d

matrix = pd.DataFrame(ddict).fillna(-1)
matrix.to_csv("/home/aimzez/DATA/phylogeny/check_monophly.csv", sep=",")
#%% Calculate distance of some kind of distance between paralogs
ddict = {}
for t in files:
    this_marker = os.path.basename(t).replace(".aa.tre.renamed","")
    tree = Tree(t)
    #l.append(len(tree.get_leaves()))
    #continue

    for tip in [n for n in tree.get_leaves()]:
        spl = tip.name.split("&")
        spl_2 = spl[1].split("|")
        tree_name = spl[0]
        isolate = spl_2[0]
        phead = spl_2[1]
        tip.add_feature("isolate", isolate)
        tip.add_feature("gene", phead)
    
    u_isolates = set([x.isolate for x in tree.get_leaves()])
    d = {k: v for k,v in zip(u_isolates, [0]*len(u_isolates))}
    for i in u_isolates:
        tips_to_distance = []
        for idx,tip in enumerate(tree.get_leaves()):
            if tip.isolate == i:
                tips_to_distance.append(tip)
        
        if len(tips_to_distance) == 1:
            continue
        first_tip = tips_to_distance[0]
        tips_to_distance = tuple(i for i in tips_to_distance[1:])
        d[i] = len(first_tip.get_common_ancestor(*tips_to_distance).get_leaves())/len(tree.get_leaves())
    
    ddict[this_marker] = d
    matrix = pd.DataFrame(ddict)
    matrix.to_csv("/home/aimzez/DATA/phylogeny/leaves_of_common_ancestor.csv", sep=",")
    
#%% Dupication level
#%% Calculate distance of some kind of distance between paralogs
ddict = {}
for t in files:
    this_marker = os.path.basename(t).replace(".aa.tre.renamed","")
    tree = Tree(t)
    #l.append(len(tree.get_leaves()))
    #continue

    for tip in [n for n in tree.get_leaves()]:
        spl = tip.name.split("&")
        spl_2 = spl[1].split("|")
        tree_name = spl[0]
        isolate = spl_2[0]
        phead = spl_2[1]
        tip.add_feature("isolate", isolate)
        tip.add_feature("gene", phead)
    
    u_isolates = set([x.isolate for x in tree.get_leaves()])
    d = {k: v for k,v in zip(u_isolates, [0]*len(u_isolates))}
    for i in u_isolates:
        tips_to_distance = []
        for idx,tip in enumerate(tree.get_leaves()):
            if tip.isolate == i:
                tips_to_distance.append(tip)
        if len(tips_to_distance) == 1:
            d[i] = np.nan
        else:
            d[i] = len(tips_to_distance)
    
    ddict[this_marker] = d
    matrix = pd.DataFrame(ddict).fillna(0)
    matrix.to_csv("/home/aimzez/DATA/phylogeny/duplication_level.csv", sep=",")
'''