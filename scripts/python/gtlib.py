#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:14:22 2020

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

    def __init__(self, tree_paths, rooted = False, suffix=".aa.tre.renamed"):
        self.trees = {}
        tree_files = tree_paths
        for t in tree_files:
            this_marker = os.path.basename(t).replace(suffix, "")
            if not rooted:
                tree = Tree(t)
            else:
                tree = Tree(t, format=1)
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
    
    def midpoint_root(self):
        for m,t in self.trees.items():
            t.set_outgroup(t.get_midpoint_outgroup())
        return None

    def outgroup_root(self, outgroup_isolates):
        new_trees = {}
        for m,t in self.trees.items():
            print(m)
            outgroup = [x for x in t.get_leaves() if x.isolate in outgroup_isolates]
            if len(outgroup) == 0:
                new_trees[m] =t
                continue
            root_node = t.get_common_ancestor(outgroup)
            if root_node is t:
                print("Root is whole tree.")
                new_trees[m] = t
                continue
            t.set_outgroup(root_node)
            new_trees[m] = t

        self.trees = new_trees
        return None


    def annotate_scores(self, metadata):
        for marker, tree in self.trees.items():
            sub = metadata[(metadata.marker == marker)].groupby("isolate")
            for isolate, group in sub:
                for tup in group.itertuples():
                    leaf = [l for l in tree.get_leaves() if l.gene ==
                            tup.gene and l.isolate == tup.isolate]
                    if(len(leaf) == 0):
                        print(f"Leaf not in tree (this is OK): {tup}")
                        continue
                    elif(len(leaf) > 1):
                        print(f"Multiple leafs corresponding to this hit (this not OK): {tup}")
                        sys.exit()
                        continue
                    else:
                        leaf = leaf[0]
                        leaf.add_feature("score", tup.score)
                        leaf.add_feature("evalue", tup.evalue)

    def max_score_difference(self):
        ddict = {}
        for marker, tree in self.trees.items():
            u_isolates = set([x.isolate for x in tree.get_leaves()])
            d = {k: v for k, v in zip(u_isolates, [0] * len(u_isolates))}

            for i in u_isolates:
                i_scores = [l.score for l in tree.get_leaves()
                            if l.isolate == i]
                d[i] = max(i_scores) - min(i_scores)

            ddict[marker] = d

        return ddict

    def monophyly_matrix(self):
        ddict = {}
        for marker, tree in self.trees.items():
            u_isolates = set([x.isolate for x in tree.get_leaves()])
            d = {k: v for k, v in zip(u_isolates, [0] * len(u_isolates))}

            for i in u_isolates:
                b, r, _ = tree.check_monophyly([i], "isolate", unrooted=True)
                if b:
                    d[i] = 1

            ddict[marker] = d

        return ddict

#%%


def filter_pipeline(trees, outdir, remove_remaining_polyphyletic_taxa=False):
    for marker, tree in trees.trees.items():
        # print(f"------{marker}")
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
                continue  # There is only one hit. Take it.

            else:  # There is more than one hit. Compare scores.
                bool_within_30p = [x >= max(i_scores) * 0.70 for x in i_scores]
                if all(bool_within_30p):

                    b, r, _ = tree.check_monophyly(
                        [i], "isolate", unrooted=True)
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
                                # All the hit scores are within 30% of the
                                # maximum score. Ask monophyly. If yes, take
                                # highest. If no, look @ gene tree.
                                continue
                else:
                    not_within_30p = [x for x in bool_within_30p if not x]
                    if len(not_within_30p) == 1:
                        low_idx = i_scores.index(min(i_scores))
                        low_name = i_names[low_idx]
                        tree_copy = tree.copy()
                        tree_copy.get_leaves_by_name(low_name)[0].delete()
                        # print(f"One tip with low score, {low_name}. Deleting
                        # it and reasking monophyly.")
                        b, r, _ = tree_copy.check_monophyly(
                            [i], "isolate", unrooted=True)

                        if b:
                            # print(f"{i} is monophyletic now. Let's take the
                            # best hit: {high_name}\n")
                            filtered_d[i] = [high_name]
                            # There is only one low hit, and removing it DOES
                            # make the isolate monophyletic. Take the highest
                            # hit.
                            continue
                        else:
                            if genus_is_monophyletic(i_genus, tree_copy):
                                filtered_d[i] = [high_name]
                                continue
                            else:
                                # print(f"{i} is still not monophyletic. Remove
                                # the low one, put all others in the tree
                                # again.")
                                i_names.pop(low_idx)
                                # Either remove this isolate from the tree ...
                                if remove_remaining_polyphyletic_taxa:
                                    continue

                                # ... Or put all the remaining hits back in the tree
                                else:
                                    filtered_d[i] = i_names
                                    # There is only one low hit, but removing
                                    # it DOESN'T make the isolate monophyletic.
                                    # Keep all high hits in tree.
                                    continue
                    else:
                        is_within_30p = [x for x in bool_within_30p if x]
                        if len(is_within_30p) == 1:
                            # print(f"{i}: singular high hit, {high_name}")
                            filtered_d[i] = [high_name]
                            # There is only one high hit. Drop lows and take
                            # it.
                            continue

                        else:  # There are multiple high hits. Drop all lows and reask monophyly. If yes, take highest. If no, look @ gene tree.
                            high_indices = [i for i in range(
                                0, len(i_scores)) if bool_within_30p[i]]
                            high_names = [i_names[i] for i in high_indices]

                            low_indices = [i for i in range(
                                0, len(i_scores)) if not bool_within_30p[i]]
                            low_names = [i_names[i] for i in low_indices]
                            tree_copy = tree.copy()
                            for leaf in low_names:
                                tree_copy.get_leaves_by_name(leaf)[0].delete()

                            b, r, _ = tree_copy.check_monophyly(
                                [i], "isolate", unrooted=True)

                            if b:
                                # print(f"{i} is monophyletic now. Let's take
                                # the best hit: {high_name}\n")
                                filtered_d[i] = [high_name]
                                continue  # There are multiple high hits, but removing all the low ones makes the isolate monophyletic. Take it!
                            else:
                                if genus_is_monophyletic(i_genus, tree_copy):
                                    filtered_d[i] = [high_name]
                                    continue

                                # Isolate is STILL polyphyletic
                                else:
                                    # Either remove this isolate from the tree
                                    # ...
                                    if remove_remaining_polyphyletic_taxa:
                                        continue

                                    # ... Or put all the remaining hits back in the tree
                                    else:
                                        filtered_d[i] = high_names
                                        # There are multiple high hits, and
                                        # removing all the low ones DOESN'T
                                        # make the isolate monophytletic. Keep
                                        # all the high hits in the tree.
                                        continue

            print(f"ERROR: Iteration for {i} in {marker} did nothing.")
            sys.exit(1)

        with open(os.path.join(outdir, f"{marker}.toget"), 'w') as out:
            for key, values in filtered_d.items():
                 #print(key, values)
                out.write("\n".join([x.split("&")[1] for x in values]))
                out.write("\n")
    return filtered_d
#%%


def final_score_filter(trees, outpath):
    for marker, tree in trees.trees.items():
        scores = [x.score for x in tree.get_leaves()]
        mean = np.mean(scores)
        std = np.std(scores)

        tree_copy = tree.copy()
        for leaf in [x.name for x in tree.get_leaves() if x.score <= mean - (std * 1.5)]:
            print(leaf)
            tree_copy.get_leaves_by_name(leaf)[0].delete()

        with open(os.path.join(outpath, f"{marker}.toget"), 'w') as out:
            for leaf in tree_copy.get_leaves():
                out.write(leaf.name.split("&")[1])
                out.write("\n")


#%%
def genus_is_monophyletic(genus, tree):
    b, r, _ = tree.check_monophyly([genus], "genus", unrooted=True)
    return b
