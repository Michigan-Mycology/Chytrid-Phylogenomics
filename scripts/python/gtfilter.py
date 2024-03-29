#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:33:50 2020

@author: aimzez
"""

import sys
import os
import argparse
import pandas as pd
import itertools
import copy
import json
import logging
import numpy as np

# os.environ["CHYTRID_PHYLO_PY"] =
# "/home/aimzez/work/Chytrid-Phylogenomics/scripts/python"
sys.path.append(
    os.path.join(os.environ.get('CHYTRID_PHYLO'), "scripts", "python"))
import gtlib
import ete3


def write_monophyly_matrix(all_trees, outpath) -> None:
    mono_ddict = all_trees.monophyly_matrix()
    frame = pd.DataFrame(mono_ddict).fillna(-1)
    frame.to_csv(os.path.join(outpath, "monophyly_matrix.csv"), sep=",")


''' Codes for compare_monophyly_matrices
11 = Monophyletic then, Monophyletic now
10 = Monophyletic then, Polyphyletic now
01 = Polyphyletic then, Monophyletic now
00 = Polyphyletic then, Polyphyletic now
-10 = Missing then, Polyphyletic now
-11 = Missing then, Monophyletic now
-1-1 = Missing then, Missing now
0-1 = Polphyletic then, Missing now
1-1 = Monophyletic then, Missing Now
'''


def compare_monophyly_matrices(all_trees, outpath, phylymat) -> None:
    this_ddict = all_trees.monophyly_matrix()
    this_frame = pd.DataFrame(this_ddict).fillna(-1)
    that_frame = pd.read_csv(phylymat, header=0).set_index("Unnamed: 0")

    ddict = {}
    for t in this_frame.iterrows():
        isolate = t[0]
        d = {}
        for marker, this_value in t[1].items():
            that_value = that_frame.loc[isolate, marker]
            d[marker] = f"{str(int(that_value))}{str(int(this_value))}"
        ddict[isolate] = d

    compare_matrix = pd.DataFrame(ddict).transpose()
    compare_matrix.to_csv(os.path.join(outpath, "compare_matrix.csv"), sep=",")


def read_metadata(hitreport_path) -> pd.DataFrame:
    metadata = pd.read_csv(hitreport_path, sep="\t", header=None)
    metadata.columns = ["gene", "marker", "evalue", "score"]
    new = metadata.gene.str.split("|", expand=True)

    if new.shape[1] != 2:
        print(f"[ERROR] Malformed hit report at: {hitreport_path}")
        print(
            "[ERROR] Make sure that the first column of your hit report follows the format: LocusTagPrefix|Protein_ID"
        )
        sys.exit(1)

    metadata["isolate"] = new[0]
    metadata["gene"] = new[1]
    return metadata


def phyly_score_filter(all_trees, outpath, hitreport_path,
                       remove_remaining_polyphyletic_taxa):

    meta = read_metadata(hitreport_path)
    all_trees.annotate_scores(meta)
    gtlib.filter_pipeline(all_trees, outpath,
                          remove_remaining_polyphyletic_taxa)


def hard_score_filter(all_trees, outpath, hitreport_path):
    meta = read_metadata(hitreport_path)
    all_trees.annotate_scores(meta)
    gtlib.final_score_filter(all_trees, outpath)


def taxon_occupancy(all_trees, isolates, outpath):

    def ltp2stl(ltp, mapper):
        return mapper[ltp]

    ltp_mapper = {}
    with open(isolates, 'r') as iso:
        for line in iso:
            spl = [x.strip() for x in line.split("\t")]
            ltp_mapper[spl[1]] = spl[0]
    mono_ddict = all_trees.monophyly_matrix()
    frame = pd.DataFrame(mono_ddict).replace(0.0, 1.0).fillna(0.0)

    summed = frame.agg("sum", axis="columns") / frame.shape[1]
    summed = summed.sort_values(ascending=False)

    summed_renamed = pd.DataFrame(summed.reset_index())
    summed_renamed.columns = ["ltp", "occupancy"]
    summed_renamed["stl"] = summed_renamed.ltp.apply(ltp2stl,
                                                     mapper=ltp_mapper)
    summed_renamed = summed_renamed[["stl", "ltp", "occupancy"]]

    summed_renamed.to_csv(os.path.join(outpath, "python_isolate_repr.tsv"),
                          sep="\t",
                          index=False)


def write_per_marker_quartet_group_monophyletic_clades(all_trees,
                                                       quartets_path,
                                                       json_path):
    # outgroup_isolates = ["GCF_000016345.1", "GCF_003065365.1", "GCF_000019745.1", "GCF_003261295.1"]
    # all_trees.root(outgroup_isolates)
    all_trees.midpoint_root()

    # for mark, tree in all_trees.trees.items():
    #    for l in tree.get_leaves():
    #        print(l.isolate)

    quartet_groups = {}
    # Read in quartet groups
    with open(quartets_path, 'r') as qf:

        for line in qf:
            spl = [x.strip() for x in line.split("\t")]
            for i in range(0, len(spl)):
                if i not in quartet_groups:
                    quartet_groups[i] = []
                if len(spl[i]) > 0:
                    quartet_groups[i].append(spl[i])

    mono = [0] * (len(quartet_groups) - 1)
    all_mono = 0
    per_marker_quartet_group_monophyletic_clades = {}
    for marker, tree in all_trees.trees.items():
        per_marker_quartet_group_monophyletic_clades[marker] = {}
        copy_qg = copy.deepcopy(quartet_groups)
        # print(copy_qg)
        mono_log = [True] * (len(quartet_groups) - 1)
        for i in range(0, len(quartet_groups) - 1):

            # Skip quartet group if none of the tips are present in the gene tree
            # ^ Necessary because TreeNode.check_monophyly(ignore_missing = True) will return True
            # if all tips are missing.
            if all([
                    x not in [l.isolate for l in tree.get_leaves()]
                    for x in copy_qg[i]
            ]):
                per_marker_quartet_group_monophyletic_clades[marker][i] = None
                continue

            else:
                per_marker_quartet_group_monophyletic_clades[marker][i] = []
                res = tree.check_monophyly(copy_qg[i],
                                           "isolate",
                                           ignore_missing=True,
                                           unrooted=True)
                if res[0]:
                    mono[i] += 1
                    per_marker_quartet_group_monophyletic_clades[marker][
                        i].append(tuple(quartet_groups[i]))
                else:
                    p = 1
                    while len(copy_qg[i]) != 0:
                        for comb in itertools.combinations(
                                copy_qg[i],
                                len(copy_qg[i]) - p):
                            # print(len(copy_qg[i]), len(comb), p)
                            subtract_res = tree.check_monophyly(
                                comb,
                                "isolate",
                                ignore_missing=True,
                                unrooted=True)
                            if subtract_res[0]:
                                # print(marker, i, p, len(comb), len(
                                #    copy_qg[i]), comb)
                                per_marker_quartet_group_monophyletic_clades[
                                    marker][i].append(comb)
                                p -= len(comb) + 1
                                for c in comb:
                                    copy_qg[i].remove(c)
                                break
                        p += 1
                    '''
                    print(list(itertools.combinations(
                        quartet_groups[i], len(quartet_groups[i]) - 1)))
                    print(quartet_groups[i])
                    print(i, [x.isolate for x in res[2]])
                    print(i, [x.isolate for x in res[2]
                              if x.isolate in quartet_groups[i]])
                    '''
                    mono_log[i] = False
        if all(mono_log):
            all_mono += 1

    with open(json_path, "w") as outfile:
        json.dump(per_marker_quartet_group_monophyletic_clades, outfile)

    # for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
    #    for i, clades in qg.items():
    #        print(marker, i, clades)

    print(mono)
    print(all_mono)

    return None


def print_nice_quartet_monophyly(json_path):
    with open(json_path, 'r') as openfile:
        per_marker_quartet_group_monophyletic_clades = json.load(openfile)
        for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
            for i, clades in qg.items():
                print(marker, i, clades)


class qNode(object):

    def __init__(self, tips):
        if isinstance(tips, qNode) or isinstance(tips, list):
            self.tips = tips
        elif isinstance(tips, int):
            self.tips = [tips]
        else:
            raise ValueError(
                "Argument `tips` to qNode is not a list, qNode, or an int")

        self.n = 0

    def __repr__(self):
        if len(self.tips) > 1:
            return str(self.tips).replace("[", "(").replace("]", ")")
        elif len(self.tips) == 1:
            return str(self.tips[0])
        else:
            raise ValueError("Tips attribute of qNode is zero-length.")

    def __iter__(self):
        self.n = 0
        return self

    def __next__(self):
        if self.n <= len(self.tips) - 1:
            res = self.tips[self.n]
            self.n += 1
            return res
        else:
            raise StopIteration

    def to_list(self):
        ltopo = []

        # Terminal qNode objects
        # where len(qNode.tips) == 1
        # and type(qNode.tips[0]) == int
        # will get caught here.
        # This block ensures that they return a nonempty list.
        if len(self.tips) == 1:
            assert isinstance(
                self.tips[0], int
            ), "1-length qNode.tips is not composed of a single integer value"
            ltopo.append(self)
        else:
            assert len(self.tips) == 2, "qNode.tips is longer than 2 elements"
            for t in self.tips:
                assert isinstance(
                    t, qNode
                ), "2-length qNode.tips is of mixed composition when it should be 2 qNode objects"

                if len(t.tips) == 1:
                    ltopo.append(t)
                else:
                    ltopo.append(t.to_list())
        # ltopo.sort()
        return ltopo

    def get_subtopologies_as_strings(self):
        lres = []
        lrepr = self.to_list()
        # print(lrepr)
        # return None
        while True:
            next_tips = []
            for t in lrepr:
                # print("T", t)
                # print("LRES", lres)
                if isinstance(t, list):
                    for tt in t:
                        # print("TT", tt)
                        # if tt not in lres:
                        next_tips.append(tt)
                    add_t = [str(x) for x in t]
                    add_t.sort()
                    lres.append(add_t)
                else:
                    lres.append(str(t))
            lrepr = next_tips
            # print("NEXT TIPS", next_tips)
            if len(next_tips) == 0:
                break
        return lres

    def get_subtopologies(self):
        lrepr = self.to_list()
        lres = [lrepr]
        # print(lrepr)
        # return None
        while True:
            next_tips = []
            for t in lrepr:
                # print("T", t)
                # print("LRES", lres)
                if isinstance(t, list):
                    for tt in t:
                        # print("TT", tt)
                        # if tt not in lres:
                        next_tips.append(tt)
                    add_t = t
                    add_t.sort(key=lambda x: str(x))
                    lres.append(add_t)
                else:
                    lres.append(t)
            lrepr = next_tips
            # print("NEXT TIPS", next_tips)
            if len(next_tips) == 0:
                break
        return lres

    def get_subtopologies_as_qNodes(self):
        remaining = self
        lres = [remaining]
        # print(lrepr)
        # return None
        while True:
            next_tips = []
            for t in remaining:
                # print("T", t)
                # print("LRES", lres)
                if isinstance(t, qNode):
                    for tt in t:
                        # print("TT", tt)
                        # if tt not in lres:
                        next_tips.append(tt)
                    lres.append(t)
                # else:
                #    lres.append(t)
            remaining = next_tips
            # print("NEXT TIPS", next_tips)
            if len(next_tips) == 0:
                break
        return lres

    def tips_included(self):
        ti = [
            x for x in self.get_subtopologies_as_strings()
            if isinstance(x, str)
        ]
        ti.sort()
        return ti

    def is_synonymous_with(self, other):
        if all([
                s in other.get_subtopologies_as_strings()
                for s in self.get_subtopologies_as_strings()
        ]) and other.tips_included() == self.tips_included():
            return True

        else:
            return False

    def leaves(self):
        return self.node.get_leaves()


class shallowest_sistership_manager(object):

    def __init__(self, quartet_group_nodes):
        self.node = None
        self.qgi = None
        self.c = None

        self.quartet_group_nodes = quartet_group_nodes

        self.sister_pair = None
        self.left_out_group = None

    def _set(self, node, qgi, c):
        self.node = node
        self.qgi = qgi
        self.c = c

    def set(self, node, quartet_group_nodes, qgi, c):
        if self.node is None:
            self._set(node, qgi, c)
        else:
            if node is quartet_group_nodes[c]:
                self._set(node, qgi, c)

            else:
                pass

    def is_set(self):
        if self.node is not None:
            return True

        else:
            return False

    def get_left_out_groups(self):

        left_out_group = [
            x for x in self.quartet_group_nodes.keys()
            if str(x) not in [str(self.c), str(self.qgi)]
        ]

        if len(left_out_group) == 0:
            self.left_out_group = None
        else:
            self.left_out_group = left_out_group

        return self.left_out_group

    def commit(self, sister_pair=None):

        if self.qgi is None or self.c is None:
            raise ValueError(
                "One of `qgi` or `c` is None. Cannot commit to a sister pair.")

        if sister_pair is None:
            self.sister_pair = qNode([qNode(self.qgi.tips), self.c])
            left_out_group = self.get_left_out_groups()


def get_quartet_topology_first(quartet_group_nodes):
    shallowest_sistership = shallowest_sistership_manager(quartet_group_nodes)
    left_out_group = None

    # Generate the combinations of possible bifurcation topologies
    # base on the length of the quartet_group_nodes dictionary
    ncomps = len(quartet_group_nodes)
    # quartet_tips = list(quartet_group_nodes.keys())

    combs = [itertools.combinations(quartet_group_nodes.keys(), r=2)]
    '''
    combs = [
        ''.join(x) for x in list(
            itertools.combinations([str(x) for x in quartet_tips], r=2))
    ]
    '''
    # qr = {x: 0 for x in combs}
    # print(qr)

    sister_pair = None
    # print(qr)
    for qgi, node in quartet_group_nodes.items():
        if node is not None:
            compare_to = list(quartet_group_nodes.keys())
            # print(compare_to)
            compare_to.remove(qgi)
            compare_to.sort(key=lambda x: len(x.tips))
            print(f"> ANCHORED ON NODE {qgi}")
            print(f"NODE {qgi} LEAVES:", node.get_leaves())
            print("> Let us compare it to these nodes:", compare_to)
            sisters = node.get_sisters()
            sn = 0
            print(f"NODE {qgi} has {len(sisters)} sister nodes")
            print(f"NODE {qgi} is the root: {node.is_root()}")

            if node.is_root() and len(compare_to) == 1:
                # print(c, qgi.tips, compare_to[0], qgi)
                # sys.exit()
                shallowest_sistership.set(node, quartet_group_nodes, qgi,
                                          compare_to[0])
                sister_pair = qNode([qNode(qgi.tips), compare_to[0]])
                # print(c)
                left_out_group = [
                    x for x in quartet_group_nodes.keys()
                    if str(x) not in [str(compare_to[0]),
                                      str(qgi)]
                ]

                if len(left_out_group) == 0:
                    left_out_group = None

                print("LEFTOUT:", left_out_group)

                break

            # print(node.get_leaf_names())
            # if len(sis) == 1:
            for sis in sisters:
                print(f"> Looking at sister node #{sn}")
                print(f"Sister #{sn} num Leaves:", len(sis.get_leaves()))
                print(f"Sister #{sn} Leaves:", sis.get_leaves())
                # print(quartet_group_nodes.keys())
                for c in compare_to:
                    if quartet_group_nodes[c] is None:
                        continue

                    # If the two clades, together, are monophyletic
                    # print(sis.get_leaves())
                    if sis is quartet_group_nodes[c]:
                        print("first branch")
                        # print(f"{qgi} is sister to {c}")

                        # Generate the qr index so that different combinations aren't
                        # duplicated
                        # qr_idx = [int(x) for x in list(f"{qgi}{c}")]
                        # qr_idx.sort()
                        # qr_idx = [str(x) for x in qr_idx]
                        # qr_idx = ''.join(qr_idx)

                        # This is the shallowest sister relationship
                        assert sis.get_ancestors()[0] == node.get_ancestors(
                        )[0], "First ancestor of shallowest sisters is not the same"
                        print(
                            f"BRANCH 1: Setting shallowest_sistership between nodes {qgi} and {c}"
                        )
                        shallowest_sistership.set(sis.get_ancestors()[0],
                                                  quartet_group_nodes, qgi, c)
                        break
                print("-------")
                print(f"Done with Sister #{sn}")
                sn += 1
            # else:
            # for s in curr_node_sis:
            # print(s.get_leaves())
            # print("---")
            # print("----------")
            # pass
            # Break if shallowest sister relationship has
            # been found
            print(
                f"FIRST PASS ({shallowest_sistership.qgi}, {shallowest_sistership.c}): {shallowest_sistership.is_set()} {sn} {len(sisters)}"
            )
            if shallowest_sistership.is_set() and sn == len(sisters):

                # shallowest_sistership.commit()

                # print("LEFTOUT:", shallowest_sistership.left_out_group)

                sister_pair = qNode([
                    qNode(shallowest_sistership.qgi.tips),
                    shallowest_sistership.c
                ])
                # print(c)
                left_out_group = [
                    x for x in quartet_group_nodes.keys() if str(x) not in [
                        str(shallowest_sistership.c),
                        str(shallowest_sistership.qgi)
                    ]
                ]
                print("LEFTOUT:", left_out_group)
                if len(left_out_group) == 0:
                    left_out_group = None

                print(
                    f">>>>>>>LOCKING IN shallowest_sistership between nodes {shallowest_sistership.qgi} and {shallowest_sistership.c}"
                )
                break

        if sister_pair is not None:
            print("breaking because sister_pair is not None")
            break
            # print(qr_idx)
            # left_out_group = ''.join([
            #    str(x) for x in compare_to if str(x) not in qr_idx
            #])
        # Get the relationship of sister-pair to next sister
    print("SISTER PAIR: ", sister_pair)
    print(f"Group left out is : {left_out_group}")
    if shallowest_sistership.is_set() and sister_pair is not None:
        # print(left_out_group)
        if left_out_group is not None:
            # print("It's not none!")
            next_qgn = {}
            for g in left_out_group:
                next_qgn[g] = quartet_group_nodes[g]
            next_qgn[sister_pair] = shallowest_sistership.node

        else:
            next_qgn = {0: shallowest_sistership.node}
            print(next_qgn)
        # print(next_qgn)
        return (sister_pair, next_qgn)
    else:
        return None


def get_quartet_topology_second(quartet_group_nodes):
    shallowest_sistership = shallowest_sistership_manager(quartet_group_nodes)
    left_out_group = None

    # Generate the combinations of possible bifurcation topologies
    # base on the length of the quartet_group_nodes dictionary
    ncomps = len(quartet_group_nodes)
    # quartet_tips = list(quartet_group_nodes.keys())

    combs = [itertools.combinations(quartet_group_nodes.keys(), r=2)]
    '''
    combs = [
        ''.join(x) for x in list(
            itertools.combinations([str(x) for x in quartet_tips], r=2))
    ]
    '''
    # qr = {x: 0 for x in combs}
    # print(qr)

    sister_pair = None
    for qgi, node in quartet_group_nodes.items():
        if node is not None:
            compare_to = list(quartet_group_nodes.keys())
            # print(compare_to)
            compare_to.remove(qgi)
            compare_to.sort(key=lambda x: len(x.tips))
            print(f"> ANCHORED ON NODE {qgi}")
            print(f"NODE {qgi} LEAVES:", node.get_leaves())
            print("> Let us compare it to these nodes:", compare_to)
            sisters = node.get_sisters()
            sn = 0
            print(f"NODE {qgi} has {len(sisters)} sister nodes")
            print(f"NODE {qgi} is the root: {node.is_root()}")

            if node.is_root() and len(compare_to) == 1:
                # print(c, qgi.tips, compare_to[0], qgi)
                # sys.exit()
                shallowest_sistership.set(node, quartet_group_nodes, qgi,
                                          compare_to[0])
                sister_pair = qNode([qNode(qgi.tips), compare_to[0]])
                # print(c)
                left_out_group = [
                    x for x in quartet_group_nodes.keys()
                    if str(x) not in [str(compare_to[0]),
                                      str(qgi)]
                ]

                if len(left_out_group) == 0:
                    left_out_group = None

                print("LEFTOUT:", left_out_group)

                break

            # print(node.get_leaf_names())
            # if len(sis) == 1:
            for sis in sisters:
                print(f"> Looking at sister node #{sn}")
                print(f"Sister #{sn} num Leaves:", len(sis.get_leaves()))
                print(f"Sister #{sn} Leaves:", sis.get_leaves())
                # print(quartet_group_nodes.keys())
                for c in compare_to:
                    if quartet_group_nodes[c] is None:
                        continue

                    if quartet_group_nodes[c] in sis.get_descendants():
                        print("second branch")
                        # break
                        other_compares = list(compare_to)
                        other_compares.remove(c)
                        print("OTHER_COMPARES", other_compares)

                        other_compare_nodes = {
                            k: v
                            for k, v in quartet_group_nodes.items()
                            if str(k) in [str(x) for x in other_compares]
                        }
                        # print(other_compare_nodes)
                        # print(curr_node.get_leaves())
                        overlap_test = [
                            v not in sis.get_descendants()
                            for k, v in other_compare_nodes.items()
                        ]
                        # print(quartet_group_nodes[c].get_leaves())
                        print("OVERLAP TEST", overlap_test)
                        if all(overlap_test) and quartet_group_nodes[
                                qgi] in sis.get_ancestors()[0].get_descendants(
                        ):
                            print(
                                f"{qgi} is a descedant of a clade containing {c} (along with other tips), with size {sis.get_ancestors()[0].get_leaf_names()}"
                            )
                            print(
                                f"BRANCH 2: Setting shallowest_sistership between nodes {qgi} and {c}"
                            )
                            shallowest_sistership.set(sis.get_ancestors()[0],
                                                      quartet_group_nodes, qgi,
                                                      c)

                            break

                    # If the the two clades, together, are paraphyletic
                    else:
                        print("third branch")
                        print(
                            f"> Sister node #{sn} is not Node {c}. Interate back over ancestors of {qgi} and look in the descendants."
                        )
                        i = 0
                        # print(
                        #    len(list(quartet_group_nodes[c].iter_ancestors())))
                        depth = 1
                        for curr_node in quartet_group_nodes[c].iter_ancestors(
                        ):
                            print(
                                f">> Looking for Node {c} in descendants {depth} nodes back from {qgi}"
                            )
                            print(sister_pair)
                            # print(i)
                            # print(curr_node.get_leaves())
                            # i = i + 1
                            curr_node_sisters = curr_node.get_sisters()
                            # print("NODE", qgi, "LEAVES:",
                            # quartet_group_nodes[qgi].get_leaves())
                            # print("CURRENT NODE LEAVES:",
                            # curr_node.get_leaves())
                            # print("CURRENT NODE SISTERS:", curr_node_sisters)
                            # print("Curr_node_sis:", curr_node_sis)
                            # if len(curr_node_sisters) == 1:
                            # for s in curr_node_sisters:
                            # for s in curr_node_sis:
                            #    print(s.get_leaves())
                            #    print("---")
                            # print("-----")
                            if curr_node is quartet_group_nodes[c]:
                                print(
                                    f"{qgi} is sister to the {c} PLUS: {curr_node.get_leaf_names()}"
                                )

                            else:
                                other_compares = list(compare_to)
                                other_compares.remove(c)
                                print("OTHER_COMPARES", other_compares)

                                other_compare_nodes = {
                                    k: v
                                    for k, v in quartet_group_nodes.items() if
                                    str(k) in [str(x) for x in other_compares]
                                }
                                # print(other_compare_nodes)
                                # print(curr_node.get_leaves())
                                overlap_test = [
                                    v not in curr_node.get_descendants()
                                    for k, v in other_compare_nodes.items()
                                ]
                                # print(quartet_group_nodes[c].get_leaves())
                                print("OVERLAP TEST", overlap_test)

                                if quartet_group_nodes[
                                        c] in curr_node.get_descendants(
                                ) and all(
                                            overlap_test
                                ) and quartet_group_nodes[
                                            qgi] in curr_node.get_descendants(
                                ):
                                    print(
                                        f"{qgi} is a descedant of a clade containing {c} (along with other tips), with size {len(curr_node.get_leaf_names())}"
                                    )
                                    print(
                                        f"BRANCH 3: Setting shallowest_sistership between nodes {qgi} and {c}"
                                    )
                                    shallowest_sistership.set(
                                        curr_node, quartet_group_nodes, qgi, c)
                                    break

                                else:
                                    print(
                                        "In the sister:",
                                        quartet_group_nodes[c] in [
                                            x.get_descendants()
                                            for x in curr_node.get_sisters()
                                        ])
                            depth += 1
                print("-------")
                print(f"Done with Sister #{sn}")
                sn += 1
            # else:
            # for s in curr_node_sis:
            # print(s.get_leaves())
            # print("---")
            # print("----------")
            # pass
            # Break if shallowest sister relationship has
            # been found
            print(
                f"SECOND PASS ({shallowest_sistership.qgi}, {shallowest_sistership.c}): {shallowest_sistership.is_set()} {sn} {len(sisters)}"
            )
            if shallowest_sistership.is_set() and sn == len(sisters):

                sister_pair = qNode([
                    qNode(shallowest_sistership.qgi.tips),
                    shallowest_sistership.c
                ])
                # print(c)
                left_out_group = [
                    x for x in quartet_group_nodes.keys() if str(x) not in [
                        str(shallowest_sistership.c),
                        str(shallowest_sistership.qgi)
                    ]
                ]
                print("LEFTOUT:", left_out_group)
                if len(left_out_group) == 0:
                    left_out_group = None
                print(
                    f">>>>>>>LOCKING IN shallowest_sistership between nodes {shallowest_sistership.qgi} and {shallowest_sistership.c}"
                )
                break

        if sister_pair is not None:
            print("breaking because sister_pair is not None")
            break

    # Get the relationship of sister-pair to next sister
    print("SISTER PAIR: ", sister_pair)
    print(f"Group left out is : {left_out_group}")
    if shallowest_sistership.is_set() and sister_pair is not None:
        # print(left_out_group)
        if left_out_group is not None:
            # print("It's not none!")
            next_qgn = {}
            for g in left_out_group:
                next_qgn[g] = quartet_group_nodes[g]
            next_qgn[sister_pair] = shallowest_sistership.node

        else:
            next_qgn = {0: shallowest_sistership.node}
            print(next_qgn)
        # print(next_qgn)
        return (sister_pair, next_qgn)
    else:
        return None


def get_all_possible_topologies(all_tip_nodes):
    tips = list(all_tip_nodes)

    possible_topologies = list(tips)
    i = 1
    while i <= len(all_tip_nodes) - 1:
        for pt in itertools.combinations(tips, 2):
            if len(pt[0].tips_included()) + len(
                    pt[1].tips_included()) > len(all_tip_nodes):
                continue
            else:
                if pt not in possible_topologies:
                    # print(pt, pt[0].tips_included(), pt[1].tips_included())
                    if not any([
                            x in pt[1].tips_included()
                            for x in pt[0].tips_included()
                    ]):
                        # print(">>>", pt, pt[0].tips_included(), pt[
                        # 1].tips_included())
                        topo_to_add = qNode([pt[0], pt[1]])
                        if not any([
                                t.is_synonymous_with(topo_to_add)
                                for t in possible_topologies
                        ]):
                            possible_topologies.append(qNode([pt[0], pt[1]]))

                        tips.append(qNode([pt[0], pt[1]]))

        i += 1
        # print(tips)
    return possible_topologies


def sisterhood_heatmap(all_tip_nodes, topologies, heatmap_outpath):
    all_possible_sisters = [
        x for x in get_all_possible_topologies(all_tip_nodes)
        if len(x.tips_included()) < len(all_tip_nodes)
    ]
    print(all_possible_sisters)
    ldict = []
    index = []

    row_col_map = {}
    for pt in all_possible_sisters:
        row_col_map[pt] = str(pt)
        ldict.append({str(t): np.nan for t in all_possible_sisters})
        index.append(str(pt))

    df = pd.DataFrame(ldict, index=index)
    print([type(x) for x in df.index])
    for marker, topo in topologies.items():
        if topo is None:
            continue
        else:
            print(marker, topo)
            for sub in [
                    s for s in topo.get_subtopologies_as_qNodes()
                    if len(s.tips) > 1
            ]:
                row = None
                col = None
                for t, s in row_col_map.items():

                    if sub.tips[0].is_synonymous_with(t):
                        col = s
                    if sub.tips[1].is_synonymous_with(t):
                        row = s

                    if row is not None and col is not None:
                        break

                assert not any([x in sub.tips[1].tips_included() for x in sub.tips[0].tips_included(
                )]) or any([x in sub.tips[0].tips_included() for x in sub.tips[1].tips_included()]), "Topologies have duplicated tips."

                if np.isnan(df.loc[row, col]) and np.isnan(df.loc[col, row]):
                    df.loc[row, col] = 1
                else:
                    if np.isnan(df.loc[row, col]) and not np.isnan(df.loc[col, row]):
                        df.loc[col, row] += 1
                    elif np.isnan(df.loc[col, row]) and not np.isnan(df.loc[row, col]):
                        df.loc[row, col] += 1
                    else:
                        raise ValueError(
                            "When trying to make the heatmap: both df.loc[row,col] and df.loc[col,row] are not NaN. This should not have happened.")
                """
                if df.loc[row, col] == 0 and df.loc[col, row] != 0:
                    df.loc[col, row] += 1
                else:
                    df.loc[row, col] += 1
                """
                # print(pt, "|", sub.tips[0], row, type(row), df.index[
                #     row], "|", sub.tips[1], col, type(col), df.columns[col])
    df.to_csv(heatmap_outpath, sep="\t")


def fill_group_ancestral_nodes(marker: str, genetree: ete3.Tree, group_tips: dict) -> dict:
    logger = logging.getLogger(f"{marker}:fill_group_ancestral_nodes")
    quartet_group_nodes = {}
    for group_id, clades in group_tips.items():
        if clades is not None:

            # Sort clades by size, just in case there are more than 1 and we
            # have to pick the largest one
            clades.sort(key=lambda c: len(c), reverse=True)
            clade_sizes = [len(x) for x in clades]

            # Get the leaf names from the largest list of tips (i.e., clade)
            # We don't know if there is more than 1 list of tips yet, but it doesn't matter
            # because we're going to take the largest clade anyways.
            # Unless we skip because there are >1 clades of max size.
            leaves = [
                x for x in genetree.get_leaves()
                if x.isolate in clades[0]
            ]
            assert len(leaves) > 0, "0-length list of leaves"

            # A clade goes through here if its tips ARE monophyletic.
            # This means that it has exactly 1 list of tips in `clades` list
            # variable.
            if len(clades) == 1:

                if len(leaves) > 1:
                    quartet_group_nodes[qNode(int(group_id))] = genetree.get_common_ancestor(
                        leaves)
                else:
                    quartet_group_nodes[qNode(int(group_id))] = leaves[0]

            # A clade goes through here if its tips ARE NOT monophyletic.
            # This means that it has >1 list of tips in `clades` list variable.
            elif len(clades) > 1:

                logger.info(f"Group {group_id} is NOT monophyletic in gene tree `{marker}`. It is spread across {len(clades)} clades.")

                # This group is split over atleast 2 clades with the same, biggest size
                # We have not decided what to do with this yet, so skip it...
                if clade_sizes.count(max(clade_sizes)) > 1:
                    logger.warning(
                        f"More than one clade for non-monophyletic group {group_id} in gene tree `{marker}` is of maximum size. Skipping... "
                    )

                # This group is split over atleast 2 clades, but there is one single largest
                # Move forward with that one
                else:
                    logger.info(f"Taking largest clade {clades[0]} from list {clades}")

                    quartet_group_nodes[qNode(int(group_id))] = genetree.get_common_ancestor(
                        leaves)

            else:
                raise ValueError(f"Length of `clades` variable is 0: {clades}")
    return quartet_group_nodes


def quartet_repr(all_trees, json_path, outpath, outputbase):
    logger = logging.getLogger("quartet_repr")
    all_trees.midpoint_root()
    with open(json_path, 'r') as openfile:
        per_marker_quartet_group_monophyletic_clades = json.load(openfile)

    topologies = {}

    # There has to be an easier way to do this, but for now...
    # Get all possible group ids for making the heatmap later
    all_tip_nodes = [
        list(v.keys()) for _k, v in per_marker_quartet_group_monophyletic_clades.items() if len(v.keys()) == max([len(v.keys()) for _k, v in per_marker_quartet_group_monophyletic_clades.items()])][0]
    print(all_tip_nodes)
    all_tip_nodes = [qNode(int(x)) for x in all_tip_nodes]

    for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
        sister_pair = None

        gt = all_trees.trees[marker]
        quartet_group_nodes = fill_group_ancestral_nodes(marker=marker,
                                                         genetree=gt, group_tips=qg)

        if len(quartet_group_nodes.keys()) == 1:
            print("CRITICAL: There is only one group present.")

        while True:
            res = get_quartet_topology_first(quartet_group_nodes)
            if res is None:
                while True:
                    res2 = get_quartet_topology_second(quartet_group_nodes)
                    if res2 is None:
                        if sister_pair is not None:
                            print(
                                f"Unable to resolve topology past newick for: {sister_pair} |{marker}"
                            )
                        break
                    else:
                        sister_pair, quartet_group_nodes = res2
                        if len(quartet_group_nodes) == 1:
                            print(marker, "Newick:", sister_pair)
                            break
                break
            else:
                sister_pair, quartet_group_nodes = res
                if len(quartet_group_nodes) == 1:
                    print(marker, "Newick:", sister_pair)
                    break

        topologies[marker] = sister_pair

    topology_counts = {}
    marker2topology = {}

    topology_counts_outpath = os.path.join(outpath, f"{outputbase}_topology_counts.tsv")
    marker2topology_outpath = os.path.join(outpath, f"{outputbase}_marker2topology.tsv")
    heatmap_outpath = os.path.join(outpath, f"{outputbase}_topology_heatmap.tsv")

    for marker, topo in topologies.items():
        if topo is None:
            marker2topology[marker] = None
            continue
        inserted = False
        for q, c in topology_counts.items():
            if topo.is_synonymous_with(q):
                inserted = True
                topology_counts[q] += 1
                marker2topology[marker] = q
                break
        if not inserted:
            marker2topology[marker] = topo
            topology_counts[topo] = 1

    for topo, count in topology_counts.items():
        print(topo, count)

    logger.info(f"Writing topology counts to: {topology_counts_outpath}")
    with open(topology_counts_outpath, 'w') as topocounts:
        for topo, count in topology_counts.items():
            topocounts.write(f"{topo}\t{count}\n")

    logger.info(f"Writing marker topologies to: {marker2topology_outpath}")
    with open(marker2topology_outpath, 'w') as mark2topo:
        for marker, topo in marker2topology.items():
            mark2topo.write(f"{marker}\t{topo}\n")

    logger.info(f"Writing sisterhood heatmap to: {marker2topology_outpath}")
    sisterhood_heatmap(all_tip_nodes, topologies, heatmap_outpath)
    #    print("--------")

    # for k, v in qr.items():
    #    print(k, v)


if __name__ == "__main__":

    import logging

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("task")
    parser.add_argument("-g",
                        "--genetrees",
                        action="store",
                        required=True,
                        help="Path to directory containing gene trees.")
    parser.add_argument(
        "-r",
        "--hitreport",
        action="store",
        required=False,
        help="ONLY REQUIRED FOR `filter-phyly`. Path to hit_report_all.csv generated by domtbl2unaln."
    )
    parser.add_argument("-o",
                        "--outpath",
                        action="store",
                        required=False,
                        default=".",
                        help="Path to where output files will be written.")
    parser.add_argument("-b",
                        "--outputbase",
                        action="store",
                        required=False,
                        default="gtfilter",
                        help="SO FAR ONLY USED BY `quartet-repr`. Base name for writing output files.")
    parser.add_argument(
        "-m",
        "--phylymat",
        action="store",
        required=False,
        help="ONLY REQUIRED FOR `matcompare`. Another monophyly matrix to compare to this one."
    )
    parser.add_argument(
        "-i",
        "--isolates",
        action="store",
        required=False,
        help="ONLY REQUIRED FOR `taxocc`. Two-column, tab-separated list that has full tree tip label in column 1 and LTP in column 2."
    )
    parser.add_argument(
        "--remove-remaining-polyphyly",
        action="store_true",
        required=False,
        help="ONLY REQUIRED FOR `filter-phyly`. This flag will remove all unresolved polyphyletic taxa in each marker tree after completion of the filtering pipeline."
    )
    parser.add_argument(
        "-q",
        "--quartets",
        action="store",
        required=False,
        help="ONLY REQUIRED FOR `quartet-repr`, A four column tab-separated text files that has the tip labels for each quartet in a column."
    )
    parser.add_argument(
        "-j",
        "--json",
        action="store",
        required=False,
        default="monophyletic_groups.json",
        help="ONLY REQUIRED FOR `quartet-monophyly` or `quartet-group-monophyly`: Path to the input (`quartet-repr`) or output (`quartet-monophyly`) json file of monophyletic groupings of different quartet groups."
    )
    parser.add_argument(
        "--suffix",
        action="store",
        required=False,
        default=".aa.tre.renamed",
        help="Suffix for gene tree files. Default: `.aa.tre.renamed`")
    parser.add_argument(
        "--rooted",
        action="store_true",
        required=False,
        help="Pass this flag if the trees are rooted. ETE3 needs to know.")
    args = parser.parse_args()

    task_cases = [
        "phylymat", "matcompare", "filter-phyly", "filter-score", "taxocc",
        "quartet-repr", "quartet-group-monophyly",
        "print-nice-quartet-monophyly"
    ]
    if args.task not in task_cases:
        print(f"Bad task selection: {args.task}. Pick from {task_cases}.")
        sys.exit(1)
    files = [
        os.path.join(args.genetrees, t) for t in os.listdir(args.genetrees)
        if t.endswith(args.suffix)
    ]
    markers = [os.path.basename(x.replace(args.suffix, "")) for x in files]

    try:
        all_trees = gtlib.MultiMarkerGeneTrees(files,
                                               suffix=args.suffix,
                                               rooted=args.rooted)
    except IOError:
        print(
            f"There are NO tree located at `{args.genetrees}` that end with the suffix `{args.suffix}`"
        )

    if args.task == "phylymat":
        write_monophyly_matrix(all_trees, args.outpath)

    elif args.task == "matcompare":
        compare_monophyly_matrices(all_trees, args.outpath, args.phylymat)

    elif args.task == "filter-phyly":
        phyly_score_filter(all_trees, args.outpath, args.hitreport,
                           args.remove_remaining_polyphyly)

    elif args.task == "filter-score":
        hard_score_filter(all_trees, args.outpath, args.hitreport)

    elif args.task == "taxocc":
        taxon_occupancy(all_trees, args.isolates, args.outpath)

    elif args.task == "quartet-repr":
        quartet_repr(all_trees, args.json, args.outpath, args.outputbase)

    elif args.task == "print-nice-quartet-monophyly":
        print_nice_quartet_monophyly(args.json)

    elif args.task == "quartet-group-monophyly":
        write_per_marker_quartet_group_monophyletic_clades(
            all_trees, args.quartets, args.json)

    else:
        # Shouldn't happen
        pass
