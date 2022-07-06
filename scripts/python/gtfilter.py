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

# os.environ["CHYTRID_PHYLO_PY"] =
# "/home/aimzez/work/Chytrid-Phylogenomics/scripts/python"
sys.path.append(
    os.path.join(os.environ.get('CHYTRID_PHYLO'), "scripts", "python"))
import gtlib


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
                                                       quartets_path):
    # outgroup_isolates = ["GCF_000016345.1", "GCF_003065365.1", "GCF_000019745.1", "GCF_003261295.1"]
    # all_trees.root(outgroup_isolates)
    all_trees.midpoint_root()

    #for mark, tree in all_trees.trees.items():
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

    mono = [0] * len(quartet_groups)
    all_mono = 0
    per_marker_quartet_group_monophyletic_clades = {}
    for marker, tree in all_trees.trees.items():
        per_marker_quartet_group_monophyletic_clades[marker] = {}
        copy_qg = copy.deepcopy(quartet_groups)
        # print(copy_qg)
        mono_log = [True] * len(quartet_groups)
        for i in range(0, len(quartet_groups)):

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

    with open(
            "per_marker_quartet_group_monophyletic_clades_b_c_sep_inclOutgroup.json",
            "w") as outfile:
        json.dump(per_marker_quartet_group_monophyletic_clades, outfile)

    # for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
    #    for i, clades in qg.items():
    #        print(marker, i, clades)

    print(mono)
    print(all_mono)

    return None


def print_nice_quartet_monophyly():
    with open("per_marker_quartet_group_monophyletic_clades_b_c_sep.json",
              'r') as openfile:
        per_marker_quartet_group_monophyletic_clades = json.load(openfile)
        for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
            for i, clades in qg.items():
                print(marker, i, clades)


class qNode(object):

    def __init__(self, tips):
        self.tips = tips

    def __repr__(self):
        return str(self.tips).replace("[", "(").replace("]", ")")


def get_quartet_topology(quartet_group_nodes):
    shallowest_sistership = None
    left_out_group = None

    # Generate the combinations of possible bifurcation topologies
    # base on the length of the quartet_group_nodes dictionary
    ncomps = len(quartet_group_nodes)
    #quartet_tips = list(quartet_group_nodes.keys())

    print(quartet_group_nodes)
    combs = [itertools.combinations(quartet_group_nodes.keys(), r=2)]
    '''
    combs = [
        ''.join(x) for x in list(
            itertools.combinations([str(x) for x in quartet_tips], r=2))
    ]
    '''
    #qr = {x: 0 for x in combs}
    #print(qr)

    sister_pair = None
    # print(qr)
    for qgi, node in quartet_group_nodes.items():
        if node is not None:
            compare_to = list(quartet_group_nodes.keys())
            print(compare_to)
            compare_to.remove(qgi)
            sis = node.get_sisters()
            # print(node.get_leaf_names())
            if len(sis) == 1:
                for c in compare_to:
                    if quartet_group_nodes[c] is None:
                        continue

                    # If the two clades, together, are monophyletic
                    if sis[0] is quartet_group_nodes[c]:
                        # print(f"{qgi} is sister to {c}")

                        # Generate the qr index so that different combinations aren't
                        # duplicated
                        #qr_idx = [int(x) for x in list(f"{qgi}{c}")]
                        #qr_idx.sort()
                        #qr_idx = [str(x) for x in qr_idx]
                        #qr_idx = ''.join(qr_idx)

                        # This is the shallowest sister relationship
                        assert sis[0].get_ancestors()[0] == node.get_ancestors(
                        )[0], "First ancestor of shallowest sisters is not the same"
                        shallowest_sistership = sis[0].get_ancestors()[0]
                        #print(qr_idx)
                        #left_out_group = ''.join([
                        #    str(x) for x in compare_to if str(x) not in qr_idx
                        #])
                        left_out_group = [
                            x for x in quartet_group_nodes.keys()
                            if str(x) not in [str(c), str(qgi)]
                        ]
                        print(left_out_group)
                        #assert len(
                        #    left_out_group
                        #) <= 1, "More than one object in left_out_group"
                        if len(left_out_group) == 0:
                            left_out_group = None
                        #else:
                        #    left_out_group = left_out_group[0]
                        sister_pair = qNode([qNode(qgi), c])
                        break

                    # If the the two clades, together, are paraphyletic
                    else:
                        #continue
                        # print(sis)
                        for curr_node in quartet_group_nodes[c].iter_ancestors(
                        ):
                            curr_node_sis = curr_node.get_sisters()
                            #print("Curr_node_sis:", curr_node_sis)
                            if len(curr_node_sis) == 1:
                                if curr_node_sis[0] is quartet_group_nodes[c]:
                                    print(
                                        f"{qgi} is sister to the {c} PLUS: {curr_node.get_leaf_names()}"
                                    )
                                    break

                                else:
                                    if node in curr_node.get_descendants():
                                        print(
                                            f"{qgi} is a descedant of a clade containing {c} (along with other tips), with size {len(curr_node.get_leaf_names())}"
                                        )
                                        break

                    # Break if shallowest sister relationship has been found
                if sister_pair is not None:
                    break

            elif len(sis) < 1:
                continue

            else:
                print(">>>>>>>")
                print(
                    f"More than one sister clade for group {qgi}! {len(sis)} of them..."
                )
                for s in sis:
                    print(s.get_leaf_names())
                print("<<<<<<<")
                #raise ValueError("More than one sister clade!")

    # Get the relationship of sister-pair to next sister
    print(f"Group left out is : {left_out_group}")
    if shallowest_sistership is not None and sister_pair is not None:
        print(left_out_group)
        if left_out_group is not None:
            print("It's not none!")
            next_qgn = {
                sister_pair: shallowest_sistership,
            }
            for g in left_out_group:
                next_qgn[g] = quartet_group_nodes[g]

        else:
            next_qgn = {0: shallowest_sistership}
        print(next_qgn)
        return (sister_pair, next_qgn)
    else:
        return None


def quartet_repr(all_trees, quartets_path):
    all_trees.midpoint_root()
    with open(
            "per_marker_quartet_group_monophyletic_clades_b_c_sep_inclOutgroup.json",
            'r') as openfile:
        per_marker_quartet_group_monophyletic_clades = json.load(openfile)
    for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
        sister_pair = None
        quartet_group_nodes = {0: None, 1: None, 2: None}
        gt = all_trees.trees[marker]
        for i, clades in qg.items():
            if clades is not None:
                if len(clades) == 1:
                    leaves = [
                        x for x in gt.get_leaves() if x.isolate in clades[0]
                    ]

                    assert all([
                        x.is_leaf() for x in leaves
                    ]), "You got nonleaves mixed in with your leaves"
                    assert len(leaves) > 0, "0-length list of leaves"
                    # print(gt)
                    # print(leaves)
                    # print(marker, i, [x.isolate for x in leaves])
                    quartet_group_nodes[int(i)] = gt.get_common_ancestor(
                        leaves)
                else:
                    print(marker, len(clades), clades)
        quartet_group_nodes = {
            qNode(k): v
            for k, v in quartet_group_nodes.items()
        }
        print(quartet_group_nodes)
        next_qgn = quartet_group_nodes
        print(marker)
        while True:
            res = get_quartet_topology(next_qgn)
            if res is None:
                if sister_pair is not None:
                    print(
                        f"Unable to resolve topology past newick: {sister_pair}"
                    )
                break
            else:
                sister_pair, next_qgn = res
                if len(next_qgn) == 1:
                    print(marker, "Newick:", sister_pair)
                    break

        #print(marker, ">>>")
        #next_qgn = get_quartet_topology(quartet_group_nodes)
        #print(marker, next_qgn)
        # if next_qgn is not None:
        #    next_qgn = get_quartet_topology(next_qgn)
        #    print(marker, next_qgn)

        print("--------")

    # for k, v in qr.items():
    #    print(k, v)


if __name__ == "__main__":

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
        help=
        "ONLY REQUIRED FOR `filter-phyly`. Path to hit_report_all.csv generated by domtbl2unaln."
    )
    parser.add_argument("-o",
                        "--outpath",
                        action="store",
                        required=True,
                        help="Path to where output files will be written.")
    parser.add_argument(
        "-m",
        "--phylymat",
        action="store",
        required=False,
        help=
        "ONLY REQUIRED FOR `matcompare`. Another monophyly matrix to compare to this one."
    )
    parser.add_argument(
        "-i",
        "--isolates",
        action="store",
        required=False,
        help=
        "ONLY REQUIRED FOR `taxocc`. Two-column, tab-separated list that has full tree tip label in column 1 and LTP in column 2."
    )
    parser.add_argument(
        "--remove-remaining-polyphyly",
        action="store_true",
        required=False,
        help=
        "ONLY REQUIRED FOR `filter-phyly`. This flag will remove all unresolved polyphyletic taxa in each marker tree after completion of the filtering pipeline."
    )
    parser.add_argument(
        "-q",
        "--quartets",
        action="store",
        required=False,
        help=
        "ONLY REQUIRED FOR `quartet-repr`, A four column tab-separated text files that has the tip labels for each quartet in a column."
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

    all_trees = gtlib.MultiMarkerGeneTrees(files,
                                           suffix=args.suffix,
                                           rooted=args.rooted)

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
        quartet_repr(all_trees, args.quartets)

    elif args.task == "print-nice-quartet-monophyly":
        print_nice_quartet_monophyly()

    elif args.task == "quartet-group-monophyly":
        write_per_marker_quartet_group_monophyletic_clades(
            all_trees, args.quartets)

    else:
        # Shouldn't happen
        pass
