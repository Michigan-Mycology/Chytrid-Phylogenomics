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
sys.path.append(os.path.join(os.environ.get(
    'CHYTRID_PHYLO'), "scripts", "python"))
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
    that_frame = pd.read_csv(phylymat,  header=0).set_index("Unnamed: 0")

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
        print("[ERROR] Make sure that the first column of your hit report follows the format: LocusTagPrefix|Protein_ID")
        sys.exit(1)

    metadata["isolate"] = new[0]
    metadata["gene"] = new[1]
    return metadata


def phyly_score_filter(all_trees, outpath, hitreport_path, remove_remaining_polyphyletic_taxa):

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
    summed_renamed["stl"] = summed_renamed.ltp.apply(
        ltp2stl, mapper=ltp_mapper)
    summed_renamed = summed_renamed[["stl", "ltp", "occupancy"]]

    summed_renamed.to_csv(os.path.join(
        outpath, "python_isolate_repr.tsv"), sep="\t", index=False)


def write_per_marker_quartet_group_monophyletic_clades(all_trees, quartets_path):
    quartet_groups = {0: [],
                      1: [],
                      2: [],
                      3: [],
                      }
    # Read in quartet groups
    with open(quartets_path, 'r') as qf:
        for line in qf:
            spl = [x.strip() for x in line.split("\t")]
            for i in range(0, 4):
                if len(spl[i]) > 0:
                    quartet_groups[i].append(spl[i])

    mono = [0, 0, 0, 0]
    all_mono = 0
    per_marker_quartet_group_monophyletic_clades = {}
    for marker, tree in all_trees.trees.items():
        per_marker_quartet_group_monophyletic_clades[marker] = {}
        copy_qg = copy.deepcopy(quartet_groups)
        # print(copy_qg)
        mono_log = [True, True, True, True]
        for i in range(0, 3):
            per_marker_quartet_group_monophyletic_clades[marker][i] = []
            res = tree.check_monophyly(
                copy_qg[i], "isolate", ignore_missing=True, unrooted=True)
            if res[0]:
                mono[i] += 1
                per_marker_quartet_group_monophyletic_clades[
                    marker][i].append(tuple(quartet_groups[i]))
            else:
                p = 1
                while len(copy_qg[i]) != 0:
                    for comb in itertools.combinations(
                            copy_qg[i], len(copy_qg[i]) - p):
                        #print(len(copy_qg[i]), len(comb), p)
                        subtract_res = tree.check_monophyly(
                            comb, "isolate", ignore_missing=True, unrooted=True)
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

    with open("per_marker_quartet_group_monophyletic_clades.json", "w") as outfile:
        json.dump(per_marker_quartet_group_monophyletic_clades, outfile)

    # for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
    #    for i, clades in qg.items():
    #        print(marker, i, clades)

    print(mono)
    print(all_mono)

    return None


def quartet_repr(all_trees, quartets_path):
    with open("per_marker_quartet_group_monophyletic_clades.json", 'r') as openfile:
        per_marker_quartet_group_monophyletic_clades = json.load(openfile)
    for marker, qg in per_marker_quartet_group_monophyletic_clades.items():
        for i, clades in qg.items():
            print(marker, i, clades)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("task")
    parser.add_argument("-g", "--genetrees", action="store",
                        required=True, help="Path to directory containing gene trees.")
    parser.add_argument("-r", "--hitreport", action="store", required=False,
                        help="ONLY REQUIRED FOR `filter-phyly`. Path to hit_report_all.csv generated by domtbl2unaln.")
    parser.add_argument("-o", "--outpath", action="store", required=True,
                        help="Path to where output files will be written.")
    parser.add_argument("-m", "--phylymat", action="store", required=False,
                        help="ONLY REQUIRED FOR `matcompare`. Another monophyly matrix to compare to this one.")
    parser.add_argument("-i", "--isolates", action="store", required=False,
                        help="ONLY REQUIRED FOR `taxocc`. Two-column, tab-separated list that has full tree tip label in column 1 and LTP in column 2.")
    parser.add_argument("--remove-remaining-polyphyly", action="store_true", required=False,
                        help="ONLY REQUIRED FOR `filter-phyly`. This flag will remove all unresolved polyphyletic taxa in each marker tree after completion of the filtering pipeline.")
    parser.add_argument("-q", "--quartets", action="store", required=False,
                        help="ONLY REQUIRED FOR `quartet-repr`, A four column tab-separated text files that has the tip labels for each quartet in a column.")
    parser.add_argument("--suffix", action="store", required=False, default=".aa.tre.renamed",
                        help="Suffix for gene tree files. Default: `.aa.tre.renamed`")
    args = parser.parse_args()

    task_cases = ["phylymat", "matcompare",
                  "filter-phyly", "filter-score", "taxocc",
                  "quartet-repr", "quartet-group-monophyly"]
    if args.task not in task_cases:
        print(f"Bad task selection: {args.task}. Pick from {task_cases}.")
        sys.exit(1)

    files = [os.path.join(args.genetrees, t) for t in os.listdir(
        args.genetrees) if t.endswith(args.suffix)]
    markers = [os.path.basename(x.replace(args.suffix, ""))
               for x in files]

    all_trees = gtlib.MultiMarkerGeneTrees(files, suffix=args.suffix)

    if args.task == "phylymat":
        write_monophyly_matrix(all_trees, args.outpath)

    elif args.task == "matcompare":
        compare_monophyly_matrices(all_trees, args.outpath, args.phylymat)

    elif args.task == "filter-phyly":
        phyly_score_filter(all_trees, args.outpath,
                           args.hitreport, args.remove_remaining_polyphyly)

    elif args.task == "filter-score":
        hard_score_filter(all_trees, args.outpath, args.hitreport)

    elif args.task == "taxocc":
        taxon_occupancy(all_trees, args.isolates, args.outpath)

    elif args.task == "quartet-repr":
        quartet_repr(all_trees, args.quartets)
    
    elif args.task == "quartet-group-monophyly":
        write_per_marker_quartet_group_monophyletic_clades(all_trees, args.quartets)

    else:
        # Shouldn't happen
        pass
