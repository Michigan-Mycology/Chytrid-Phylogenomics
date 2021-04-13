#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 18:06:50 2020

@author: Kevin Amses (amsesk@umich.edu)
"""
#%%
# Requires PyVCF and plotly
import vcf
import pandas as pd
import numpy as np
import sys
import argparse
import os
import subprocess

#%%


def in_select_contigs(contig: str, select_contigs: list):
    if contig in select_contigs:
        return True
    else:
        return False


def calc_allele_freqs(sample):
    ad = sample["AD"]
    p = ad[0]
    q = ad[1]

    try:
        p_freq = p / (p + q)
    except ZeroDivisionError:
        p_freq = 0.0

    try:
        q_freq = q / (p + q)
    except ZeroDivisionError:
        q_freq = 0.0

    return p_freq, q_freq


def q1_q3_dpfilt(df: pd.DataFrame):
    summary = df.depth.describe()
    print(summary)
    q1, q3 = (summary.loc["25%"], summary.loc["75%"])
    print(f"[Depth Filter] Filtering SNPs by depth. SNPs outside the inclusive range of ({q1}, {q3}) will be excluded.")
    return df[(df.depth >= q1) & (df.depth <= q3)]


def stdev_mqrsfilt(df: pd.DataFrame):
    summary = df.mqrs.describe()
    upper = summary["mean"] + summary["std"]
    lower = summary["mean"] - summary["std"]
    print(f"[MQRankSum Filter] Filtering SNPs by MQRankSum. SNPs outside the inclusive range of ({lower}, {upper}) will be excluded.")
    return df[(df.mqrs >= lower) & (df.mqrs <= upper)]


def strict_mqrsfilt(df: pd.DataFrame):
    print(f"[MQRankSum Filter] Filtering SNPs by MQRankSum strictly.")
    return df[df.mqrs == 0.0]

# Location of Chytrid Phylo directory
CHTRYIDPHYLO = os.path.join(os.environ.get('CHYTRID_PHYLO'))


#%% Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--select_contigs", action="store", required=False, default=None,
                    metavar="str", help="A list of contig headers to keep. SNPs on the rest will not be output.")
parser.add_argument("--skip_na_mqrs", action="store_true", required=False,
                    help="Skips contigs that don't have a value in the MQRankSum field.")
parser.add_argument("--dpfilt", action="store_true",
                    help="Filter SNPs by depth so that only SNPs with depths >= Q1 and <= Q3 are included.")
parser.add_argument("--mqrsfilt", action="store_true",
                    help="Filter SNPs by MQRankSum. Values outside 1 standard deviation from the mean will be excluded.")
parser.add_argument("-s", "--strain", action="store", required=True,
                    help="Strain name to display at the top of the final plots.")
parser.add_argument("-p", "--plot", action="store_true",
                    help="Print an allele frequency histogram to PDF using R.")
parser.add_argument("--position", action="store_true",
                    help="Print SNP position to snp_stats.tsv")
parser.add_argument("--kmercounts", action="store", required=False,
                    help="Optionally provide the path to and prefix of khist and peaks files from kmercountexact for incorporation into plots.")
parser.add_argument("vcf", metavar="PATH", help="Path to VCF file.")
args = parser.parse_args()

if args.select_contigs is not None:
    select_contigs = [x.strip() for x in open(args.select_contigs).readlines()]
    select_contigs = [x.split(" ")[0] for x in select_contigs]

#%% Variables
snp_stats_out = f"{args.strain}.snp_stats.tsv"

plot_script_path = os.path.join(
    CHYTRIDPHYLO, "ploidy", "scripts", "plot_af_hist.R")
#plot_script_path = "/home/aimzez/work/Chytrid-Phylogenomics/ploidy/scripts/plot_af_hist.R"

#%% Generate SNP stats
if not os.path.isfile(snp_stats_out):
    vcf_reader = vcf.Reader(open(args.vcf))

    llist = []
    for record in vcf_reader:

        # Skip contigs that aren't in the select_contigs list.
        # e.g., for L50-only selection
        if args.select_contigs is not None:
            if not in_select_contigs(record.CHROM, select_contigs):
                continue

        # Skip SNPs that don't have a MQRankSum value
        if "MQRankSum" not in record.INFO and args.skip_na_mqrs:
            continue

        for sample in record.samples:
            if len(record.get_hets()) != 0:
                ad = sample["AD"]
                dp = sample["DP"]

                if len(sample["AD"]) > 3 or len(record.REF) > 1 or len(record.ALT[0]) > 1:
                    continue
                else:

                    # Skip SNPs that don't have a MQRankSum value
                    if "MQRankSum" not in record.INFO:
                        if args.skip_na_mqrs:
                            continue
                        else:
                            mqrs = np.nan
                    else:
                        mqrs = record.INFO["MQRankSum"]

                    p, q = calc_allele_freqs(sample)

                    newrow = [record.CHROM, p, q, dp, mqrs]

                    if (args.position):
                        newrow.insert(1, record.POS)
                    llist.append(newrow)
                    # print(newrow)
    if args.position:
        df = pd.DataFrame(
            llist, columns=["contig", "position", "p", "q", "depth", "mqrs"])
    else:
        df = pd.DataFrame(llist, columns=["contig", "p", "q", "depth", "mqrs"])

    if df.shape[0] == 0:
        print("[ERROR] All SNPs have been filtered out by --skip_na_mqrs.")
        sys.exit(1)

    # Filter SNPs to only include SNPs with depths that are >= Q1 and <= Q3
    if args.dpfilt:
        df = q1_q3_dpfilt(df)

    # Filter SNPs to only include SNPs that MQRankSum that are within one standard
    # deviation of the mean.
    if args.mqrsfilt:
        df = strict_mqrsfilt(df)

    df.to_csv(snp_stats_out, sep="\t", index=False)

else:
    print(f"[WARNING] SNP stat file already exists at: {snp_stats_out}. If you asked for plot, we'll make it from the older file now. Otherwise we're done.")

#%% Generate plots in R
if args.plot:
    plotcmd = [
        "Rscript",
        "--vanilla",
        plot_script_path,
        snp_stats_out,
        args.strain
    ]

    # Add path/prefix for kmercount outputs if provided
    if args.kmercounts is not None:
        plotcmd.append(args.kmercounts)

    p = subprocess.Popen(plotcmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    out, err = p.communicate()
    if (p.returncode != 0):
        print("\n[[ R ENCOUNTERED AN ERROR, SEE BELOW ]]\n")
        print(err.decode("utf-8"))


#%%
'''
fig = go.Figure()
fig.add_trace(go.Histogram(x=df.p, nbinsx=100))
fig.add_trace(go.Histogram(x=df.q, nbinsx=100))
fig.update_layout(barmode="overlay")
fig.update_traces(opacity=0.75)
fig.show(renderer="svg")

#%%
freqs,depths = vcf_reader_get_freqs_and_depths(vcf_reader)

#%%
df = pd.DataFrame(freqs, columns = ["p", "q"])
df["dp"] = depths

df.to_csv("/home/aimzez/DATA/pursuit/ploidy/ARG085_freq_dp.csv", sep=",", index=False)

#%%
fig = go.Figure()
fig.add_trace(go.Scatter(x=df.p, y=df.dp, mode="markers", marker=dict(size=[1], color=["black"])))
fig.show(renderer="svg")
'''
