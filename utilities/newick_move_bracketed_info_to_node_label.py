#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 15:15:29 2021

@author: aimzez
"""

#%%
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--treefile", action="store", required=True, help = "Path to annotated newick tree file.")
parser.add_argument("-o", "--output", action="store", required = True, help = "Output path for reformatted tree.")
args = parser.parse_args()

pat = re.compile("([0-9]+[:])([0-9.]+)([\[].*[\]])(.*)")

tree = open(args.treefile).readline()
output = args.output

out_list = []
for i in tree.split(")"):
    s = re.search(pat, i)
    if s is not None:
        bootstrap = s.group(1).replace(":", "")
        brlen = s.group(2)
        bracket_info = s.group(3).replace("[", "").replace("]","")
        bracket_info = bracket_info.replace("eqp-ic:", "/eqp-ic=")
        bracket_info = bracket_info.replace("qp-ic:", "/qp-ic=")
        bracket_info = bracket_info.replace("lq-ic:", "/iq-c=")
        bracket_info = bracket_info.replace(";", "")
        trailing = s.group(4)
        out = f"{bootstrap}{bracket_info}:{brlen}{trailing}"
        #print (i)
        #print(out)
        #print("-------")
        
        out_list.append(out)
    else:
        #pass
        out_list.append(i)
with open(output, 'w') as f:
    f.write(")".join(out_list))
