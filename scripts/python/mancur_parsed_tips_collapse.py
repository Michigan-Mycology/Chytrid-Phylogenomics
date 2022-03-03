#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 10:21:02 2020

@author: aimzez
"""

#%%
import pandas as pd
import numpy as np
import sys

tsv = pd.read_csv(sys.argv[1], sep="\t")
tsv = tsv[tsv.category == "SALVAGEABLE"]
tsv_g = tsv.groupby("marker")

#%%
tsv = tsv_g.agg({
    "marker": lambda x: x[0],
    "category": lambda x: x[0],
    "delete_tips": lambda x: ";".join(x)
        })
tsv

#%%
tsv.to_csv(sys.argv[1].replace('.tsv','_collapsed.tsv'), sep="\t", index=False)
