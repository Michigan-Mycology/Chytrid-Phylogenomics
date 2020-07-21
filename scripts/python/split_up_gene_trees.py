#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:50:22 2020

@author: aimzez
"""

#%%
import numpy as np
import os

#%%

names = ["Rabern", "Tim", "Gus"]

files = np.array(os.listdir("/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/pdf"))
spl = np.array_split(files, len(names))

i=0
for subarray in spl:
    with open(f"{names[i]}_markers.txt", 'w') as f:
        for filename in subarray:
            marker = filename.split(".")[0]
            f.write(f"{marker}\n")
    i+=1
    

#%%
all_files = []
for i in names:
    all_files = all_files + os.listdir(f"{i}")
    
#%%
print(len(all_files))
print(len(set(all_files)))
all_markers = ([x.split(".")[0] for x in all_files])
print(len(all_markers))
print(len(set(all_markers)))