#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:13:55 2020

@author: aimzez
"""

import os
import argparse
from scgid.sequence import AASequenceCollection, AASequence

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignments", action = "store", required=True, help="Directory that contains the multiple sequence alignments to concatenate.")
parser.add_argument("-o", "--output", action = "store", required=False, default="concatenated_alignment.fasta", help = "Output filename and path.")
parser.add_argument("-s", "--suffix", action="store", required=False, default=".fasta", help="MSA filename suffix.")
args = parser.parse_args()

msa_list = [os.path.join(args.alignments,x) for x in os.listdir(args.alignments) if x.endswith(args.suffix)]

concat_aln_len = 0
concat_aln = AASequenceCollection()

partitions = []

for msa_path in msa_list:
    added_this_round = []
    this_msa = AASequenceCollection().from_fasta(msa_path)
    this_msa_seqs = this_msa.seqs()
    this_msa_len = len(next(iter(this_msa_seqs)).string)
    marker = os.path.basename(msa_path).split('.')[0]
    for s in this_msa_seqs:
        ltp = s.header.split("|")[0]
        
        assert ltp not in added_this_round, "Multiple sequences corresponding to `{ltp}` in `{msa_path}`."
        assert len(s.string) == this_msa_len, "Sequences in `{msa_path}` are not all the same length."
        
        if ltp in concat_aln.index:
            concat_aln.index[ltp].string += s.string
        
        else:
            new_sequence = f"{'-'*concat_aln_len}{s.string}"
            concat_aln.index[ltp] = AASequence(ltp, new_sequence)
        
        added_this_round.append(ltp)
    
    partitions.append( (marker, concat_aln_len+1, concat_aln_len+this_msa_len))

    concat_aln_len += this_msa_len
    
    for missing_ltp in [x for x in concat_aln.index.keys() if x not in added_this_round]:
        concat_aln.index[missing_ltp].string += '-'*this_msa_len
    

with open(args.output, 'w') as out:
    sorted_ltp = list([x.header for x in concat_aln.seqs()])
    sorted_ltp.sort()
    for ltp in sorted_ltp:
        out.write(f"{concat_aln.index[ltp].to_fasta()}\n")

with open(f"{args.output}.partitions.nex", 'w') as parts:
    parts.write("#nexus\n")
    parts.write("begin sets;\n")
    for marker,start,end in partitions:
        parts.write(f"\tcharset {marker} = {start}-{end};\n")
    parts.write(f"\tcharpartition mine = {None};\n")
    parts.write("end;\n")
        

