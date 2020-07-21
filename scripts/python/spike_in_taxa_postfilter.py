import sys
import os
import pandas as pd
import argparse
from scgid.sequence import AASequenceCollection, AASequence

FASTA_SUFFIX = "fasta"

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sheet", required=True, help="Sheet that says whether or not to spike in a each taxon.")
parser.add_argument("-r", "--raw_unaln", required=True, help="Directory that contains the RAW unaln fastas.")
parser.add_argument("-f", "--filt_unaln", required=True, help="Directory that contains the filtered unaln fastas.")
args = parser.parse_args()

sheet = pd.read_excel(args.sheet)
sheet = sheet[sheet["Spike.In"] == True]
taxa_to_spike = sheet.LTP.to_list()

for filt_fname in [x for x in os.listdir(args.filt_unaln) if x.endswith(FASTA_SUFFIX)]:
    path = os.path.join(args.filt_unaln, filt_fname)
    this_marker = filt_fname.split(".")[0]
    with open(path) as rf:
        headers = [x for x in rf.readlines() if x.startswith(">")]
        LTP_present = set([x.split("|")[0].replace(">","") for x in headers])
        if any([x not in taxa_to_spike for x in LTP_present]): 
            this_marker_filtered_fasta = AASequenceCollection().from_fasta(os.path.join(args.filt_unaln, f"{this_marker}.fasta"))
            this_marker_raw_fasta = AASequenceCollection().from_fasta(os.path.join(args.raw_unaln, f"{this_marker}.fasta"))
            for ltp in taxa_to_spike:
                if ltp in LTP_present:
                    continue
                for seq in this_marker_raw_fasta.seqs():
                    if seq.header.split("|")[0] == ltp:
                        this_marker_filtered_fasta.index[seq.header] = AASequence(seq.header, seq.string)
                        print(this_marker, f"Added {seq.header}")


            #Write the "spiked-in" fasta
            with open(f"{filt_fname}.spiked", 'w') as f:
                for seq in this_marker_filtered_fasta.seqs():
                    f.write(f"{seq.to_fasta()}\n")
 
