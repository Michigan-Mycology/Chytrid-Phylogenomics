import argparse
import os
import sys
from scgid.sequence import AASequenceCollection

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--delete_taxa", action = "store", required = True, help="Newline separated list of taxon LTPs to delete from unaligned fastas.")
parser.add_argument("-u", "--unaln_dir", action = "store", required = True, help="Directory containing unaligned fasta files.")
parser.add_argument("-o", "--outdir", action = "store", required = False, default=".", help="Output directory.")
args = parser.parse_args()

delete_taxa = []
with open(args.delete_taxa, 'r') as dt:
    for line in dt:
        delete_taxa.append(line.strip())

for fasta in [x for x in os.listdir(args.unaln_dir) if x.endswith(".fasta")]:
    seqs = AASequenceCollection().from_fasta(os.path.join(args.unaln_dir, fasta))
    fasta_outpath = os.path.join(args.outdir, fasta.replace(".fasta", ".fasta"))
    if os.path.isfile(fasta_outpath):
        print(f"ERROR: Writing filtered fasta files to `{os.path.abspath(args.outdir)}` would overwrite an existing fasta file. Please specify an outpath that does not contain existing fasta files!")
        sys.exit(1)
    with open(fasta_outpath, 'w') as outfasta:
        for seq in seqs.seqs():
            header_spl = seq.header.split("|")
            ltp = header_spl[0]
            if ltp not in delete_taxa:
                outfasta.write(f"{seq.to_fasta()}\n")
        
