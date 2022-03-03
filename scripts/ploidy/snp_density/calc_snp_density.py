import argparse
import os
from scgid.sequence import DNASequenceCollection

parser = argparse.ArgumentParser()
parser.add_argument("--assembly", action="store", required=True, help = "Path to assembly.")
parser.add_argument("--parsed_vcf", action="store", required=True, help = "Path to partsed VCF.")
parser.add_argument("--strain", action="store", required=True, help = "Strain name.")
args = parser.parse_args()

asmlen = sum([len(x.string) for x in DNASequenceCollection().from_fasta(args.assembly).seqs()])
nsnps = len(open(args.parsed_vcf).readlines()) - 1 #Minus 1 for column headers
snp_density = float(nsnps)/float(asmlen)

print (','.join( [str(x) for x in [args.strain, os.path.basename(args.assembly), asmlen, nsnps, snp_density]] ))
