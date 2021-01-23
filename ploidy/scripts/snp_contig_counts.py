import argparse
from scgid.sequence import DNASequenceCollection

def asm_to_contig_lengths_dict(assembly):
    d = {}
    for i in assembly.seqs():
        d[i.header.split(" ")[0]] = len(i.string)

    return d

parser = argparse.ArgumentParser()
parser.add_argument("snp_stats_tsv", action="store", help="The *_snp_stats.tsv output file from running the `vcf_to_af.py` script included in this repo.")
parser.add_argument("-a", "--l50_assembly", action="store", required=True, help="Path to the L50 assembly - to get contig lengths.")
parser.add_argument("-o", "--output", action="store", required=True, help="Where to write counts.")
args = parser.parse_args()

tsv = open(args.snp_stats_tsv)

#Skip header line
tsv.readline()

asm = DNASequenceCollection().from_fasta(args.l50_assembly)
l50_headers = [x.header.split(" ")[0] for x in asm.seqs()]
out_tsv_path = args.output

contig_lengths = asm_to_contig_lengths_dict(asm)
contig_counts = {}

print(f"Working on {args.l50_assembly}")

for line in tsv:
    spl = [x.strip() for x in line.split("\t")]
    contig = spl[0]
    if contig in contig_counts:
        contig_counts[contig] += 1
    else:
        contig_counts[contig] = 1

for h in l50_headers:
    if h not in contig_counts:
        contig_counts[h] = 0

with open(out_tsv_path, 'w') as o:
    for k,v in contig_counts.items():
        snp_density = float(v)/float(contig_lengths[k])
        o.write(f"{k}\t{v}\t{snp_density}\n")

