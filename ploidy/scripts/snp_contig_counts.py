import argparse

parser = argparse.ArgumentParser()
parser.add_argument("snp_stats_tsv", action="store", help="The *_snp_stats.tsv output file from running the `vcf_to_af.py` script included in this repo.")
parser.add_argument("-o", "--output", action="store", required=True, help="Where to write counts.")
args = parser.parse_args()

tsv = open(args.snp_stats_tsv)
out_tsv_path = args.output

contig_counts = {}

for line in tsv:
    spl = [x.strip() for x in line.split("\t")]
    contig = spl[0]
    if contig in contig_counts:
        contig_counts[contig] += 1
    else:
        contig_counts[contig] = 1

with open(out_tsv_path, 'w') as o:
    for k,v in contig_counts.items():
        o.write(f"{k}\t{v}\n")

