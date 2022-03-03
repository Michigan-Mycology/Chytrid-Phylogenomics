import argparse
import re
from scgid.sequence import DNASequenceCollection

WINDOW_SIZE = 500

parser = argparse.ArgumentParser()
parser.add_argument("--l50fasta", action="store", help="Path to L50 fasta.")
parser.add_argument("--snp_stats", action="store",
                    help="Path to snp_stats.tsv")
parser.add_argument("--strain", action="store",
                    help="Prefix for output file.")
args = parser.parse_args()

contig_order = []
with open(args.l50fasta, 'r') as l50fasta:
    for line in l50fasta:
        if line.startswith(">"):
            line = re.sub("^>", "", line)
            line = line.strip()
            line = line.split(" ")[0]
            contig_order.append(line)

l50fasta = DNASequenceCollection().from_fasta(args.l50fasta)
contig_lengths_d = {}
for s in l50fasta.seqs():
    contig_lengths_d[s.header.split(" ")[0]] = len(s.string)

llist = []
with open(args.snp_stats) as snp_stats:
    snp_stats.readline()
    for line in snp_stats:
        spl = [x.strip() for x in line.split("\t")]
        contig = spl[0]
        position = spl[1]
        llist.append([contig, contig_lengths_d[contig], position])

unique_contigs = []
for c in [x[0] for x in llist]:
    if c not in unique_contigs:
        unique_contigs.append(c)

for idx, contig in enumerate(contig_order):
    if contig not in unique_contigs:
        unique_contigs.insert(idx, contig)

assert len(unique_contigs) == len(
    l50fasta.index), "Number of unique contigs != Number of contigs in FASTA"

running_position = 0
with open(f"{args.strain}_snps_in_windows.tsv", 'w') as out:
    for c in unique_contigs:
        sub_llist = [x for x in llist if x[0] == c]

        if len(sub_llist) > 0:
            this_contig_length = set([x[1] for x in sub_llist])
            assert len(this_contig_length) == 1, f"More than one unique length for contig: {c}"
            this_contig_length = int(list(this_contig_length)[0])

        else:
            this_contig_length = contig_lengths_d[c]

        these_snp_positions = [int(x[2]) for x in sub_llist]

        start_points = [x for x in range(
            0, this_contig_length) if x % WINDOW_SIZE == 0]

        for start in start_points:
            end = start + WINDOW_SIZE
            if (end > this_contig_length):
                end = this_contig_length
            snps_in_window = [
                x for x in these_snp_positions if x > start and x <= end]

            out.write('\t'.join([c, str(start), str(
                end), str(len(snps_in_window)), str(running_position)]))
            out.write("\n")

        running_position += this_contig_length
