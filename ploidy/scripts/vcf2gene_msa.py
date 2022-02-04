import vcf
import argparse
import scgid.sequence
import sys

parser = argparse.ArgumentParser()
parser.add_argument("vcf", action = "store", help = "Path to VCF file.")
parser.add_argument("--gene", action = "store", help = "Name of gene.")
parser.add_argument("-g", "--gff", action = "store", help = "Path to GFF3 gene map.")
parser.add_argument("-f", "--fna", action = "store", help = "Path to assembly fasta.")
parser.add_argument("--flanking", action = "store", help = "Length upstream and downstream of gene to include.", required = False, default = 0)
parser.add_argument("--samples", action = "store", help = "Samples to print sequences for.")
args = parser.parse_args()

flanking = int(args.flanking)

gene_coords = None
for line in open(args.gff, 'r'):
	if line.startswith("#"):
		continue
	else:
		spl = [x.strip() for x in line.split('\t')]
		inner_spl = spl[8].split(';')
		ident = inner_spl[1].split("=")[1]
		feature_type = spl[2]
		
		if ident == args.gene:
			if feature_type == "gene":
				gene_coords = (spl[0], int(spl[3]), int(spl[4]))

print(gene_coords)
if gene_coords is None:
	print("Gene not found.")
	sys.exit(1)

start = gene_coords[1] - flanking
end = gene_coords[2] + flanking
reference_sequence = None
print(start, end)
fasta = scgid.sequence.DNASequenceCollection().from_fasta(args.fna)
for s in fasta.seqs():
	if s.header.split(" ")[0] == gene_coords[0]:
		reference_sequence = s.string[start-1:end]

relative_start = start - start
relative_end = end - start - 1

samples = [x.strip() for x in args.samples.split(',')]
outseqs = {x:list(reference_sequence) for x in samples}
records_in_range = []
in_correct_chrom = False
for record in vcf.Reader(open(args.vcf, 'r')):
	if record.CHROM == gene_coords[0]:
		in_correct_chrom = True
		if record.POS >= start and record.POS <= end:

			assert len(record.ALT) == 1, "More than 1 alternative allele"
			assert len(record.ALT[0]) == 1, "Alternate allele length > 1"

			for s in record.samples:
				if s.sample in samples:
					if s.gt_type in [1,2]:
						if len(record.REF) > 1:
							print("Multiple character REF not implemented")
							sys.exit(1)
						else:
							idx = record.POS - start
							#print(outseqs[s.sample][idx-3:idx+3], record.REF, record.ALT)
							assert outseqs[s.sample][idx] == record.REF, "Sequence nucleotide not the same as REF"
							outseqs[s.sample][idx] = record.ALT[0]

	if in_correct_chrom:
		if record.CHROM != gene_coords[0]:
			break

for k,v in outseqs.items():
	seq = ''.join([str(x) for x in v])
	print(f">{k}\n{seq}")

