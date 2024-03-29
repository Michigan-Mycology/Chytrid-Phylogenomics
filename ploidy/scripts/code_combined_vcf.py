import vcf
import sys
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("vcf", action = "store", help = "Path to vcf file.")
parser.add_argument("--write_snp_ident", action = "store_true", help = "Pass this flag to output identifying information of filtered SNPs.")
parser.add_argument("--mostly_dryrun", action = "store_true", help = "Pass this flag to process 10,000 VCf records and then exit.")
#parser.add_argument("--imp_mean", action = "store_true", help = "Pass this flag to immpute Call means for NA values.")
#parser.add_argument("--assume_ref", action = "store_true", help = "Pass this flag to assume REF allele for NA values.")

args = parser.parse_args()

GQCUTOFF = 99
PNANCUTOFF = 0.05

gt2code = {
	0: "R",
	1: "H",
	2: "A",
	None: "N"
}

reader = vcf.Reader(open(args.vcf, 'r'))

if args.write_snp_ident:
    snp_idents = []

sample_genotype_codes = []

first_iteration = True
i = 0
for record in reader:
	if len(record.ALT) > 1 or any([len(x) > 1 for x in record.ALT]):
		continue

	this_call_GT = []
	this_call_GQ = []
	sample_order = []

	if not "MQRankSum" in record.INFO:
		continue
	else:
		mqrs = record.INFO["MQRankSum"]
		if mqrs != 0.0:
			continue

	for srec in record.samples:

		if srec.gt_type is None or srec["GQ"] < GQCUTOFF:
			this_call_GT.append(gt2code[None])

		else:
			this_call_GT.append(gt2code[srec.gt_type])

		if first_iteration:
			sample_order.append(srec.sample)

	if first_iteration:
		SAMPLE_ORDER = sample_order

	nsamples = len(this_call_GT)

	nnan = len([x for x in this_call_GT if x == "N"])

	# If this Call's sample genotypes are mostly NA, skip over it
	if float(nnan)/float(nsamples) > PNANCUTOFF:
		pass

	else:

		sample_genotype_codes.append(this_call_GT)

		if args.write_snp_ident:
			snp_idents.append( (len(sample_genotype_codes), record.CHROM, record.POS) )

	i += 1
	if i % 5000 == 0:
		print(f"Processed {i} records.")
		print(f"Kept {len(sample_genotype_codes)} SNPs so far. ({round((float(len(sample_genotype_codes))/float(i))*100,2)}%).")
		print(f"-----")

	if args.mostly_dryrun:
		if i == 10000:
			break

	first_iteration = False

df = pd.DataFrame(sample_genotype_codes).transpose()
df.index = SAMPLE_ORDER

df.to_csv(f"pd_parsed_vcf_gq{GQCUTOFF}_pnan{PNANCUTOFF}.csv", sep = ",")

if args.write_snp_ident:
    with open(f"pd_parsed_vcf_gq{GQCUTOFF}_pnan{PNANCUTOFF}.snpidents.csv", 'w') as snpidout:
        for t in snp_idents:
            snpidout.write(f"{t[0]},{t[1]},{t[2]}\n")
