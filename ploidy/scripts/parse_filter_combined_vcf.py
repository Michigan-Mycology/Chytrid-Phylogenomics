import vcf
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("vcf", action = "store", help = "Path to vcf file.")
args = parser.parse_args()

reader = vcf.Reader(open(args.vcf, 'r'))


sample_genotype_codes = {}

first_iteration = True
i=0
for record in reader:
	if  len(record.ALT) > 1 or any([len(x) > 1 for x in record.ALT]):
		continue
	for srec in record.samples:

		if first_iteration:
			sample_genotype_codes[srec.sample] = ""

		#genotype = srec.data.GT
		#genoqual = srec.data.GQ

		# hom_ref = 0
		# het = 1
		# hom_alt = 2
		if srec.gt_type is None:
			sample_genotype_codes[srec.sample] = sample_genotype_codes[srec.sample] + "-"

		else:
			sample_genotype_codes[srec.sample] = sample_genotype_codes[srec.sample] + str(srec.gt_type)

	first_iteration = False
	i+=1

	if i == 1000:
		break


for key,v in sample_genotype_codes.items():
	print(f">{key}\n{v}")


