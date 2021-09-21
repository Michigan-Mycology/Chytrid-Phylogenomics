import sys
import vcf
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s",
                    "--samples",
                    action="store",
                    required=False,
                    help="Samples to retain from full VCF.")
parser.add_argument("-g",
                    "--genes",
                    action="store",
                    required=False,
                    help="Table of genes of interest.")
parser.add_argument("vcf", action="store", help="Path to VCF.")
args = parser.parse_args()

GQCUTOFF = 99
gt2code = {0: "R", 1: "H", 2: "A", None: "N"}

goi = [[p.strip() for p in x.split('\t')]
       for x in open(args.genes).readlines()]
goi_ids = {x[1]: x[0] for x in goi}

reader = vcf.Reader(open(args.vcf, 'r'))

first_iteration = True
sample_order = []
rows = []
for record in reader:

    this_row = []

    for i in record.INFO["ANN"]:
        spl = [x.strip() for x in i.split("|")]
        annot = spl[1]
        impact = spl[2]
        geneid = spl[4]
        if geneid in goi_ids.keys():
            this_row.append(goi_ids[geneid])
            this_row.append(record.CHROM)
            this_row.append(record.POS)
            this_row.append(record.REF)
            this_row.append(record.ALT)
            this_row.append(annot)
            this_row.append(impact)

            for srec in record.samples:

                if srec.gt_type is None or srec["GQ"] < GQCUTOFF:
                    this_row.append(gt2code[None])

                else:
                    this_row.append(gt2code[srec.gt_type])

                if first_iteration:
                    sample_order.append(srec.sample)

            if first_iteration:
                SAMPLE_ORDER = sample_order

            rows.append(this_row)

            first_iteration = False

df = pd.DataFrame(rows).transpose()
df.index = ["GeneID", "Chrom", "Pos", "Ref", "Alt", "Annotation", "Impact"
            ] + SAMPLE_ORDER

df.to_csv("goi_snp_dist.csv", sep=",")
