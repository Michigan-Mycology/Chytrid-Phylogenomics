import sys
import vcf
import argparse
from collections import namedtuple
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
parser.add_argument(
    "-f",
    "--features_tab",
    action="store",
    required=False,
    help=
    "Path to chromosomal features tabular, with geneID in column 1 and the information you want to map in column 2."
)
parser.add_argument("-o",
                    "--outfile",
                    action="store",
                    required=True,
                    help="Name of output file.")
parser.add_argument("vcf", action="store", help="Path to VCF.")
args = parser.parse_args()

features = {}

geneinfo = namedtuple("geneinfo", ["name", "desc"])
if args.features_tab is not None:
    with open(args.features_tab, 'r') as ft:
        for line in ft:
            if line.startswith("!"):
                continue

            spl = [x.strip() for x in line.split("\t")]
            inner_spl = [x.strip() for x in spl[1].split("|")]
            gene_name = inner_spl[0]
            gene_desc = inner_spl[1]

            if len(gene_name) == 0:
                gene_name = None

            features[spl[0]] = geneinfo(gene_name, gene_desc)

GQCUTOFF = 99
gt2code = {0: "R", 1: "H", 2: "A", None: "N"}

if args.genes is not None:
    goi = [[p.strip() for p in x.split('\t')]
           for x in open(args.genes).readlines()]
    goi_ids = {x[1]: x[0] for x in goi}

if args.samples is not None:
    samples_to_keep = args.samples.split(',')

reader = vcf.Reader(open(args.vcf, 'r'))

first_iteration = True
sample_order = []
rows = []
writer = vcf.Writer(open(f"{args.outfile.replace('.csv', '.vcf')}", 'w'),
                    reader)
for record in reader:

    gt_counts = {"R": 0, "H": 0, "A": 0, "N": 0}

    for i in record.INFO["ANN"]:

        this_row = []

        spl = [x.strip() for x in i.split("|")]
        annot = spl[1]
        impact = spl[2]
        geneid = spl[4]

        if args.genes is not None and args.samples is None:
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

                print(this_row)
                print(len(this_row))
                rows.append(this_row)

                first_iteration = False

        if args.genes is None and args.samples is not None:

            output_fname = "sample_snps.csv"

            if annot in [
                    "synonymous_variant", "upstream_gene_variant",
                    "downstream_gene_variant", "intergenic_region",
                    "intron_variant"
            ]:
                continue

            this_row.append(geneid)

            if args.features_tab is not None:
                this_row.append(features[geneid].name)
                this_row.append(features[geneid].desc)

            this_row.append(record.CHROM)
            this_row.append(record.POS)
            this_row.append(record.REF)
            this_row.append(record.ALT)
            this_row.append(annot)
            this_row.append(impact)

            for srec in record.samples:

                gt_to_append = None

                if srec.gt_type is None or srec["GQ"] < GQCUTOFF:
                    gt_to_append = gt2code[None]

                else:
                    gt_to_append = gt2code[srec.gt_type]

                gt_counts[gt_to_append] += 1

                if srec.sample in samples_to_keep:
                    this_row.append(gt_to_append)

                    if first_iteration:
                        sample_order.append(srec.sample)

            insert_pos = len(this_row) - len(samples_to_keep)
            this_row.insert(insert_pos, gt_counts["N"])
            this_row.insert(insert_pos, gt_counts["A"])
            this_row.insert(insert_pos, gt_counts["H"])
            this_row.insert(insert_pos, gt_counts["R"])

            ### Only record sample cluster SNPs that are not all identical to REF
            if not all([x == "R" for x in this_row[13:]]):
                writer.write_record(record)
                rows.append(this_row)

            if first_iteration:
                SAMPLE_ORDER = sample_order

            first_iteration = False

df = pd.DataFrame(rows)

columns = [
    "GeneID",
    "Chrom",
    "Pos",
    "Ref",
    "Alt",
    "Annotation",
    "Impact",
    "nR",
    "nH",
    "nA",
    "nN",
] + SAMPLE_ORDER

if args.features_tab is not None:
    columns.insert(1, "GeneDesc")
    columns.insert(1, "GeneName")

df.columns = columns
df.to_csv(args.outfile, sep="\t")
