import vcf
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samples", action="store", required=False,
                    help="Samples to retain from full VCF.")
parser.add_argument("vcf", action="store", help="Path to VCF.")
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

writer = vcf.Writer(open(f"{args.vcf}.filt", 'w'), reader)

sample_genotype_codes = []

first_iteration = True
i = 0
for record in reader:
    if len(record.ALT) > 1 or any([len(x) > 1 for x in record.ALT]):
        continue

    this_call_GT = []
    this_call_GQ = []
    #sample_order = []

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

        # if first_iteration:
        #	sample_order.append(srec.sample)

    # if first_iteration:
    #	SAMPLE_ORDER = sample_order

    nsamples = len(this_call_GT)
    nnan = len([x for x in this_call_GT if x == "N"])

    # If this Call's sample genotypes are mostly NA, skip over it
    if float(nnan) / float(nsamples) > PNANCUTOFF:
        pass

    else:
        writer.write_record(record)
        sample_genotype_codes.append(this_call_GT)

    i += 1
    if i % 5000 == 0:
        print(f"Processed {i} records.")
        print(f"Kept {len(sample_genotype_codes)} SNPs so far. ({round((float(len(sample_genotype_codes))/float(i))*100,2)}%).")
        print(f"-----")

    # if i == 10000:
    #   break

    first_iteration = False
