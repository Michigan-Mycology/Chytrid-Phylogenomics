import vcf
import sys

reader = vcf.Reader(open(sys.argv[1], 'r'))

N_MULTI_ALLELIC = 0
N_DI_ALLELIC = 0

NSNPS = 0
NDELS = 0
NINSE = 0

i = 0
for record in reader:
    i += 1
    if  len(record.ALT) > 1:
        N_MULTI_ALLELIC += 1

    else:
        N_DI_ALLELIC += 1

    alt_lens = [len(x) for x in record.ALT]
    ref_len = len(record.REF[0])
    assert isinstance(record.REF, str), f"Ref is not str"
    if any([x > 1 for x in alt_lens]) or ref_len > 1:
        for a in alt_lens:
            if a > ref_len:
                NINSE += 1
            elif a < ref_len:
                NDEL += 1
            else:
                NSNPS += 1

    else:
        for a in alt_lens:
            assert a == 1, "woops"
            NSNPS += 1

    if i % 5000 == 0:
        print(f"Processed {i} lines.")
        print(f"So far, {N_MULTI_ALLELIC} multi and {N_DI_ALLELIC} diallelic variants.")

print(f"""
>Diallelic variants: {N_MULTI_ALLELIC}
Diallelic variants: {N_DI_ALLELIC}
Insertion Alleles: {NINSE}
Deletion Alleles: {NDELS}
Single nucleotide variants: {NSNPS}
""")

