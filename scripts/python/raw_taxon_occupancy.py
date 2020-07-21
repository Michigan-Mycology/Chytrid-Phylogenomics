import pandas as pd
import os

def count_in_unaln(row):
    unaln_dir = "/home/aimzez/DATA/phylogeny/unaln"
    for f in os.listdir(unaln_dir):
        with open(f, 'r') as fasta:
            headers = [l for l in fasta.readlines() if l.startswith(">")]
            for h in headers:
                if row.LTP in h:
                    row.occurrences_in_unaln += 1
                    break
    print (row)
    return row

isolates = pd.read_excel("/home/aimzez/work/pursuit/sheets/Pursuit_Isolates.xlsx")
isolates = isolates[["SPECIES.TREE.LABEL", "LTP"]]
isolates["occurrences_in_unaln"] = 0
print(isolates)

isolates = isolates.apply(count_in_unaln, axis=1)
isolates.to_csv("taxon_repr_pre_trimal.tsv", sep="\t", index=False)
