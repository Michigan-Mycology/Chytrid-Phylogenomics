from scgid.sequence import AASequenceCollection
import sys
import os

SUFFIX = ".toget"

toget_dir = os.path.abspath(sys.argv[1])
pep_combined = AASequenceCollection().from_fasta(sys.argv[2])

for f in [x for x in os.listdir(toget_dir) if x.endswith(f"{SUFFIX}")]:
    path = os.path.join(toget_dir, f)
    wanted = [x.strip() for x in open(path).readlines()]
    marker = f.replace(f"{SUFFIX}","")
    out = [x for x in pep_combined.seqs() if x.header.split(" ")[0] in wanted]
    with open(f"{marker}.fasta", "w") as out_fasta:
        for s in out:
            out_fasta.write(s.to_fasta())
            out_fasta.write("\n")
