from scgid.sequence import DNASequenceCollection
import sys
import os

wanted = sys.argv[1]
if os.path.isfile(wanted):
    wanted = open(sys.argv[1]).readlines()
    wanted = [x.strip() for x in wanted]
else:
    wanted = wanted.split(',')

allseqs = DNASequenceCollection().from_fasta(sys.argv[2], spades=False)

out = [x for x in allseqs.seqs() if x.header.split(" ")[0] in wanted]
for i in out:
    print (i.to_fasta())
