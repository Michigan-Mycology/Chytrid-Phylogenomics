import scgid.sequence
import sys
import numpy as np

full = scgid.sequence.DNASequenceCollection().from_fasta(sys.argv[1])

seqlist = list(full.index.values())
seqlist.sort(key = lambda x: len(x.string), reverse=True)
seqlens = np.array([len(s.string) for s in seqlist])

total = seqlens.sum()
halfway = np.ceil(total/2)
cumsum = seqlens.cumsum()

l50_idx = 0

for i,cs in enumerate(cumsum):
    if cs > halfway:
        l50_idx = i
        break

for s in seqlist[0:l50_idx+1]:
    print (s.to_fasta())
