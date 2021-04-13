import sys
from scgid.sequence import AASequenceCollection

SEQ_GAP_THRESHOLD = 0.75

aln = AASequenceCollection().from_fasta(sys.argv[1])
for s in aln.seqs():
    if float(s.string.count("-"))/float(len(s.string)) <= SEQ_GAP_THRESHOLD:
        print(s.to_fasta())

    else:
        # Remove sequence
        pass
