import sys
sys.path.append("/home/amsesk/scripts/")
from vcf_record import VCFRecord
from scgid.sequence import DNASequenceCollection

l50_fasta = sys.argv[1]
vcf = sys.argv[2]

headers = [s.header for s in DNASequenceCollection().from_fasta(l50_fasta).index.values()]

with open(vcf, 'r') as v:
    for line in v:
        spl = [x.strip() for x in line.split("\t")]
        if spl[0] in headers:
            print (line.strip())

