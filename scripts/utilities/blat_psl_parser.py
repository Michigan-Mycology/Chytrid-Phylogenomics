import argparse
import logging
import sys
from collections import namedtuple
class SyntenicCoverage(object):
    def __init__(self, total_length):
        self.total_length = total_length
        self.covered_ranges = {}

    def add_range(self, target, start, end):
        start = start-1
        end = end-1

        if target not in self.covered_ranges:
            self.covered_ranges[target] = [0]*self.total_length

        for i in range(start,end):
            if self.covered_ranges[target][i] == 0:
                self.covered_ranges[target][i] = 1

    def calc_syntenic_percentage (self):
        pass

class PslRecord(object):
    def __init__(self, line):

        self.raw = line

        spl = [x.strip() for x in line.split('\t')]
        self.matches = int(spl[0])
        self.misMatches = int(spl[1])
        self.repMatches = int(spl[2])
        self.nCount = spl[3]
        self.qNumInsert = int(spl[4])
        self.qBaseInsert = int(spl[5])
        self.tNumInsert = int(spl[6])
        self.tBaseInsert = int(spl[7])
        self.strand = spl[8]
        self.qName = spl[9]
        self.qSize = int(spl[10])
        self.qStart = int(spl[11])
        self.qEnd = int(spl[12])
        self.tName = spl[13]
        self.tSize = int(spl[14])
        self.tStart = int(spl[15])
        self.tEnd = int(spl[16])
        self.blockCount = spl[17]
        self.blockSizes = spl[18]
        self.qStarts = spl[19]
        self.tStarts = spl[20]

        # Make sure that alignment lengths are equal after considering gap bases
        if not self.assert_aln_lengths_match():
            logging.critical(f"Malformed PSL file format: Alignments lengths, minus gaps,  are not equal. Offending Entry: {line}")
            sys.exit()
    def q_aln_len_minus_gaps(self):
        return abs(self.qStart - (self.qEnd - self.qBaseInsert))

    def t_aln_len_minus_gaps(self):
        return abs(self.tStart - (self.tEnd -self.tBaseInsert))

    def assert_aln_lengths_match (self) -> bool:
        if self.q_aln_len_minus_gaps() != self.t_aln_len_minus_gaps():
            return False

        else:
            return True

    def percent_identity(self):
        gaps = self.qBaseInsert + self.tBaseInsert
        pident = (self.matches - self.misMatches) / (self.q_aln_len_minus_gaps() + gaps)
        pident_t = (self.matches - self.misMatches) / (self.t_aln_len_minus_gaps() + gaps)

        assert pident == pident_t, "Internal calculation error. Percent identities do not match."

        return pident

    def query_coverage(self):
        return self.q_aln_len_minus_gaps() / self.qSize

    def target_coverage(self):
        return self.t_aln_len_minus_gaps() / self.tSize

parser = argparse.ArgumentParser()
parser.add_argument("--psl", action = "store", required = True, help = "PSL file from blat")
args = parser.parse_args()

syntenic_pairs = {}
with open(args.psl) as psl:
    # Skip first 5 header lines
    for i in range(0,5):
        psl.readline()

    for line in psl:

        record = PslRecord(line)

        if record.q_aln_len_minus_gaps() < 250 or record.qName == record.tName:
            continue

        if record.qName not in syntenic_pairs:
            syntenic_pairs[record.qName] = SyntenicCoverage(record.qSize)

        if record.percent_identity() >= 0.95:
            syntenic_pairs[record.qName].add_range(record.tName, record.qStart, record.qEnd)
rm_len = 0
for contig,syncov in syntenic_pairs.items():
    for target,coverage in syncov.covered_ranges.items():
        duplicated_percent = coverage.count(1)/syncov.total_length
        if duplicated_percent > 0.50:
            rm_len += syncov.total_length
            print( '\t'.join([contig, target, str(duplicated_percent), str(syncov.total_length)]) )

print(f"TOTAL_REMOVED_LENGTH: {rm_len}")
