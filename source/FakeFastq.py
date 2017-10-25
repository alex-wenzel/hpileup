"""
This script takes a bed file and reference genome as input and produces a fastq file with simulated
reads based on the input bed regions.  The reads will be of length readlen and will share
a fraction overlap of their basepairs with the previous and next reads.
"""

#Global
from Bio import SeqIO

#Repos
import tools.formats.Bed as Bed

#Local

class FakeFastq:
    """Implements simulations of a fastq file with reads based on
    input regions with configurable length and overlap wtih each other"""
    def __init__(self, bedfile, refpath, readlen=1000, overlap=0.75):
        """Using the regions in bedfile, creates a fastq file with 'reads' generated
        from refpath of length readlen and overlap fraction overlap"""
        print("FakeFastq: Loading input bed...")
        self.bedfile = bedfile
        print("FakeFastq: Loading reference genome...")
        self.reference = SeqIO.index(refpath, "fasta")
        self.readlen = readlen
        self.overlap = overlap
        self.fqrecs = []

class FakeFastqRec:
    """A python representation of a simple fastq record containing the input
    fake seq, a unique recid, and a fake uniform high quality qual string"""
    def __init__(self, seq, recid):
        """Stores the read seq and unique recid and generates the uniform
        fake qual string"""
        self.seq = seq
        self.recid = recid
        self.qual = "~"*len(self.seq)

    def __str__(self):
        """Returns the string representation of the FakeFastqRec - the four
        lines with minimum data necessary to represent a fastq line"""
        ret_s = ""
        ret_s += "@"+self.recid+"\n"
        ret_s += self.seq+"\n"
        ret_s += "+\n"
        ret_s += self.qual+"\n"
        return ret_s

if __name__ == "__main__":
    print("FakeFastq.py")
    ff = FakeFastq("configs/input.bed", "../NA12878_PMS2/ncbi37.fa")
