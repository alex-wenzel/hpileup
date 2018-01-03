"""
This script takes a bed file and reference genome as input and produces a fastq file with simulated
reads based on the input bed regions.  The reads will be of length readlen and will share
a fraction overlap of their basepairs with the previous and next reads.
"""

#Global
from Bio import SeqIO
import sys

#Repos
import tools.formats.Bed as Bed
from tools.io.Log import Log

#Local

class FakeFastq:
    """Implements simulations of a fastq file with reads based on
    input regions with configurable length and overlap wtih each other"""
    def __init__(self, bedfile, refpath, readlen=1000, overlap=0.75):
        """Using the regions in bedfile, creates a fastq file with 'reads' generated
        from refpath of length readlen and overlap fraction overlap"""
        self.l = Log()
        self.l.log("FakeFastq: Loading input bed...", tabs=1)
        self.bedfile = bedfile
        self.l.log("FakeFastq: Loading reference genome...", tabs=1)
        self.reference = SeqIO.index(refpath, "fasta")
        self.readlen = readlen
        self.overlap = overlap
        self.fqrecs = []

        self.l.log("FakeFastq: Generating reads...", tabs=1)
        self.build_recs()

    """
    Record Building
    """

    def build_recs(self):
        """Iterates through the bedfile and populates fqrecs with 
        FakeFastqRec objects corresponding to the bed regions"""
        for bedline in self.bedfile:
            try:
                chromseq = str(self.reference[bedline.chromosome].seq)
            except KeyError:
                self.l.error("FakeFastq: Chromosome '"+bedline.chromosome+"' not found in reference fasta", 
                                tabs=1, die=True, code=1)
            seq = chromseq[bedline.start-1:bedline.end]
            self.l.log("FakeFastq: Building reads for "+str(bedline)[:-1]+"...", tabs=1)
            self.build_reads(seq)
        self.l.log("FakeFastq: Done with all regions", tabs=1)
    
    def build_reads(self, seq):
        """Takes the sequence corresponding to the current bedfile in build_recs and populates
        self.fqrecs with FakeFastqRec objects"""
        indx = 0  #current position in seq (reads start from here)
        rid = 1  #current unique record id for building FakeFastqRec objects
        window_incr = self.readlen - int(self.readlen*self.overlap)  #the amount to increase indx on each iteration
        while True:
            subseq = seq[indx:indx+self.readlen]  #a section of the sequence of readlen length
            if len(subseq) == 0:
                break  #we've read the entire sequence and generated all necessary reads
            #save a FakeFastqRec object with the subsequence and the next unique identifier
            self.fqrecs.append(FakeFastqRec(subseq, "r"+str(rid)))
            if len(subseq) < self.readlen:
                break  #this is the last read
            rid += 1  #increment the unique record identifier
            indx += window_incr  #move the start of the next read

    """
    Saving reads
    """

    def save(self, outpath):
        open(outpath, 'w').write(str(self))
        self.l.log("FakeFastq: All regions saved to "+outpath, tabs=1)

    """
    Operators
    """

    def __str__(self):
        return ''.join([str(read) for read in self.fqrecs])

class FakeFastqRec:
    """A python representation of a simple fastq record containing the input
    fake seq, a unique recid, and a fake uniform high quality qual string"""
    def __init__(self, seq, recid):
        """Stores the read seq and unique recid and generates the uniform
        fake qual string"""
        self.seq = seq
        self.recid = recid
        self.qual = "~"*len(self.seq)

    """
    Operators
    """

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
