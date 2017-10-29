"""
This class takes a FASTA file as input and a bed file and
produces a subset fasta file with only the sequences from the regions
in the bedfile
"""

#Global
from Bio import SeqIO
import sys

#Repos
import tools.formats.Bed as Bed

#Local

class FastaSubset:
    """A class that builds a representation of a FASTA file containing only regions
    from an input bed file"""
    def __init__(self, fastapath, bedpath):
        print("FastaSubset: Loading input fasta file...")
        self.fastapath = fastapath
        self.bedpath = bedpath
        self.input_fasta = SeqIO.index(self.fastapath, "fasta")
        self.input_bed = Bed.Bed(self.bedpath)
        self.lines = []
        
        self.parse()
    
    """
    Sequence retrieval
    based on bed lines
    """
        
    def parse(self):
        """Iterates through self.input_bed and retrieves the relevant
        sequences from self.input_fasta"""
        recid = 1
        for bedline in self.input_bed:
            chrom = bedline.chromosome
            try:
                seq = self.input_fasta[chrom].seq
            except KeyError:  #likely naming convention, try a couple possibilities
                try:
                    if "chr" in chrom:  #check if input fasta convention is without chr
                        try:
                            seq = self.input_fasta[chrom.replace("chr", "")].seq
                        except KeyError:
                            print("ERROR: FastaSubset: Could not find "+chrom+" in "+self.fastapath)
                            print("\tCouldn't find "+chrom.replace("chr", "")+" either.")
                            sys.exit(1)
                    else:  #check if input fasta convention uses chr
                        try:
                            seq = self.input_fasta["chr"+chrom].seq
                        except KeyError:
                            print("ERROR: FastaSubset: Could not find "+chrom+" in "+self.fastapath)
                            print("\tCoudln't find chr"+chrom+" either.")

if __name__ == "__main__":
    print("FastaSubset.py")
