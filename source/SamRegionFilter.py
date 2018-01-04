"""
This script defines a class to take a SAM file and BED
file as input and filter the SAM file such that it only contains
reads that start and/or end within regions defined in the bed file
"""

#Global

#Repos
import tools.io.file_iterator as file_iterator
from tools.formats.Bed import Bed

#Local

class SamRegionFilter:
    """This class filters an input SAM file such that it will only contain reads that start
    and/or end within regions of the input Bed file"""
    def __init__(self, sampath, bedpath):
        """Loads the input sampath and bedpath and starts filtering"""
        self.sampath = sampath
        self.bedpath = bedpath
        self.bed = Bed(bedpath)
        self.samlines = []

        self.filter_sam()

    """
    Filtering reads
    """

    def filter_sam(self):
        """Iterates through the input sam file and saves any lines
        that pass filters - that is, they either partially or fully overlap
        with the input bed regions"""
        for line in file_iterator.iterate(open(self.sampath)):
            if line[0] == "@":
                self.samlines.append(line)
                continue
            linevals = line.strip('\n').split('\t')
            chrom = linevals[2]
            if "chr" not in chrom:
                chrom = "chr"+chrom
            start = int(linevals[3])
            readlen = len(linevals[9])
            end = start + readlen
            if self.overlaps_bed(chrom, start, end):
                self.samlines.append(line)
            
    def overlaps_bed(self, chrom, start, end):
        """Checks if the supplied chromosome, start, and end
        positions overlap with self.bed"""
        for region in self.bed:
            bed_chrom = region.chromosome
            if "chr" not in bed_chrom:
                bed_chrom = "chr"+bed_chrom
            bed_start = region.start
            bed_end = region.end
            if chrom != bed_chrom:
                continue
            if start > bed_end or end < bed_start:
                continue
            return True
        return False
    
    """
    Saving filtered SAM lines
    """
    
    def save(self, outpath):
        """Saves the contents of self.samlines, which should be filtered,
        to outpath"""
        open(outpath, 'w').writelines(self.samlines)

if __name__ == "__main__":
    print("SamRegionFilter.py")
