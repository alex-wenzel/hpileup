"""
This class implements the methods to determine the simplest
and most complete expression of homolog complexity given the input
regions and the aligned artificial reads

See prototype in homolog_mapping3.ipynb
"""

#Global
import pysam

#Repos
import tools.formats.Bed as Bed
from tools.io.Log import Log

#Local

class QnameMaps:
    """Takes a SAM file as input and parses the reads into
    a dictionary of the format
    {qname: {"input": (chrom, start, end), "homs": [(chrom, start, end)]}}"""
    def __init__(self, sampath, bedpath):
        """Loads the bed and sam and begins parsing according to the format above"""
        self.sampath = sampath
        self.bedpath = bedpath
        self.sam = pysam.AlignmentFile(self.sampath, mode='r')
        self.input_regions = Bed.Bed(self.bedpath)
        self.maps = {}

        self.parse()

    """
    SAM Parsing
    """

    def parse(self):
        """Iterates through self.sam and populates self.maps"""
        unmapped = []  #[(qname, chrom, start, end)]
        for rec in self.sam.fetch():
            chrom, start, end = rec.reference_name, rec.reference_start, rec.reference_end
            qname = rec.query_name
            try:
                self.maps[qname]['homs'].append((chrom, start, end))
            except KeyError:
                if self.input_regions.overlap(chrom, start, end):
                    self.maps[qname] = {"input": (chrom, start, end), "homs": []}
                else:
                    unmapped.append((qname, chrom, start, end))
        for rec in unmapped:
            qname, chrom, start, end = rec
            try:
                self.maps[qname]['homs'].append((chrom, start, end))
                unmapped.remove(rec)
            except KeyError:
                pass

    """
    Operators
    """

    def __getitem__(self, key):
        return self.maps[key]

    def keys(self):
        return self.maps.keys()

class HomMap:
    """This class represents a single input region and a set of homologous regions"""
    def __init__(self, chrom, start, end, homs):
        """Store initialization parameters"""
        self.chrom = chrom
        self.start = start
        self.end = end
        self.homs = homs

    """
    Check if this HomMap can merge with another HomMap
    """

    def can_merge(self, hm):
        """Checks the three rules and returns a new HomMap object if this HomMap
        object can merge with another input"""
        new_hm = self.region_overlap(hm)
        if new_hm == None:
            return None
        if len(self.homs) != len(hm.homs):
            return None
        merged_homs = self.homs_merge(self.homs+hm.homs)
        if len(merged_homs) != len(self.homs):
            return None
        new_hm.homs = merged_homs
        return new_hm

    def region_overlap(self, hm):
        """Checks if hm1 and hm2 overlap.  If they do, return a new HomMap object
        with the correct coordinates and an empty homs list, if not, return None."""
        new_hm = HomMap(None, None, None, [])
        if hm.chrom != self.chrom:
            return None
        else:
            new_hm.chrom = self.chrom
        if self.start >= hm.start and self.start <= hm.end:
            new_hm.start = hm.start
            new_hm.end = self.end
            return new_hm
        elif self.end >= hm.start and self.end <= hm.end:
            new_hm.start = self.start
            new_hm.end = hm.end
            return new_hm
        elif self.start <= hm.start and self.end >= hm.end:
            new_hm.start = self.start
            new_hm.end = self.end
            return new_hm
        elif self.start >= hm.start and self.end <= hm.end:
            new_hm.start = hm.start
            new_hm.end = hm.end
            return new_hm
        else:
            return None

    """
    Homs merging functions (original tuple recs)
    """

    def homs_merge(self, homs):
        """Given a list [(chrom, start, end)], sorts and merges the records
        such that they are increasing alphabetically by chromosome and internally
        by start, and none are directly adjacent or overlapping"""
        rec_l = homs
        if len(rec_l) <= 1:
            return rec_l
        rec_l = self.homs_sort(rec_l)
        indx = 0
        while True:
            merged = self.merge(rec_l[indx], rec_l[indx+1])
            if len(merged) == 2:
                indx += 1
            else:
                rec_l = rec_l[:indx]+merged+rec_l[indx+2:]
            if indx == len(rec_l)-1:
                break
        return rec_l
        
    def homs_sort(self, homs):
        """Takes a list [(chrom, start, end)] and returns the list sorted overall
        by alphabetically-increasing chrom and internally by increasing start"""
        rec_d = {}  #{chromosome: [(chromosome, start, end), ...]}
        #populate rec_d
        for rec in homs:
            if rec[0] in rec_d.keys():
                rec_d[rec[0]].append(rec)
            else:
                rec_d[rec[0]] = [rec]
        #internally sort every list in rec_d by start coordinate
        for chrom in rec_d.keys():
            rec_d[chrom] = sorted(rec_d[chrom], key=lambda x: x[1])
        #rebuild a list by sorting the keys (chromosomes) of rec_d
        out_l = []
        for chrom in sorted(rec_d.keys()):
            out_l += rec_d[chrom]
        return out_l
    
    def merge(self, rec1, rec2):
        """Checks if rec1 and rec2 overlap - checks if chromosomes match and
        if rec1.end >= rec2.start or rec2.end >= rec1.start"""
        if rec1[0] != rec2[0]:  #chrom1 != chrom2
            return [rec1, rec2]
        elif rec1[1] <= rec2[1] and rec1[2] >= rec2[1]:  #start1 <= start2 and end1 >= end2
            return [(rec1[0], rec1[1], rec2[2])]
        else:
            return [rec1, rec2]

class MergedMaps:
    """This class takes in a QnameMaps object and
    initlializes a set of HomMap objects with each individual entry
    from the QnameMaps object and performs merging according to the rules
    listed above"""
    def __init__(self, qm, filt_len=1000):
        """Saves the qm and calls HomMap loading functions and merging functions"""
        self.qm = qm
        self.filt_len = filt_len
        self.hms = []
        
        self.load_hms()
        self.merge()
        self.hm_filter()
        
    """
    HomMap loading
    """
    
    def load_hms(self):
        """Parses each entry of self.qm into a HomMap object
        and populates self.hms"""
        for qname in self.qm.keys():
            chrom, start, end = self.qm[qname]['input']
            if chrom != "7":  #Debug
                continue  #Debug
            homs = self.qm[qname]['homs']
            self.hms.append(HomMap(chrom, start, end, homs))
            
    """
    Merging
    """
    
    def merge(self):
        """Iteratively updates self.hms until it can no longer
        be reduced according to the 3 rules above"""
        self.hms = self.sort(self.hms)
        hms = self.hms[:]
        if len(hms) == 1:
            return hms
        indx = 0
        while True:
            new_hms = hms[indx].can_merge(hms[indx+1])
            if new_hms == None:
                indx += 1
            else:
                hms = hms[:indx]+[new_hms]+hms[indx+2:]
            if indx == len(hms)-1:
                break
        self.hms = hms
        
    def sort(self, hm_l):
        """Takes a list of HomMap objects and sorts them according to 
        alphabetical chromosome and then increasing start position"""
        hm_d = {}  #{chromosome: [HomMap]}
        #populate hm_d
        for hm in hm_l:
            if hm.chrom in hm_d.keys():
                hm_d[hm.chrom].append(hm)
            else:
                hm_d[hm.chrom] = [hm]
        #internally sort every list in hm_d by start coordinate
        for chrom in hm_d.keys():
            hm_d[chrom] = sorted(hm_d[chrom], key=lambda x: x.start)
        #rebuild a list by sorting the keys (chromosomes) of hm_d
        out_l = []
        for chrom in sorted(hm_d.keys()):
            out_l += hm_d[chrom]
        return out_l
    
    """
    Filtering
    """
    
    def hm_filter(self):
        """If any HMs are equal to or shorter than self.filt_len,
        remove them from the results"""
        self.hms = filter(lambda x: x.end-x.start > self.filt_len, self.hms)
    
    """
    Saving to bed file
    """
    
    def save(self, outpath):
        """Writes self.hms post-merge to a bed file"""
        outlines = []
        for hm in self.hms:
            outlines.append('\t'.join(map(str, [hm.chrom, hm.start, hm.end, len(hm.homs)]))+'\n')
        open(outpath, 'w').writelines(outlines)

if __name__ == "__main__":
    print("HomologMapping.py")
    qm = QnameMaps("PMS2/artificial_aligned.sam", "PMS2/PMS2_input.bed")
    mm = MergedMaps(qm, filt_len=1000)
    for hm in mm.hms:
        print(hm.chrom, hm.start, hm.end, hm.homs)
