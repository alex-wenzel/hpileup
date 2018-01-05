"""
This script defines a class to call GATK on a SAM file
assuming that different regions of the file have different ploidies.
The class uses ploidy information from previous pipeline steps
"""

#Global
import subprocess

#Repos
from tools.formats.Bed import Bed
from tools.io.Log import Log

#Local

class SegGATK:
    """Wrapper around GATK to call it incrementally on different
    regions of the same SAM file that have different ploidies"""
    def __init__(self, sampath, args):
        """Saves sampath and args (argparse object from
        hpileup) and starts processing on each region"""
        self.sampath = sampath
        self.args = args
        self.ploidy_bed = Bed(self.args.outdir+"ploidy.bed")
        self.l = Log()

        self.iterate_gatk()

    """
    GATK iteration
    """

    def iterate_gatk(self):
        """Calls GATK successively on the input sam file,
        using each region and ploidy from self.ploidy"""
        gatk = self.args.gatk
        ref = self.args.ref
        self.l.log("SegGATK: Calling GATK on all ploidy regions for "+self.sampath+"...")
        for line in self.ploidy_bed:
            ploidy = str(2+int(line.data[0])*2)
            reg_str = '-'.join([line.chromosome, str(line.start), str(line.end)])
            outpath = '.'.join(self.sampath.split('.')[:-1])+"_"+reg_str+".vcf"
            cmd = "java -jar "+gatk+" -T UnifiedGenotyper "
            cmd += "-R "+ref+" -I "+self.sampath+" -o "+outpath
            cmd += " -glm BOTH -L "
            cmd += line.chromosome+":"+str(line.start)+"-"+str(line.end)
            cmd += " -ploidy "+ploidy
            self.l.log("Calling GATK on "+self.sampath+" region "+reg_str+" with ploidy "+ploidy)
            self.l.log('\t'+cmd)
            subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    print("SegGATK.py")
