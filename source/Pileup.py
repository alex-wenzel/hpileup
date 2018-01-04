"""
This script defines a class to run the pileup section of the
pipeline.  See pileup2.ipynb for prototype/reference
"""

#Global
import subprocess

#Repos
from tools.io.Log import Log

#Local
from Alignment import Alignment
from SamRegionFilter import SamRegionFilter as SRF

class Pileup:
    """Takes any number of SAM/BAM files and uses samtools and Bowtie2
    to collapse all reads from various homologous regions to the
    user supplied input regions"""
    def __init__(self, args):
        """Iterates through the args.samples (sam/bam files) and implements collapsing
        of all reads from homologous regions onto the user input regions.
        args is the argparse object from hpileup.py"""
        self.args = args
        self.l = Log()
        
        self.run()

    """
    Pipeline driver
    """

    def run(self):
        """Calls the pipeline steps for each file in self.sam_paths"""
        self.l.log("Piling up reads from the input SAM/BAM files...")
        for sample_path in self.args.samples:
            self.l.log("Pileup: Processing "+sample_path+"...", tabs=1)
            ##convert from bam to sam if needed
            sample_path = self.check_for_bam(sample_path)
            outbase = ".".join(sample_path.split('.')[:-1])

            ##Filter out regions that are not in input or homologous regions
            self.keep_input_homolog(outbase)

            ##Revert to FASTQ
            self.revert(outbase)

            ##Realign filtered to all
            self.realign(outbase)

            ##Filter to only input regions
            self.keep_input(outbase)

    """
    BAM -> SAM
    """

    def check_for_bam(self, sambam_path):
        """Checks if the input file is a bam file.  If so, convert it and return
        a new path, otherwise return the input path"""
        if sambam_path.split('.')[-1] == ".bam":
            self.l.log("Converting "+sample_path+" to uncompressed SAM format...", tabs=2)
            return bam2sam(sambam_path)
        return sambam_path

    def bam2sam(self, bam_path):
        """Converts a BAM file to SAM. Saves the result in the same
        directory as the BAM input.  Returns the path to the new SAM file"""
        outpath = '.'.join(bam_path.split('.')[:-1])+'.sam'
        samtools = self.args.samtools_loc
        cmd = samtools+" view -h -o "+outpath+" "+bam_path
        subprocess.call(cmd, shell=True)
        return outpath

    """
    Filtering functions
    """

    def keep_input_homolog(self, outbase):
        """Takes a path to a sam file and filters it so that only the regions
        in the input and homologous regions are kept"""
        sampath = outbase+".sam"
        self.l.log("Filtering "+sampath+" to contain reads in input or homlogs...", tabs=2)
        input_hom_bed_path = self.args.outdir+"input_homolog.bed"
        srf = SRF(sampath, input_hom_bed_path)
        srf.save(outbase+"_input-homolog.sam")

    def keep_input(self, outbase):
        """Takes a path to a sam file and filters it so that only the regions
        in the input regions are kept"""
        sampath = outbase+"_input-homolog_realigned.sam"
        self.l.log("Filtering "+sampath+" to contain reads in input regions only...", tabs=2)
        input_bed_path = self.args.outdir+"input.bed"
        srf = SRF(sampath, input_bed_path)
        srf.save(outbase+"_realigned_input.sam")

    """
    Revert SAM --> FASTQ
    """

    def revert(self, outbase):
        """Calls samtools collate and fastq to revert sampath to a 
        fastq file saved using the base prefix outbase"""
        samtools = self.args.samtools_loc
        sampath = outbase+"_input-homolog.sam"
        collated_out = outbase+"_input-homolog_collated"
        fastq_out = outbase+"_input-homolog.fq"
        self.l.log("Reverting "+sampath+" to FASTQ...", tabs=2)

        cmd = samtools+" collate -u -n 1 -l 1 --output-fmt SAM "+sampath+" "+collated_out
        self.l.log(cmd, tabs=3)
        subprocess.call(cmd, shell=True)

        cmd = samtools+" fastq "+collated_out+".sam > "+fastq_out
        self.l.log(cmd, tabs=3)
        print("")
        subprocess.call(cmd, shell=True)
        print("")

    """
    Alignment wrapper
    """

    def realign(self, outbase):
        """Calls the Alignment class, which calls Bowtie2 to realign reads"""
        input_fq = outbase+"_input-homolog.fq"
        output_sam = outbase+"_input-homolog_realigned.sam"
        self.l.log("Aligning "+input_fq+ " to the genome (output at "+output_sam+")...", tabs=2)
        Alignment(self.args, input_fq, output_sam)

if __name__ == "__main__":
    print("Pileup.py")
