"""
A wrapper class to accomodate different genome sequence aligners.

Aligners currently supported (config file aligner name):
    --Bowtie (bowtie)
    --BWA-MEM (bwamem)
"""

#Global
import os
import subprocess
import sys

#Repos
from tools.io.Log import Log

#Local

class Alignment:
    """Wrapper class for calling different aligners in the
    pipeline"""
    def __init__(self, args, fastq_path, out_sam_path):
        """Takes the argparse object from Hpileup (for bowtie2 configs, and an
        input fastq path and output sam path, and performs bowtie2 alignment.
        Hardcoded to use the -k 10 alignment option"""
        self.args = args
        self.fastq_path = fastq_path
        self.out_sam_path = out_sam_path
        self.l = Log()
        
        self.l.log("Alignment: Preparing to align "+self.fastq_path+" to "+args.bowtie2_ref)
        self.l.log("Alignment: Checking file locations...")
        self.check_files()
        self.l.log("Alignment: Reference files ready, preparing to call Bowtie...")
        self.call_aligner()

    """
    Pre-processing
    """

    def check_files(self):
        """Checks if all the files specified for alignment are correct
        and in the right place"""
        if not os.path.isfile(self.args.bowtie2_ref+".1.bt2"):  #check if user has generated bowtie ref files
            self.l.error("Alignment: Valid bowtie reference files not found at "+self.args.bowtie2_ref,
                            die=True, code=1)

    """
    Alignment Wrapping
    """

    def call_aligner(self):
        """Uses subprocess to make calls to the aligner based on specified aligner name"""
        cmd = self.args.bowtie2_loc+" -x "+self.args.bowtie2_ref
        cmd += " -p "+str(self.args.threads)+" -k 10 -U "
        cmd += self.fastq_path+" -S "+self.out_sam_path
        self.l.log("Calling Bowtie 2 with the following command...\n\t"+cmd)
        subprocess.call(cmd, shell=True)
        self.l.log("Bowtie 2 finished")

if __name__ == "__main__":
    print("Alignment.py")
