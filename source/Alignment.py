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

#Local

SUPPORTED = "\n\t--Bowtie (bowtie)\n"
SUPPORTED += "\t--BWA-MEM (bwamem)\n"

class Alignment:
    """Wrapper class for calling different aligners in the
    pipeline"""
    def __init__(self, aligner_name, aligner, fastq_path, refpath, out_sam_path):
        """Takes the aligner config information, an input fastq, and a path
        for sam output, and runs the specified aligner on the input data and
        saves to the out_sam_path"""
        self.aligner_name = aligner_name
        self.aligner = aligner
        self.fastq_path = fastq_path
        self.refpath = refpath
        self.out_sam_path = out_sam_path
        
        print("Alignment: Preparing to align "+fastq_path+" to "+refpath)
        print("Alignment: Checking file locations...")
        self.check_files()
        print("Alignment: Reference files ready, preparing to call the aligner...")
        self.call_aligner()

    """
    Pre-processing
    """

    def check_files(self):
        """Checks if all the files specified for alignment are correct
        and in the right place"""
        if not os.path.isfile(self.aligner):  #check if user specified a real file for the aligner
            print("ERROR: Alignment: Aligner file not found at "+self.aligner)
            sys.exit(1)
        if self.aligner_name == "bowtie":  #user chose bowtie as their aligner
            print("Alignment: Chosen aligner: bowtie")
            if not os.path.isfile(self.refpath+".1.ebwt"):  #check if user has generated bowtie ref files
                print("ERROR: Alignment: Valid bowtie reference files not found at "+self.refpath)
                sys.exit(1)
        elif self.aligner_name == "bwamem":
            print("Alignment: Chosen aligner: bwamem")
            if not os.path.isfile(self.refpath):  #check if user specified a real file for the fasta
                print("ERROR: Alignment: Reference file "+self.refpath+" not found")
                sys.exit(1)
        else:
            print("Alignment: ERROR: Use of unsupported aligner "+self.aligner_name)
            print("\tPlease use one of the following alignment tools:")
            print(SUPPORTED)
            sys.exit(1)

    """
    Alignment Wrapping
    """

    def call_aligner(self):
        """Uses subprocess to make calls to the aligner based on specified aligner name"""
        if self.aligner_name == "bowtie":
            print("Alignment: Passing over execution to bowtie...")
            subprocess.call(self.aligner+" -a -S "+self.refpath+" "+self.fastq_path+" > "+self.out_sam_path, shell=True)
        elif self.aligner_name == "bwamem":
            print("Alignment: Passing over execution to BWA-MEM...")
            subprocess.call(self.aligner+" mem -a -t 2 "+self.refpath+" "+self.fastq_path+" > "+self.out_sam_path, shell=True)

if __name__ == "__main__":
    print("Alignment.py")
