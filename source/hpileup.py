"""
Welcome to the hpileup pipeline.

Developed by Alexander Wenzel under the mentorship of Dr. Vikas Bansal (UCSD)

contact: atwenzel@eng.ucsd.edu

The hpileup pipline leverages existing alignment and variant calling tools
to enable more informed variant calling in parts of the genomes with homologous
regions.

This script defines the driver file for the pipeline
"""

#Global
import os
import subprocess
import sys

#Repos
import tools.converters.sam2bed as sam2bed
import tools.formats.Bed as Bed
import tools.io.configs as configs

#Local
import Alignment
import FakeFastq

class Hpileup:
    """Driver class for the pipeline, executes all analysis"""
    def __init__(self, inputpath, configpath):
        """Takes the input bed and configuration path with pipeline options in key=value
        format and executes the pipeline"""
        if not os.path.exists("data/"):
            subprocess.call("mkdir data", shell=True)
        self.inputbed = Bed.Bed(inputpath)
        self.config = configs.read_config_kv(configpath)
        
        ###Pipeline workflow starts here###
        
        ##Build a fastq file with reads from the reference genome
        ##corresponding to the input bed file
        """self.ffq = FakeFastq.FakeFastq(self.inputbed, self.config['reference'], 
                                        readlen=int(self.config['readlen']),
                                        overlap=float(self.config['overlap']))
        self.ffq.save("data/input_ref_reads.fq")"""

        ##Run the user-specified aligner on the FakeFastq output
        self.aligner = Alignment.Alignment(self.config['aligner_name'], self.config['aligner'],
                                            "data/input_ref_reads.fq", self.config['aligner_ref'],
                                            "data/input_reads_aligned.sam")
        
        ##Convert the aligned file to a set of bed regions
        print("Hpileup: Converting data/input_reads_aligned.sam to bed format")
        s2b = sam2bed.Sam2Bed("data/input_reads_aligned.sam")

        print("Hpileup: Saving regions to data/input_reads_aligned.bed")
        s2b.save("data/input_reads_aligned.bed")

        print("Hpileup: Calling bedtools2 to sort data/input_reads_aligned.bed")
        subprocess.call(self.config['bedtools2']+"sortBed -i data/input_reads_aligned.bed > data/input_reads_aligned_sorted.bed", shell=True)

        print("Hpileup: Calling bedtools2 to merge data/input_reads_aligned.bed")
        subprocess.call(self.config['bedtools2']+"mergeBed -i data/input_reads_aligned_sorted.bed > data/input_reads_aligned.bed", shell=True)

        print("Hpileup: Cleaning up intermediate files...")
        subprocess.call("rm data/input_reads_aligned_sorted.bed", shell=True)
        

if __name__ == "__main__":
    print("\n===============")
    print("=   hpileup   =")
    print("===============\n")
    try:
        inputpath = sys.argv[1]
        configpath = sys.argv[2]
    except IndexError:
        print("Usage: python hpileup.py <input_bed_path> <config_path>")
        sys.exit(1)
    Hpileup(inputpath, configpath)
