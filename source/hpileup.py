"""
This file defines driver code for the hpileup pipeline for variant calling
on regions with complex genetic homologies

Developed by Alexander Wenzel under the mentorship of Dr. Vikas Bansal (UCSD)

contact: atwenzel@eng.ucsd.edu

The hpileup pipline leverages existing alignment and variant calling tools
to enable more informed variant calling in parts of the genomes with homologous
regions.

This script defines the driver file for the pipeline
"""

#Global
import argparse
import os
import subprocess
import sys

#Repos
import tools.converters.sam2bed as sam2bed
import tools.formats.Bed as Bed
from tools.io.Log import Log

#Local
import Alignment
import FakeFastq
import FastaSubset
import HomologMapping

class Hpileup:
    """Driver class for the pipeline, executes all analysis"""
    #def __init__(self, inputpath, configpath):
    def __init__(self, args):
        """Takes the argparse object from main and executes the full pipeline"""
        self.l = Log()
        self.l.log("Loading "+args.input+"...")
        self.input_bed = Bed.Bed(args.input)
        
        ##create the output directory, if it does not exist
        self.l.log("Checking output directory...")
        if not os.path.exists(args.outdir):
            self.l.log("Directory "+args.outdir+" does not exist, making it...", tabs=1)
            subprocess.call("mkdir "+args.outdir, shell=True)
        if args.outdir[-1] != "/":
            args.outdir = args.outdir+"/"

        ##Generate FakeFastq file
        #self.l.log("Generating the artificial FASTQ file for "+args.input+"...")
        #ffq = FakeFastq.FakeFastq(self.input_bed, args.ref)
        #ffq.save(args.outdir+"artificial_reads.fq")

        ##Run the Bowtie 2 aligner on the artificial reads
        #Alignment.Alignment(args, args.outdir+"artificial_reads.fq", 
        #                    args.outdir+"artificial_aligned.sam")

        ##Run the HomologMapping scripts
        self.l.log("Compiling all homologous reads...")
        qm = HomologMapping.QnameMaps(args.outdir+"artificial_aligned.sam", args.input)
        self.l.log("Merging homologous reads...")
        mm = HomologMapping.MergedMaps(qm, filt_len=1000)
        self.l.log("Saving ploidy info to "+args.outdir+"ploidy.bed")
        mm.save(args.outdir+"ploidy.bed")

if __name__ == "__main__":
    print("\n===============")
    print("=   hpileup   =")
    print("===============\n")
    
    p = argparse.ArgumentParser(description="hpileup: A pipeline for variant calling in segdups")

    p.add_argument("-i", "--input", required=True, 
                    help="The path to a BED-format file with regions of interest for variant calling")
    p.add_argument("-r", "--ref", required=True,
                    help="The path to the reference genome FASTA (not pre-generated Bowtie2 files)")
    p.add_argument("-s", "--samples", nargs='+', required=True,
                    help="Paths to one or more SAM/BAM files (BAMs will be converted to SAMs)")
    p.add_argument("-g", "--gatk", required=True,
                    help="Path to the GenomeAnalysisTK jar")
    p.add_argument("--outdir", default="./", help="The directory to use for results")
    p.add_argument("--bowtie2_loc", default="bowtie2",
                    help="If Bowtie2 is not in your PATH, use this option to specify its location")
    p.add_argument("--bowtie2_ref", required=True,
                    help="The location of the Bowtie2 reference files")
    p.add_argument("--samtools_loc", default="samtools",
                    help="If samtools is not in your PATH, use this option to specify its location")
    p.add_argument("--threads", type=int, default=1,
                    help="The number of threads to use for multi-threaded components (Bowtie2 and GATK)")
    
    args = p.parse_args()
    Hpileup(args)
