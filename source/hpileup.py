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
import subprocess
import sys

#Repos
import tools.formats.Bed as Bed
import tools.io.configs as configs

#Local
import FakeFastq

class Hpileup:
    """Driver class for the pipeline, executes all analysis"""
    def __init__(self, inputpath, configpath):
        """Takes the input bed and configuration path with pipeline options in key=value
        format and executes the pipeline"""
        subprocess.call("mkdir data", shell=True)
        self.inputbed = Bed.Bed(inputpath)
        self.config = configs.read_config_kv(configpath)
        
        self.ffq = FakeFastq.FakeFastq(self.inputbed, self.config['reference'], 
                                        readlen=int(self.config['readlen']),
                                        overlap=float(self.config['overlap']))
        self.ffq.save("data/input_ref_reads.fq")

if __name__ == "__main__":
    print("===============")
    print("=   hpileup   =")
    print("===============\n")
    try:
        inputpath = sys.argv[1]
        configpath = sys.argv[2]
    except IndexError:
        print("Usage: python hpileup.py <input_bed_path> <config_path>")
        sys.exit(1)
    Hpileup(inputpath, configpath)
