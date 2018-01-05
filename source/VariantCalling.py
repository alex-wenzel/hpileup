"""
This script implements the variant calling portion of the pipeline
For reference, see variant_calling.ipynb
"""

#Global
import subprocess

#Repos
import tools.io.file_iterator as fi
from tools.io.Log import Log

#Local
from SegGATK import SegGATK

class VariantCalling:
    """This class implements the variant calling stage of the
    pipeline"""
    def __init__(self, args):
        """Saves args (argparse object from hpileup) and executes
        the variant calling pipeline steps"""
        self.args = args
        self.l = Log()

        self.run()
        
    """
    Pipeline driver
    """

    def run(self):
        """Pipeline driver function"""
        self.l.log("VariantCalling: Calling variants for all samples...")
        for sample_path in self.args.samples:
            outbase = ".".join(sample_path.split('.')[:-1])
            self.l.log("VariantCalling: Calling variants for "+outbase+"...")

            ##set all mapping qualities to 60
            self.reset_mapq(outbase)

            ##convert sam to bam
            self.sam2bam(outbase)

            ##add read groups
            self.add_read_groups(outbase)
            
            ##sort + index
            self.sort_index(outbase)

            ##call incremental GATK
            SegGATK(outbase+"_reset-mapq_rg_sorted.bam", self.args)

    """
    Reset mapping quality
    """
    
    def reset_mapq(self, outbase):
        """Takes the location of an input SAM file
        and iterates through it, setting all mapping quality
        scores to 60"""
        outlines = []
        self.l.log("VariantCalling: Setting mapping qualities for "+outbase+" to 60...")
        for line in fi.iterate(open(outbase+'_realigned_input.sam')):
            if line[0] == "@":
                outlines.append(line)
                continue
            linevals = line.strip('\n').split('\t')
            linevals[4] = "60"
            if linevals[1] == "256":
                linevals[1] = "0"
            elif linevals[1] == "272":
                linevals[1] = "16"
            outlines.append('\t'.join(linevals)+'\n')
        open(outbase+"_reset-mapq.sam", 'w').writelines(outlines)

    """
    SAM to BAM
    """

    def sam2bam(self, outbase):
        """Converts the input sam filed to compressed bam format"""
        input_sam = outbase+"_reset-mapq.sam"
        output_bam = outbase+"_reset-mapq.bam"
        samtools = self.args.samtools_loc
        cmd = samtools+" view -S -b "+input_sam+" > "+output_bam
        self.l.log("VariantCalling: Converting "+outbase+" to compressed BAM format...")
        self.l.log("\t"+cmd)
        subprocess.call(cmd, shell=True)

    """
    Picard Read Groups
    """

    def add_read_groups(self, outbase):
        """Calls the Picard jar to add read groups to the
        input SAM (necessary for sorting/indexing"""
        input_bam = outbase+"_reset-mapq.bam"
        output_bam = outbase+"_reset-mapq_rg.bam"
        picard = self.args.picard
        cmd = "java -jar "+picard+" AddOrReplaceReadGroups "
        cmd += "I="+input_bam
        cmd += " O="+output_bam
        cmd += " RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 "
        cmd += "VALIDATION_STRINGENCY=LENIENT"
        self.l.log("VariantCalling: Adding read groups with the following command...")
        self.l.log("\t"+cmd)
        subprocess.call(cmd, shell=True)

    """
    Sort/index bam
    """

    def sort_index(self, outbase):
        """Calls samtools sort and samtools index on the input bam"""
        input_bam = outbase+"_reset-mapq_rg.bam"
        output_bam = outbase+"_reset-mapq_rg_sorted.bam"
        samtools = self.args.samtools_loc

        cmd = samtools+" sort "+input_bam+" -o "+output_bam
        self.l.log("VariantCalling: Sorting bam with the following command...")
        self.l.log("\t"+cmd)
        subprocess.call(cmd, shell=True)

        input_bam = output_bam
        cmd = samtools+" index "+input_bam
        self.l.log("VariantCalling: Indexing bam with the following command...")
        self.l.log("\t"+cmd)
        subprocess.call(cmd, shell=True)

    """
    VCF Merging
    """

    def merge_vcf(self, outbase):
        """Following multi-ploidy calls from SegGATK, merge
        all VCFs into a single output VCF with variants from
        regions with multiple ploidies (depends on implementation
        of SegGATK)"""


if __name__ == "__main__":
    print("VariantCalling.py")
