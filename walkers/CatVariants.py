import sys, os, glob
import random
import subprocess
import re

from utils.config import CONFIG
from common import atoi, natural_keys
from common import slurm_header




def build_CatVariants_sbatch(working_dir, variants_dir,  scratch=False):
    """Builds the sbatch file in order to combine genomics.vcf divided up in chr into a single one
    
    :param str working_dir: directory where files will be created
    :param str variants_dir: directory where the vcf to be merged are present
    :param bool scratch: if True works on scratch
    
    :returns: path to the sbatch file
    """
    job_name      = "CatVariants"
    output_file   = "{}_joincalled.g.vcf.gz".format(CONFIG["output_header"])
    #create the sbatch file to merge all varaints or to copy the already single one
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as CatVariants:
        slurm = slurm_header(CONFIG["uppmax_project"], working_dir, job_name)
        CatVariants.write(slurm)
        CatVariants.write("\n")
        
        if len(CONFIG["intervals_list"]) == 0:
            #in this case I need only to copy the already single file
            source = os.path.join(variants_dir, "{}_joincalled.g.vcf.gz".format(CONFIG["output_header"]) )
            dest = os.path.join(working_dir, "VCF",  "{}_joincalled.g.vcf.gz".format(CONFIG["output_header"]) )
            CatVariants.write("cp {} {}\n".format(source, dest))
        else:
            if scratch:
                CatVariants.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
                CatVariants.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
            #now cycle over the intervals and build the GATK command
            catvariants_string_input = ""
            # this sorts intervals created given that they have some number in their name specifing the order
            CONFIG["intervals_list"].sort(key=natural_keys)
            for interval in CONFIG["intervals_list"]:
                interval_name = os.path.basename(interval).split(".")[0]
                vcf_interval = os.path.join(variants_dir, "{}_joincalled_{}.g.vcf.gz".format(CONFIG["output_header"], interval_name))
                if scratch:
                    CatVariants.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(vcf_interval, job_name))
                    vcf_interval_name =  os.path.basename(vcf_interval)
                    vcf_interval = "$SNIC_TMP/{}/{}".format(job_name, vcf_interval_name)
                catvariants_string_input += "-V {} \\\n".format(vcf_interval)
        
            GATK_command = "java -cp {} org.broadinstitute.gatk.tools.CatVariants \\\n".format(CONFIG["GATK"])
            for option in CONFIG["walkers"]["CatVariants"]:
                GATK_command += "{} \\\n".format(option)
            GATK_command += "{} ".format(catvariants_string_input)
            if scratch:
                GATK_command +=  "-out $SNIC_TMP/{}/VCF/{}\n".format(job_name, output_file)
                #once this is done rsync back to lupus
                GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/\n".format(job_name, output_file , working_dir)
            else:
                GATK_command += "-out {}/VCF/{}\n\n".format(working_dir, output_file)
            CatVariants.write(GATK_command)
         #return path to sbach file
    return sbatch_file

    



def CatVariants():
    """Merges variants divided into a chr by chr fashion into a single one.
    :param list intervals: intervals to be used in case we want to split analysis chr by chr
    
    :returns: the sbatch_file to be executed
    """
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "02_GenotypeGVCFs")):
        sys.exit("Directory 02_GenotypeGVCFs does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "03_CatVariants")):
        print "WARNING: 03_CatVariants already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "03_CatVariants"))
        os.mkdir(os.path.join(cwd, "03_CatVariants", "sbatch"))
        os.mkdir(os.path.join(cwd, "03_CatVariants", "std_err"))
        os.mkdir(os.path.join(cwd, "03_CatVariants", "std_out"))
        os.mkdir(os.path.join(cwd, "03_CatVariants", "VCF"))
    #Build the sbatch files for merging step
    working_dir  = os.path.join(cwd, "03_CatVariants")
    variants_dir = os.path.join(cwd, "02_GenotypeGVCFs", "VCF")
    sbatch_file = build_CatVariants_sbatch(working_dir, variants_dir, CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    return sbatch_files

