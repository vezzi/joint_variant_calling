import sys, os, glob
import random
import subprocess
import re

from utils.config import CONFIG
from common import atoi, natural_keys
from common import slurm_header




def build_SelectVariants_sbatch(working_dir, variant_file, scratch=False):
    """Builds the sbatch file in order to combine genomics.vcf divided up in chr into a single one
    
    :param str working_dir: directory where files will be created
    :param str variants_dir: directory where the vcf to be merged are present
    :param bool scratch: if True works on scratch
    
    :returns: path to the sbatch file
    """
    job_name       = "SelectVariants"
    output_file_snp        = "{}_joincalled.snp.g.vcf.gz".format(CONFIG["output_header"])
    output_file_snp_eval   = "{}_joincalled.snp.eval".format(CONFIG["output_header"])
    output_file_indel      = "{}_joincalled.indel.g.vcf.gz".format(CONFIG["output_header"])
    output_file_indel_eval = "{}_joincalled.indel.eval".format(CONFIG["output_header"])
    #create the sbatch file to merge all varaints or to copy the already single one
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as SelectVariants:
        slurm = slurm_header(CONFIG["uppmax_project"], working_dir, job_name)
        SelectVariants.write(slurm)
        SelectVariants.write("\n")
        
        if scratch:
            SelectVariants.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            SelectVariants.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        
        GATK_input = "-V {} \\\n".format(variant_file)
        if scratch:
            SelectVariants.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(variant_file, job_name))
            variant_file_name = os.path.basename(variant_file)
            GATK_input  = "-V $SNIC_TMP/{}/{} \\\n".format(job_name, variant_file_name)
        
        GATK_command  = "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
        for option in CONFIG["walkers"]["SelectVariants"]:
            GATK_command += "{} \\\n".format(option)
        GATK_command += GATK_input
        #create command for SNPs
        GATK_command_snp = GATK_command
        GATK_command_snp += "-selectType SNP \\\n"
        #create command for indels
        GATK_command_indel = GATK_command
        GATK_command_indel += "-selectType INDEL \\\n"
        if scratch:
            GATK_command_snp +=  "-o $SNIC_TMP/{}/VCF/{}\n".format(job_name, output_file_snp)
            GATK_command_snp += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/\n".format(job_name, output_file_snp , working_dir)
            GATK_command_indel +=  "-o $SNIC_TMP/{}/VCF/{}\n".format(job_name, output_file_indel)
            GATK_command_indel += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/\n".format(job_name, output_file_indel , working_dir)
        else:
            GATK_command_snp   += "-o {}/VCF/{}\n\n".format(working_dir, output_file_snp)
            GATK_command_indel += "-o {}/VCF/{}\n\n".format(working_dir, output_file_indel)
        SelectVariants.write(GATK_command_snp)
        SelectVariants.write("\n\n")
        SelectVariants.write(GATK_command_indel)
        #now we can tun EVAL
        GATK_command = "java -Xmx250g -jar {} -T VariantEval -nt 16 \\\n".format(CONFIG["GATK"])
        for option in CONFIG["walkers"]["VariantEval"]:
            GATK_command += "{} \\\n".format(option)
        GATK_command_snp = GATK_command + "--eval {}/VCF/{} \\\n".format(working_dir, output_file_snp)
        GATK_command_snp += "-o {}/VCF/{} \n".format(working_dir, output_file_snp_eval)
        GATK_command_indel = GATK_command + "--eval {}/VCF/{} \\\n".format(working_dir, output_file_indel)
        GATK_command_indel += "-o {}/VCF/{} \n".format(working_dir, output_file_indel_eval)
        SelectVariants.write(GATK_command_snp)
        SelectVariants.write("\n\n")
        SelectVariants.write(GATK_command_indel)
        #return path to sbach file
    return sbatch_file

    



def SelectVariants():
    """Selects variants based on type, e.g. SNPs or Indels and performs also Evaluation of them
    
    :returns: the sbatch_file to be executed
    """
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "03_CatVariants")):
        sys.exit("Directory 03_CatVariants does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "04_SelectVariants")):
        print "WARNING: 04_SelectVariants already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "04_SelectVariants"))
        os.mkdir(os.path.join(cwd, "04_SelectVariants", "sbatch"))
        os.mkdir(os.path.join(cwd, "04_SelectVariants", "std_err"))
        os.mkdir(os.path.join(cwd, "04_SelectVariants", "std_out"))
        os.mkdir(os.path.join(cwd, "04_SelectVariants", "VCF"))
    #Build the sbatch files for merging step
    working_dir  = os.path.join(cwd, "04_SelectVariants")
    variant_file   = os.path.join(cwd, "03_CatVariants", "VCF", "{}_joincalled.g.vcf.gz".format(CONFIG["output_header"]))
    sbatch_file = build_SelectVariants_sbatch(working_dir, variant_file, CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    return sbatch_files

