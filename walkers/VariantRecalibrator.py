import sys, os, glob
import random
import subprocess
import re

from utils.config import CONFIG
from common import atoi, natural_keys
from common import slurm_header




def build_VariantRecalibrator_sbatch(working_dir,  variant_raw, type, scratch=False):
    """Builds the sbatch file in order to run VQSR
    
    :param str working_dir: directory where files will be created
    :param str variant_raw: vcf containing the raw variants
    :param str type: can be SNP or INDEL and specifies which options need to be used
    :param bool scratch: if True works on scratch
    
    :returns: path to the sbatch file
    """
    job_name           = "VQSR_{}".format(type)
    if type == "SNP":
        racal_file_name    = "{}_joincalled.snp.recal".format(CONFIG["output_header"])
        tranches_file_name = "{}_joincalled.snp.tranches".format(CONFIG["output_header"])
    else:
        racal_file_name    = "{}_joincalled.indel.recal".format(CONFIG["output_header"])
        tranches_file_name = "{}_joincalled.indel.tranches".format(CONFIG["output_header"])
    #create the sbatch file to merge all varaints or to copy the already single one
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as VariantRecalibrator:
        slurm = slurm_header(CONFIG["uppmax_project"], working_dir, job_name)
        VariantRecalibrator.write(slurm)
        VariantRecalibrator.write("\n")
        
        if scratch:
            VariantRecalibrator.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            VariantRecalibrator.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        
        GATK_input = "-input {} \\\n".format(variant_raw)
        if scratch:
            VariantRecalibrator.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(variant_raw, job_name))
            variant_raw_name = os.path.basename(variant_raw)
            GATK_input  = "-input $SNIC_TMP/{}/{} \\\n".format(job_name, variant_raw_name)
        
        GATK_command  = "java -Xmx64g -jar {} -T VariantRecalibrator  \\\n".format(CONFIG["GATK"])
        #add standard options
        for option in CONFIG["walkers"]["VariantRecalibrator"]:
            if isinstance(option, basestring):
                GATK_command += "{} \\\n".format(option)
        #now add specifc option for type
        added = False
        for option in CONFIG["walkers"]["VariantRecalibrator"]:
            if not isinstance(option, basestring) and type in option:
                specific_options = option[type]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"


        GATK_command += GATK_input
        if scratch:
            GATK_command += "-recalFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, racal_file_name)
            GATK_command += "-tranchesFile $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, tranches_file_name)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, racal_file_name , working_dir)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, tranches_file_name , working_dir)
        else:
            GATK_command += "-recalFile {}/VCF/{} \\\n".format(working_dir, racal_file_name)
            GATK_command += "-tranchesFile {}/VCF/{} \n\n".format(working_dir, tranches_file_name)

        VariantRecalibrator.write(GATK_command)
        #return path to sbach file
    return sbatch_file

    



def VariantRecalibrator():
    """Run VQSR
    
    :returns: the sbatch_file to be executed
    """
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "04_SelectVariants")):
        sys.exit("Directory 04_SelectVariants does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "05_VariantRecalibrator")):
        print "WARNING: 05_VariantRecalibrator already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "05_VariantRecalibrator"))
        os.mkdir(os.path.join(cwd, "05_VariantRecalibrator", "sbatch"))
        os.mkdir(os.path.join(cwd, "05_VariantRecalibrator", "std_err"))
        os.mkdir(os.path.join(cwd, "05_VariantRecalibrator", "std_out"))
        os.mkdir(os.path.join(cwd, "05_VariantRecalibrator", "VCF"))
    #Build the sbatch files for merging step
    working_dir       = os.path.join(cwd, "05_VariantRecalibrator")
    variant_raw_snp   = os.path.join(cwd, "04_SelectVariants", "VCF", "{}_joincalled.snp.g.vcf.gz".format(CONFIG["output_header"]))
    variant_raw_indel = os.path.join(cwd, "04_SelectVariants", "VCF", "{}_joincalled.indel.g.vcf.gz".format(CONFIG["output_header"]))
    sbatch_file = build_VariantRecalibrator_sbatch(working_dir, variant_raw_snp, "SNP", CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    sbatch_file = build_VariantRecalibrator_sbatch(working_dir, variant_raw_indel, "INDEL", CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    return sbatch_files

