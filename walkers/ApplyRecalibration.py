import sys, os, glob
import random
import subprocess
import re

from utils.config import CONFIG
from common import atoi, natural_keys
from common import slurm_header




def build_ApplyRecalibration_sbatch(working_dir, variant_raw,  recal, tranches, type, scratch=False):
    """Builds the sbatch file in order to run VQSR
    
    :param str working_dir: directory where files will be created
    :param str variant_raw: vcf containing the raw variants
    :param str type: can be SNP or INDEL and specifies which options need to be used
    :param bool scratch: if True works on scratch
    
    :returns: path to the sbatch file
    """
    
    job_name           = "ApplyRecalibration_{}".format(type)
    if type == "SNP":
        output_file_name    = "{}_joincalled.snp.recalibrated.filtered.vcf.gz".format(CONFIG["output_header"])
    else:
        output_file_name    = "{}_joincalled.indel.recalibrated.filtered.vcf.gz".format(CONFIG["output_header"])
    #create the sbatch file to merge all varaints or to copy the already single one
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as ApplyRecalibration:
        slurm = slurm_header(CONFIG["uppmax_project"], working_dir, job_name)
        ApplyRecalibration.write(slurm)
        ApplyRecalibration.write("\n")
        
        if scratch:
            ApplyRecalibration.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            ApplyRecalibration.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        
        GATK_command  = "java -Xmx64g -jar {} -T ApplyRecalibration  \\\n".format(CONFIG["GATK"])
        GATK_input = "-input {} \\\n".format(variant_raw)
        GATK_recal = "-recalFile {} \\\n".format(recal)
        GATK_tranches = "-tranchesFile {} \\\n".format(tranches)
        if scratch:
            ApplyRecalibration.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(variant_raw, job_name))
            variant_raw_name = os.path.basename(variant_raw)
            ApplyRecalibration.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(recal, job_name))
            recal_name = os.path.basename(recal)
            ApplyRecalibration.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(tranches, job_name))
            tranches_name = os.path.basename(tranches)
            GATK_input    = "-input $SNIC_TMP/{}/{} \\\n".format(job_name, variant_raw_name)
            GATK_recal    = "-recalFile $SNIC_TMP/{}/{} \\\n".format(job_name ,recal_name)
            GATK_tranches = "-tranchesFile $SNIC_TMP/{}/{} \\\n".format(job_name, tranches_name)
        GATK_command += GATK_input
        GATK_command += GATK_recal
        GATK_command += GATK_tranches
        #add standard options
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if isinstance(option, basestring):
                GATK_command += "{} \\\n".format(option)
        #now add specifc option for type
        added = False
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if not isinstance(option, basestring) and type in option:
                specific_options = option[type]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"

        if scratch:
            GATK_command += "-o $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, output_file_name)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, output_file_name , working_dir)
        else:
            GATK_command += "-o {}/VCF/{} \n\n".format(working_dir, output_file_name)

        ApplyRecalibration.write(GATK_command)
        #return path to sbach file
    return sbatch_file

    



def ApplyRecalibration():
    """Run VQSR -- ApplyRecalibration
    
    :returns: the sbatch_file to be executed
    """
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "05_VariantRecalibrator")):
        sys.exit("Directory 05_VariantRecalibrator does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "06_ApplyRecalibration")):
        print "WARNING: 06_ApplyRecalibration already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "06_ApplyRecalibration"))
        os.mkdir(os.path.join(cwd, "06_ApplyRecalibration", "sbatch"))
        os.mkdir(os.path.join(cwd, "06_ApplyRecalibration", "std_err"))
        os.mkdir(os.path.join(cwd, "06_ApplyRecalibration", "std_out"))
        os.mkdir(os.path.join(cwd, "06_ApplyRecalibration", "VCF"))
    #Build the sbatch files for merging step
    working_dir      = os.path.join(cwd, "06_ApplyRecalibration")
    variant_raw_snp  = os.path.join(cwd, "04_SelectVariants", "VCF", "{}_joincalled.snp.g.vcf.gz".format(CONFIG["output_header"]))
    recal_snp        = os.path.join(cwd, "05_VariantRecalibrator", "VCF", "{}_joincalled.snp.recal".format(CONFIG["output_header"]))
    tranches_snp     = os.path.join(cwd, "05_VariantRecalibrator", "VCF", "{}_joincalled.snp.tranches".format(CONFIG["output_header"]))
    #now indels
    variant_raw_indel  = os.path.join(cwd, "04_SelectVariants", "VCF", "{}_joincalled.indel.g.vcf.gz".format(CONFIG["output_header"]))
    recal_indel        = os.path.join(cwd, "05_VariantRecalibrator", "VCF", "{}_joincalled.indel.recal".format(CONFIG["output_header"]))
    tranches_indel     = os.path.join(cwd, "05_VariantRecalibrator", "VCF", "{}_joincalled.indel.tranches".format(CONFIG["output_header"]))


    sbatch_file = build_ApplyRecalibration_sbatch(working_dir, variant_raw_snp,  recal_snp, tranches_snp, "SNP", CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    sbatch_file = build_ApplyRecalibration_sbatch(working_dir, variant_raw_indel,  recal_indel, tranches_indel, "INDEL", CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    return sbatch_files

