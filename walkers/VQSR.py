import sys, os, glob
import random
import subprocess
import re

from utils.config import CONFIG
from common import atoi, natural_keys
from common import slurm_header


#### this applies VQSR as explained in
#### http://gatkforums.broadinstitute.org/gatk/discussion/2805/howto-recalibrate-variant-quality-scores-run-vqsr
#### solves problem with mixed positions


def build_VQSR_sbatch(working_dir,  variant_raw, scratch=False):
    """Builds the sbatch file in order to run VQSR
    
    :param str working_dir: directory where files will be created
    :param str variant_raw: vcf containing the raw variants
    :param bool scratch: if True works on scratch
    
    :returns: path to the sbatch file
    """
    
    job_name           = "VQSR"
    #first build the model for SNPS
    racal_file_name_snps    = "{}_joincalled.snp.recal".format(CONFIG["output_header"])
    tranches_file_name_snps = "{}_joincalled.snp.tranches".format(CONFIG["output_header"])
    #apply the model to SNPs only
    variant_recal_snp_raw_indels    = "{}_joincalled.recal_snp_raw_indels.vcf.gz".format(CONFIG["output_header"])
    #and then build the model for INDELS
    racal_file_name_indels    = "{}_joincalled.indel.recal".format(CONFIG["output_header"])
    tranches_file_name_indels = "{}_joincalled.indel.tranches".format(CONFIG["output_header"])
    variant_recal_snp_recal_indels    = "{}_joincalled.recal_snp_recal_indels.vcf.gz".format(CONFIG["output_header"])
    #create the sbatch file to merge all varaints or to copy the already single one
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as VQSR:
        slurm = slurm_header(CONFIG["uppmax_project"], working_dir, job_name)
        VQSR.write(slurm)
        VQSR.write("\n")
        ##############################################
        #### compute recalibration tables for SNPs ###
        ##############################################
        if scratch:
            VQSR.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            VQSR.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        GATK_input = "-input {} \\\n".format(variant_raw)
        if scratch:
            VQSR.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(variant_raw, job_name))
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
            if not isinstance(option, basestring) and "SNP" in option:
                specific_options = option["SNP"]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"
        GATK_command += GATK_input
        if scratch:
            GATK_command += "-recalFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, racal_file_name_snps)
            GATK_command += "-tranchesFile $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, tranches_file_name_snps)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, racal_file_name_snps , working_dir)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, tranches_file_name_snps , working_dir)
        else:
            GATK_command += "-recalFile {}/VCF/{} \\\n".format(working_dir, racal_file_name_snps)
            GATK_command += "-tranchesFile {}/VCF/{} \\\n".format(working_dir, tranches_file_name_snps)
        VQSR.write(GATK_command)
        VQSR.write("\n")
        ##########################################
        ##### now apply recalibration for SNPs ###
        ##########################################
        GATK_command  = "java -Xmx64g -jar {} -T ApplyRecalibration  \\\n".format(CONFIG["GATK"])
        #### GATK_input is the same
        if scratch:
            GATK_command += "-recalFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, racal_file_name_snps)
            GATK_command += "-tranchesFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, tranches_file_name_snps)
        else:
            GATK_command += "-recalFile {}/VCF/{} \\\n".format(working_dir, racal_file_name_snps)
            GATK_command += "-tranchesFile {}/VCF/{} \n\n".format(working_dir, tranches_file_name_snps)
        GATK_command += GATK_input
        #add standard options
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if isinstance(option, basestring):
                GATK_command += "{} \\\n".format(option)
        #now add specifc option for type
        added = False
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if not isinstance(option, basestring) and "SNP" in option:
                specific_options = option["SNP"]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"

        if scratch:
            GATK_command += "-o $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, variant_recal_snp_raw_indels)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, variant_recal_snp_raw_indels , working_dir)
        else:
            GATK_command += "-o {}/VCF/{} \n\n".format(working_dir, variant_recal_snp_raw_indels)
        VQSR.write(GATK_command)
        VQSR.write("\n")
        ################################################
        #### compute recalibration tables for INDELS ###
        ################################################
        GATK_input = "-input {}/VCF/{} \n\n".format(working_dir, variant_recal_snp_raw_indels)
        if scratch:
            GATK_input  = "-input $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, variant_recal_snp_raw_indels)
        GATK_command  = "java -Xmx64g -jar {} -T VariantRecalibrator  \\\n".format(CONFIG["GATK"])
        #add standard options
        for option in CONFIG["walkers"]["VariantRecalibrator"]:
            if isinstance(option, basestring):
                GATK_command += "{} \\\n".format(option)
        #now add specifc option for type
        added = False
        for option in CONFIG["walkers"]["VariantRecalibrator"]:
            if not isinstance(option, basestring) and "INDEL" in option:
                specific_options = option["INDEL"]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"
        GATK_command += GATK_input
        if scratch:
            GATK_command += "-recalFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, racal_file_name_indels)
            GATK_command += "-tranchesFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, tranches_file_name_indels)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, racal_file_name_indels , working_dir)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, tranches_file_name_indels , working_dir)
        else:
            GATK_command += "-recalFile {}/VCF/{} \\\n".format(working_dir, racal_file_name_indels)
            GATK_command += "-tranchesFile {}/VCF/{} \\\n".format(working_dir, tranches_file_name_indels)
        VQSR.write(GATK_command)
        VQSR.write("\n")
        ############################################
        ##### now apply recalibration for INDELS ###
        ############################################
        GATK_command  = "java -Xmx64g -jar {} -T ApplyRecalibration  \\\n".format(CONFIG["GATK"])
        #### GATK_input is the same
        if scratch:
            GATK_command += "-recalFile $SNIC_TMP/{}/VCF/{} \\\n".format(job_name, racal_file_name_indels)
            GATK_command += "-tranchesFile $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, tranches_file_name_indels)
        else:
            GATK_command += "-recalFile {}/VCF/{} \\\n".format(working_dir, racal_file_name_indels)
            GATK_command += "-tranchesFile {}/VCF/{} \n\n".format(working_dir, tranches_file_name_indels)
        GATK_command += GATK_input
        #add standard options
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if isinstance(option, basestring):
                GATK_command += "{} \\\n".format(option)
        #now add specifc option for type
        added = False
        for option in CONFIG["walkers"]["ApplyRecalibration"]:
            if not isinstance(option, basestring) and "INDEL" in option:
                specific_options = option["INDEL"]
                added = True
                for specific_option in specific_options:
                    GATK_command += "{} \\\n".format(specific_option)
        if not added:
            print "WARNING: I did not inserted any specifc option in VQSR step, there should be either a SNP or an INDEL specific option"

        if scratch:
            GATK_command += "-o $SNIC_TMP/{}/VCF/{} \n\n".format(job_name, variant_recal_snp_recal_indels)
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/ \n".format(job_name, variant_recal_snp_recal_indels , working_dir)
        else:
            GATK_command += "-o {}/VCF/{} \n\n".format(working_dir, variant_recal_snp_recal_indels)

        VQSR.write(GATK_command)
        VQSR.write("\n")

    return sbatch_file

    



def VQSR():
    """Run VQSR
    
    :returns: the sbatch_file to be executed
    """
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "03_CatVariants")):
        sys.exit("Directory 03_CatVariants does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "04_VQSR")):
        print "WARNING: 04_VQSR already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "04_VQSR"))
        os.mkdir(os.path.join(cwd, "04_VQSR", "sbatch"))
        os.mkdir(os.path.join(cwd, "04_VQSR", "std_err"))
        os.mkdir(os.path.join(cwd, "04_VQSR", "std_out"))
        os.mkdir(os.path.join(cwd, "04_VQSR", "VCF"))
    #Build the sbatch files for merging step
    working_dir   = os.path.join(cwd, "04_VQSR")
    variant_file  = os.path.join(cwd, "03_CatVariants", "VCF", "{}_joincalled.g.vcf.gz".format(CONFIG["output_header"]))
    sbatch_file   = build_VQSR_sbatch(working_dir, variant_file, CONFIG["scratch"])
    sbatch_files.append(sbatch_file)
    return sbatch_files

