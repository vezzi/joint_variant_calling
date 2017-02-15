import sys, os, glob
import re

from utils.config import CONFIG
from common import slurm_header


def build_GenotypeGVCFs_sbatch(working_dir, combined_gvcf_files, scratch=False, interval=None):
    """Builds the sbatch file in order to combine genomics.vcf samples contained in current_batch in a single one.
    
    :param str working_dir: directory where files will be created
    :param int batch: batch number, and incremental number specifing which batch lot are we processing
    :param list current_batch: list containing the samples to be combined
    :param bool scratch: if True works on scratch
    :param string interval: if not none specifies a file containing the interval(s) to be combined
    
    :returns: path to the sbatch file 
    
    """
    
    name_batch1   = os.path.basename([item for item in combined_gvcf_files if "batch1" in item][0])
    interval_name = ""
    #there must be at least one batch so look for it, not elegant but works
    if name_batch1.split("batch1") != ".g.vcf.gz":
        interval_name = name_batch1.split("batch1")[1].split(".")[0]
    job_name      = "GenotypeGVCFs{}".format(interval_name)
    output_file   = "{}_joincalled{}.g.vcf.gz".format(CONFIG["output_header"], interval_name)
    #create the sbatch file to analyse the current batch of samples
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as GenotypeGVCFs:
        slurm = slurm_header(CONFIG["uppmax_project"],  job_name, working_dir)
        GenotypeGVCFs.write(slurm)
        GenotypeGVCFs.write("\n")
        #rsync to scratch all samples
        if scratch:
            GenotypeGVCFs.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            GenotypeGVCFs.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        #now cycle over the samples, build the GATK command
        combined_gvcf_string_input = ""
        for combined_gvcf in combined_gvcf_files:
            combined_gvcf_path_dir = combined_gvcf
            if scratch:
                GenotypeGVCFs.write("rsync -rptoDLv {}* $SNIC_TMP/{}/\n".format(combined_gvcf, job_name))
                combined_gvcf_name = os.path.basename(combined_gvcf)
                combined_gvcf_path_dir = "$SNIC_TMP/{}/{}".format(job_name, combined_gvcf_name)
            combined_gvcf_string_input += "-V {} \\\n".format(combined_gvcf_path_dir)

        GATK_command= "java -Xmx250g -jar {} -T GenotypeGVCFs  \\\n".format(CONFIG["GATK"])
        for option in CONFIG["walkers"]["GenotypeGVCFs"]:
            GATK_command += "{} \\\n".format(option)
        GATK_command += "{} ".format(combined_gvcf_string_input)
        if interval is not None:
            GATK_command += "-L {} \\\n".format(interval)

        if scratch:
            GATK_command +=  "-o $SNIC_TMP/{}/VCF/{}\n".format(job_name, output_file)
            #once this is done rsync back to lupus
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/\n".format(job_name, output_file , working_dir)
        else:
            GATK_command += "-o {}/VCF/{}\n\n".format(working_dir, output_file)
        GenotypeGVCFs.write(GATK_command)
    #return path to sbach file
    return sbatch_file




def GenotypeGVCFs():
    """Runs  GenotypeGVCFs on all combined files produced previosuly (assumes folder structure)
    
    :param list intervals: intervals to be used in case we want to split analysis chr by chr
    :param int batch_numebr: numbers of batchs
    :param bool scratch: if True works on scratch
    
    :returns: list sbatch_files: list containing the sbatch files to be started
    
    """
    #creates sbatch files to merge batches of batch_size genomics vcf
    cwd = os.getcwd()
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "01_CombineGVCFs")):
        sys.exit("Directory 01_CombineGVCFs does not exits exists, something went wrong here.")
    if os.path.isdir(os.path.join(cwd, "02_GenotypeGVCFs")):
        print "WARNING: 02_GenotypeGVCFs already present, assuming this step has been completed with success."
        return sbatch_files
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "02_GenotypeGVCFs"))
        os.mkdir(os.path.join(cwd, "02_GenotypeGVCFs", "sbatch"))
        os.mkdir(os.path.join(cwd, "02_GenotypeGVCFs", "std_err"))
        os.mkdir(os.path.join(cwd, "02_GenotypeGVCFs", "std_out"))
        os.mkdir(os.path.join(cwd, "02_GenotypeGVCFs", "VCF"))
    #Build the sbatch files for the join calling step
    working_dir = os.path.join(cwd, "02_GenotypeGVCFs")
    #now retrive the VCF stored in 01_CombineGVCFs/VCF/
    combined_gvcfs_to_process = []
    if len(CONFIG["intervals_list"]) == 0:
        #no intervals, I have one file for each batch
        combined_gvcf_files = []
        for current_batch in range(1, CONFIG["batch_number"] +1):
            # for each batch create the vcf file that need to be created by combine step
            combined_gvcf_name      = "{}_batch{}.g.vcf.gz".format(CONFIG["output_header"], current_batch)
            combined_gvcf_full_path = os.path.join(cwd, "01_CombineGVCFs", "VCF", combined_gvcf_name)
            combined_gvcf_files.append(combined_gvcf_full_path)
        combined_gvcfs_to_process.append(combined_gvcf_files)
    else:
        for interval in CONFIG["intervals_list"]:
            interval_name = os.path.basename(interval).split(".")[0]
            combined_gvcf_files = []
            for current_batch in range(1, CONFIG["batch_number"] +1):
                # for each batch create the vcf file that need to be created by combine step
                combined_gvcf_name  = "{}_batch{}_{}.g.vcf.gz".format(CONFIG["output_header"], current_batch, interval_name)
                combined_gvcf_full_path = os.path.join(cwd, "01_CombineGVCFs", "VCF", combined_gvcf_name)
                combined_gvcf_files.append(combined_gvcf_full_path)
            #now ceate a list with interval file and all gvcf to be combines
            interval_plus_gvcfs = [interval ,combined_gvcf_files]
            combined_gvcfs_to_process.append(interval_plus_gvcfs)
    for interval_plus_gvcfs in combined_gvcfs_to_process:
        interval = interval_plus_gvcfs[0]
        combined_gvcf_files = interval_plus_gvcfs[1]
        sbatch_file = build_GenotypeGVCFs_sbatch(working_dir, combined_gvcf_files, CONFIG["scratch"], interval)
        sbatch_files.append(sbatch_file)
    return sbatch_files



