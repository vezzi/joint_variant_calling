import sys, os, glob
import random
import re

from utils.config import CONFIG


def build_CombineGVCFs_sbatch(working_dir, batch, current_batch, scratch=False, interval=None):
    """Builds the sbatch file in order to combine genomics.vcf samples contained in current_batch in a single one.
    
    :param str working_dir: directory where files will be created
    :param int batch: batch number, and incremental number specifing which batch lot are we processing
    :param list current_batch: list containing the samples to be combined
    :param bool scratch: if True works on scratch
    :param string interval: if not none specifies a file containing the interval(s) to be combined
    
    :returns: path to the sbatch file 
    
    """
    
    job_name      = "CombineGVCFs_batch{}".format(batch)
    output_file   = "{}_batch{}.g.vcf.gz".format(CONFIG["output_header"], batch)
    interval_name = ""
    if interval is not None:
            interval_name = os.path.basename(interval).split(".")[0] # store the interval name
            job_name    = "CombineGVCFs_batch{}_{}".format(batch, interval_name)
            output_file = "{}_batch{}_{}.g.vcf.gz".format(CONFIG["output_header"], batch, interval_name)
    #create the sbatch file to analyse the current batch of samples
    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as CombineGVCFsFile:
        CombineGVCFsFile.write("#!/bin/bash -l\n")
        CombineGVCFsFile.write("#SBATCH -A ngi2016003\n")
        CombineGVCFsFile.write("#SBATCH -p node\n")
        CombineGVCFsFile.write("#SBATCH -n 8\n")
        CombineGVCFsFile.write("#SBATCH -t 10-00:00:00\n")
        CombineGVCFsFile.write("#SBATCH -J {}\n".format(job_name))
        CombineGVCFsFile.write("#SBATCH -o {}/std_out/{}.out\n".format(working_dir, job_name ))
        CombineGVCFsFile.write("#SBATCH -e {}/std_err/{}.err\n".format(working_dir, job_name ))
        CombineGVCFsFile.write("module load bioinfo-tools\n")
        CombineGVCFsFile.write("module load GATK/3.5.0\n")
        #rsync to scratch all samples
        if scratch:
            CombineGVCFsFile.write("mkdir -p $SNIC_TMP/{} \n".format(job_name)) # create tmp directory
            CombineGVCFsFile.write("mkdir -p $SNIC_TMP/{}/VCF/ \n".format(job_name)) # create tmp directory
        #now cycle over the samples, build the GATK command
        samples_string_input = ""
        for sample in current_batch:
            sample_path_dir = sample
            if scratch:
                CombineGVCFsFile.write("rsync -rptoDLv {} $SNIC_TMP/{}/\n".format(sample, job_name))
                CombineGVCFsFile.write("rsync -rptoDLv {}.tbi $SNIC_TMP/{}/\n".format(sample, job_name))
                sample_name = os.path.basename(sample)
                sample_path_dir = "$SNIC_TMP/{}/{}".format(job_name, sample_name)
            samples_string_input += "-V {} \\\n".format(sample_path_dir)
        GATK_command= "java -Xmx200g -jar {} -T CombineGVCFs \\\n".format(CONFIG["GATK"])
        for option in CONFIG["walkers"]["CombineGVCFs"]:
            GATK_command += "{} \\\n".format(option)
        #attach the samples I am going to work with
        GATK_command += "{} ".format(samples_string_input)
        if interval is not None:
            GATK_command += "-L {} \\\n".format(interval)
        if scratch:
            GATK_command += "-o $SNIC_TMP/{}/VCF/{}\n".format(job_name, output_file)
            #once this is done rsync back to lupus
            GATK_command += "rsync $SNIC_TMP/{}/VCF/{}* {}/VCF/\n".format(job_name, output_file , working_dir)
        else:
            GATK_command += "-o {}/VCF/{}\n\n".format(working_dir, output_file)
        CombineGVCFsFile.write(GATK_command)
    #return path to sbach file
    return sbatch_file



def CombineGVCFs_resume():
    """
    Checks previously created batches, adds new samples to the last one and might create more batches
    
    :returns: list sbatch_files: list containing the number of batches and sbatch files to be started
    
    """
    #start with some sanity checks
    cwd = os.getcwd()
    working_dir = os.path.join(cwd, "01_CombineGVCFs")
    sbatch_files = []
    if not os.path.isdir(os.path.join(cwd, "01_CombineGVCFs")):
        sys.exit("ERROR: 01_CombineGVCFs must be present to resume JointCalling.")
    #check if other folders are present
    FOLDERS = ['02_GenotypeGVCFs', '03_CatVariants','04_SelectVariants','05_VariantRecalibrator','06_ApplyRecalibration',]
    for folder in FOLDERS:
        if os.path.isdir(os.path.join(cwd, folder)):
            sys.exit("ERROR: folder {} must be moved/remoed to resume JointCalling.".format(folder))
    if not os.path.isdir(os.path.join(cwd, "01_CombineGVCFs", "batches")):
        sys.exit("ERROR: 01_CombineGVCFs/batches must be present to resume JointCalling.")
    #now collect information about batches and their samples. Assume that they are called batch_00NM.txt
    batches_directory = os.path.join(cwd, "01_CombineGVCFs", "batches")
    batches = {}
    for batch in glob.glob("{}/batch*.txt".format(batches_directory)):
        CONFIG["batch_number"] += 1 # add a new batch
        batch_number = int(os.path.basename(batch).split("_")[1].split(".")[0]) # this is a batch id
        batches[batch_number] = []
        with open(batch, "r") as batch_file:
            for sample in batch_file:
                sample = sample.rstrip()
                batches[batch_number].append(sample)
    #now I know how many batches I have
    new_samples = []
    last_batch = CONFIG["batch_number"] - 1
    samples_in_last_batch = batches[last_batch]

    for sample in CONFIG["samples_JC"]:
        sample_found = False
        for batch in batches:
            if sample in batches[batch]:
                sample_found = True
        if not sample_found: #then this is a new sample
            new_samples.append(sample)
    #I need to create 00_samples.txt from here summing batches and new_samples
    with open("00_samples.txt", "w") as samplesFile:
        for batch in batches:
            for sample in batches[batch]:
                samplesFile.write("{}\n".format(sample))
        for sample in new_samples:
            samplesFile.write("{}\n".format(sample))
    #now I have new samples and I have the last batch
    if len(samples_in_last_batch) != CONFIG["batch_size"]:
        #in this case I need to delete the last batch and redo it
        CONFIG["batch_number"] -= 1
        new_samples += samples_in_last_batch #add last batch
        batch_file_to_delete = os.path.join(batches_directory, "batch_{0:04d}.txt".format(last_batch))
        os.remove(batch_file_to_delete)
        for combined_file_to_delete in glob.glob("{}/VCF/*batch{}_*".format(working_dir, last_batch)):
            os.remove(combined_file_to_delete)
    #Otherwise I have only to create new batches
    return create_batches(new_samples)



def CombineGVCFs():
    """
    Divides samples to be joint called in batches of batch_size samples. For each batch the CombineGVCF step is run.
    all inputs are taken from CONFIG. Stores in a file the batches in order to resume this step if more samples are added.
    
    :returns: list sbatch_files: list containing the number of batches and sbatch files to be started
    
    """
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "01_CombineGVCFs")):
        print "WARNING: 01_CombineGVCFs already present, assuming this step has been completed with success."
        #I need to set batch number
        number_of_processed_samples = 0
        for sample in CONFIG["samples_JC"]:
            number_of_processed_samples += 1
            if number_of_processed_samples % CONFIG["batch_size"] == 0:
                CONFIG["batch_number"] += 1
        if number_of_processed_samples % CONFIG["batch_size"] == 0:
            CONFIG["batch_number"] -= 1
        return []
    else:
        #create the folder structure
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs"))
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs", "sbatch"))
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs", "std_err"))
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs", "std_out"))
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs", "VCF"))
        os.mkdir(os.path.join(cwd, "01_CombineGVCFs", "batches")) #stores batches
    #Build the sbatch files for merging step
    random.shuffle(CONFIG["samples_JC"])
    return create_batches(CONFIG["samples_JC"])
    
    
def create_batches(samples):
    """
    takes as input the samples that need to be batched and creates the jobs
    """
    cwd = os.getcwd()
    sbatch_files = []
    working_dir = os.path.join(cwd, "01_CombineGVCFs")
    number_of_processed_samples = 0
    current_batch = [] # this will contain samples to be merged in current batch
    #shuffle the samples

    for sample in samples:
        current_batch.append(sample)
        number_of_processed_samples += 1
        if number_of_processed_samples % CONFIG["batch_size"] == 0:
            #I can build the sbatch file now, I need to check if I am working with genomics intervals or not
            if len(CONFIG["intervals_list"]) == 0:
                sbatch_file = build_CombineGVCFs_sbatch(working_dir, CONFIG["batch_number"], current_batch, CONFIG["scratch"], interval=None)
                sbatch_files.append(sbatch_file)
            else:
                for interval in CONFIG["intervals_list"]:
                    sbatch_file = build_CombineGVCFs_sbatch(working_dir, CONFIG["batch_number"], current_batch, CONFIG["scratch"], interval)
                    sbatch_files.append(sbatch_file)
            #Store the samples forming this batch
            sbatch_file = os.path.join(cwd, "01_CombineGVCFs", "batches", "batch_{0:04d}.txt".format(CONFIG["batch_number"]))
            with open(sbatch_file, "w") as batch_file:
                for element in current_batch:
                    batch_file.write("{}\n".format(element))
            CONFIG["batch_number"] += 1
            current_batch = []
    if len(current_batch) != 0:
        if len(CONFIG["intervals_list"]) == 0:
            sbatch_file = build_CombineGVCFs_sbatch(working_dir, CONFIG["batch_number"], current_batch, CONFIG["scratch"], interval=None)
            sbatch_files.append(sbatch_file)
        else:
            for interval in CONFIG["intervals_list"]:
                sbatch_file = build_CombineGVCFs_sbatch(working_dir, CONFIG["batch_number"], current_batch, CONFIG["scratch"], interval)
                sbatch_files.append(sbatch_file)
        sbatch_file = os.path.join(cwd, "01_CombineGVCFs", "batches", "batch_{0:04d}.txt".format(CONFIG["batch_number"]))
        with open(sbatch_file, "w") as batch_file:
            for element in current_batch:
                batch_file.write("{}\n".format(element))
    else:
        #decrement by one, special case if last batch has exactly batch_size samples
        CONFIG["batch_number"] -= 1

    return sbatch_files

