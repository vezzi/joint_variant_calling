import sys, os, glob
import argparse
from utils import config as conf

from common import submit_jobs, find_VCF
from walkers.CombineGVCF import  CombineGVCFs, CombineGVCFs_resume
from walkers.GenotypeGVCFs import  GenotypeGVCFs
from walkers.CatVariants import   CatVariants
from walkers.SelectVariants import  SelectVariants
from walkers.VariantRecalibrator import VariantRecalibrator
from walkers.ApplyRecalibration import ApplyRecalibration

from utils.config import CONFIG

"""
 this is the main file, it checks that input provided makes sense and tries to run one after the other
 the jobs. The wrokflow is fixed and cannot be changed but the parameters can be changed via input parametrs
 and/or via the yaml config file
"""

def check_configuration():
    """This function checks that the configuraiton file provided makes sense.
    
    :param dict configuration: configuration file that is going to be checked
    
    :returns: True if all checks succed, false otherwise
    """
    mandatory_args = ('dry_run', 'scratch', 'batch_size', 'output_header', 'intervals', 'walkers', 'samples', 'uppmax_projects', 'projects' , 'GATK')
    for mandatory_arg in mandatory_args:
        if mandatory_arg not in CONFIG:
            print "ERROR: argument {} is mandatory. If you think it does not apply leave it empty".format(mandatory_arg)
            return False
    #check tha all walkers are present
    walkers = ('CombineGVCFs','GenotypeGVCFs','CatVariants','SelectVariants','VariantEval','VariantRecalibrator', 'ApplyRecalibration')
    for walker in walkers:
        if walker not in CONFIG["walkers"]:
            print "ERROR: walker {} is mandatory. please specify it".format(walker)
            return False

    #now at least I know all keys are present
    if CONFIG["samples"] is None and CONFIG["projects"] is None:
        print "ERROR: at least one one between projects and/or  samples must contain a list"
        return False
    #if project is not none than uppmax-project need do be specifed
    if CONFIG["projects"] is not None and CONFIG["uppmax_projects"] is None:
        print "ERROR: if projects are specified their uppmax-project need to be specified as well"
        return False
    #check tha provided uppmax projects exists (if provided)
    if CONFIG["uppmax_projects"] is not None:
        for uppmax_project in CONFIG["uppmax_projects"]:
            if not os.path.exists("/proj/{}".format(uppmax_project)):
                print "ERROR: uppmax project {} does not exists.".format(uppmax_project)
                return False
    #do the same sanity check for samples and create the to be join called list
    CONFIG["samples_JC"] = []
    if CONFIG["samples"] is not None:
        for sample in CONFIG["samples"]:
            if not  os.path.exists("{}".format(sample)):
                print "ERROR: sample  {} does not exists.".format(sample)
                return False
            CONFIG["samples_JC"].append(sample)
    if CONFIG["projects"] is not None:
        for project in CONFIG["projects"]:
            for uppmax_project in CONFIG["uppmax_projects"]:
                #for each project check the directory that is supposed to contain the abalysis
                CONFIG["samples_JC"] += find_VCF(project, uppmax_project)
    #check for duplicates, if identified stop computation
    samples_to_be_joint_called_hash = {}
    for sample in CONFIG["samples_JC"]:
        if sample not in samples_to_be_joint_called_hash:
            samples_to_be_joint_called_hash[sample] = 0
        else:
            print "ERROR: sample {} found twice".format(sample)
            return False
    #create list of the intervals if intervals is specified (if it is not valid proceed anyway without intervals)
    CONFIG["intervals_list"] = []
    if os.path.exists(CONFIG["intervals"]):
        for interval in glob.glob("{}/*.intervals".format(CONFIG["intervals"])):
            CONFIG["intervals_list"].append(interval)
        if len(CONFIG["intervals_list"]) == 0:
            print "WARNING: --intervals specified, but no interval was detected. Proceeding with whole genome, check if the inteval folder contains intervals in the form of .intervals"
    #initialise batch number to 1
    CONFIG["batch_number"] = 1
    return True




def main(args):
    config = conf.load_yaml_config(args.configuration)
    if not check_configuration():
        sys.exit("ERROR: configuraiton file was malformed, please edit it and retry")
    #store in a file path to vcf that are going to be analysed
    if args.resume and os.path.exists("00_samples.txt"):
        sys.exit("ERROR: --resume  specified, however 00_samples.txt found. Please if you want to resume analysis, remove/move 00_samples.txt, 02_GenotypeGVCF, 03_ ... ")
    if os.path.exists("00_samples.txt"):
        print "WARNING: file 00_samples.txt exists, I will replace samples in config file with those stored in this file."
        #delete samples sotred in previous step
        CONFIG["samples_JC"] = []
        with open("00_samples.txt", "r") as samplesFile:
            for sample in samplesFile:
                CONFIG["samples_JC"].append(sample)
    elif not args.resume:
        #otherwise create the file and store the samples
        with open("00_samples.txt", "w") as samplesFile:
            for sample in CONFIG["samples_JC"]:
                samplesFile.write("{}\n".format(sample))
    ## IMPORTANT: samples_JC contains samples to be JointCalled.
    ######### START JOIN CALLING OF THE VARIANTS #####################
    ##### https://www.broadinstitute.org/gatk/guide/article?id=3893
    ##################################################################
    #now join batches of batch_size samples
    if args.resume:
        #recompute only last batch of sample and, in case the extra ones
        sbatch_files = CombineGVCFs_resume()
    else:
        #start from scratch
        sbatch_files = CombineGVCFs()
    slurm_jobs_id = None
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files)
    #now perform the GenotypeGVCF step
    sbatch_files = GenotypeGVCFs()
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files, slurm_jobs_id)
    #at this point merge the chr into a single one
    sbatch_files = CatVariants()
    #and execute
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files, slurm_jobs_id)
    #now perofmr VQSR
    #start with select variants and variant evaluation
    sbatch_files = SelectVariants()
    #and execute
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files, slurm_jobs_id)
    #then perfomr VQSR
    sbatch_files = VariantRecalibrator()
    #and execute
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files, slurm_jobs_id)
    #than ApplyRecalibration
    sbatch_files = ApplyRecalibration()
    if not CONFIG["dry_run"]:
        slurm_jobs_id = submit_jobs(sbatch_files, slurm_jobs_id)




if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Scripts performs join variant calling on all samples provided, or on all samples analysed belonging to the projects specified""")
    parser.add_argument('--configuration', help="configuration file, give a look to the example one to see what to do", type=str)
    parser.add_argument('--resume', action='store_true', help="if this option is specified 01_Combine_GVCFs needs to be created and populated, all other folders need to be deleted/moved. In this modality the script deletes only the last batch (if needed) and creates a n new (or more) extra batches to account for new samples",  default=False )
    args = parser.parse_args()
    main(args)



