import sys, os, glob
import argparse
import re
import subprocess



def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



def find_VCF(project, uppmax_project):
    """Finds all VCF files belonging to a project (NGI) analysed in uppmax_project (assumes NGI/irma organisation).
    
    :param str project: NGI-project (in case of Stcokholm P123)
    :param srt uppmax_project: the uppmax project name (e.g., ngi2016003)
    
    :returns: list of all vcf files matvhing the search pattern
    
    """
    samples= []
    #where data is supposed to be
    analysis_dir = "/proj/{}/nobackup/NGI/ANALYSIS/{}/piper_ngi/07_variant_calls/".format(uppmax_project, project)
    for sample in glob.glob("{}*clean.dedup.recal.bam.genomic.vcf.gz".format(analysis_dir)):
        samples.append(sample)
    return samples




def submit_jobs(sbatch_files, pending_jobs = None):
    """Submits to slurm queue the sbatch files contained in the list. Retuns jobs id.
    
    :param list sbatchfiles
    
    :retuns list: slurm ids
    """
    slurm_jobs_id = None
    for sbatch_file in sbatch_files:
        #submit the sbatch file
        command = ["sbatch", "{}".format(sbatch_file)]
        if pending_jobs is not None:
            dependency = "--dependency=afterany"
            for pending_job in pending_jobs:
                dependency+=":{}".format(pending_job)
            command = ["sbatch", dependency, "{}".format(sbatch_file)]
        p_handle = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        p_out, p_err = p_handle.communicate()
        try:
            slurm_job_id = re.match(r'Submitted batch job (\d+)', p_out).groups()[0]
            if slurm_jobs_id is None:
                slurm_jobs_id = []
            slurm_jobs_id.append(slurm_job_id)
        except AttributeError:
            raise RuntimeError('Could not submit sbatch job {}'
                           '{}'.format(sbatch_file, p_err))
    #return the list of submitted jobs
    return slurm_jobs_id


