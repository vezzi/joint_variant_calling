import sys, os, glob
import argparse
import config as conf


from config import CONFIG


def main(args):
    config = conf.load_yaml_config(args.config)
    if "GATK" not in config:
        sys.exit("ERROR: GATK must be present in config file, pointing to GATK path")
    sbatch_file = build_VCF2PED_sbatch(args.vcf, args.fam)
    submit_jobs([sbatch_file])


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


def build_VCF2PED_sbatch(vcf,  fam):
    working_dir = os.getcwd()
    job_name  = "VCF2PED"
    vcf_file  = os.path.basename(vcf)
    output_file_header = vcf_file.split(".vcf.gz")[0]
    
    os.mkdir(os.path.join(working_dir, "sbatch"))
    os.mkdir(os.path.join(working_dir, "std_err"))
    os.mkdir(os.path.join(working_dir, "std_out"))
    os.mkdir(os.path.join(working_dir, "PED"))

    sbatch_file = os.path.join(working_dir, "sbatch", "{}.sbatch".format(job_name))
    with open(sbatch_file, "w") as VCF2PED:
        VCF2PED.write("#!/bin/bash -l\n")
        VCF2PED.write("#SBATCH -A ngi2016003\n")
        VCF2PED.write("#SBATCH -p node\n")
        VCF2PED.write("#SBATCH -n 16\n")
        VCF2PED.write("#SBATCH -t 10-00:00:00\n")
        VCF2PED.write("#SBATCH -J {}\n".format(job_name))
        VCF2PED.write("#SBATCH -o {}/std_out/{}.out\n".format(working_dir, job_name ))
        VCF2PED.write("#SBATCH -e {}/std_err/{}.err\n".format(working_dir, job_name ))
        VCF2PED.write("module load bioinfo-tools\n")
        VCF2PED.write("module load GATK/3.5.0\n")

        GATK_command  = "java -Xmx64g -jar {} -T VariantsToBinaryPed  \\\n".format(CONFIG["GATK"])
        for option in CONFIG["walkers"]["VariantsToBinaryPed"]:
            GATK_command += "{} \\\n".format(option)
        GATK_command += "-V {} \\\n".format(vcf)
        GATK_command += "-m {} \\\n".format(fam)
        output_file = os.path.join(working_dir, "PED", "{}".format(output_file_header))
        GATK_command += "-bed {} \\\n".format(output_file)
        GATK_command += "-bim {} \\\n".format(output_file)
        GATK_command += "-fam {} \n\n".format(output_file)
        
        VCF2PED.write("{}".format(GATK_command))

    return sbatch_file
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Utility script to convert a VCF into a binary PED file to be used for analysis in PLINK""")
    parser.add_argument('--config', help="Configuration file (with GATK path)", type=str, required=True)
    parser.add_argument('--vcf', help="path to vcf file, in vcf.gz format)", type=str, required=True)
    parser.add_argument('--fam', help="fam file, produce it in this way 'for SAMPLE in `bcftools query -l file.vcf.gz `; do echo $SAMPLE $SAMPLE 0 0 0 0 ; done' ")

    args = parser.parse_args()
    main(args)

