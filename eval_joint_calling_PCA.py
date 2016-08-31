import sys, os, glob
import argparse
from utils import config as conf

from common import submit_jobs
from utils.config import CONFIG



def sbatch_header(time="1-00:00:00", uppmax_project="ngi2016003", sbatch_name="SGP", cwd=""):
    sbatch_header = ""
    sbatch_header += "#!/bin/bash -l\n"
    sbatch_header += "#SBATCH -A {}\n".format(uppmax_project)
    sbatch_header += "#SBATCH -p node\n"
    sbatch_header += "#SBATCH -n 16\n"
    sbatch_header += "#SBATCH -t {}\n".format(time)
    sbatch_header += "#SBATCH -J {}\n".format(sbatch_name)
    sbatch_header += "#SBATCH -o {}/{}.out\n".format(cwd,sbatch_name)
    sbatch_header += "#SBATCH -e {}/{}.err\n".format(cwd,sbatch_name)
    sbatch_header += "\n"
    sbatch_header += "module load bioinfo-tools\n"
    sbatch_header += "module load vcftools\n"
    sbatch_header += "module load GATK/3.5.0\n"
    sbatch_header += "\n"
    return sbatch_header


def merge_snps_and_indels():
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "01_merge_snp_indels")):
        print "WARNING: 01_merge_snp_indels already present, assuming this step has been completed with success."
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, "01_merge_snp_indels"))

    sbatch_command = sbatch_header(sbatch_name="merge_snp_indels", cwd=os.path.join(cwd, "01_merge_snp_indels"))
    sbatch_command += "java -Xmx250g -jar {} -T CombineVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "--variant:snps {} \\\n".format(CONFIG["popSNPs"])
    sbatch_command += "--variant:indels {} \\\n".format(CONFIG["popINDELs"])
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "01_merge_snp_indels", "SRG_joincalled.snp.indels.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "-genotypeMergeOptions PRIORITIZE \\\n"
    sbatch_command += "-priority snps,indels \\\n"
    sbatch_command += "\n"
    with open(os.path.join(cwd, "01_merge_snp_indels", "00_merge_snp_indels.sbatch"), "w") as MERGE:
        MERGE.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, "01_merge_snp_indels", "00_merge_snp_indels.sbatch")])
    return slurm_jobs_id

def select_1KGP():
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "02_select_variants_1KGP")):
        print "WARNING: 02_select_variants_1KGP already present, assuming this step has been completed with success."
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, "02_select_variants_1KGP"))

    sbatch_command = sbatch_header(sbatch_name="select_variants_1KGP", cwd=os.path.join(cwd, "02_select_variants_1KGP"))
    sbatch_command += "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-V {} \\\n".format(CONFIG["1KGP_VCF"])
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "02_select_variants_1KGP", "1KGP_selected.snp.indels.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "\n"
    with open(os.path.join(cwd, "02_select_variants_1KGP", "00_select_variants_1KGP.sbatch"), "w") as SELECT:
        SELECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, "02_select_variants_1KGP", "00_select_variants_1KGP.sbatch")])
    return slurm_jobs_id


def select_SGP(jobs_id):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "03_select_variants_SGP")):
        print "WARNING: 03_select_variants_SGP already present, assuming this step has been completed with success."
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, "03_select_variants_SGP"))

    sbatch_command = sbatch_header(sbatch_name="select_variants_SGP", cwd=os.path.join(cwd, "03_select_variants_SGP"))
    sbatch_command +=  "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-V {} \\\n".format(os.path.join(cwd, "01_merge_snp_indels", "SRG_joincalled.snp.indels.vcf.gz"))
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "03_select_variants_SGP", "SRG_joincalled.filtered.snp.indels.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "--restrictAllelesTo BIALLELIC \\\n"
    sbatch_command += "--excludeFiltered \\\n"

    sbatch_command += "\n"
    with open(os.path.join(cwd, "03_select_variants_SGP", "00_select_variants_SGP.sbatch"), "w") as SELECT:
        SELECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, "03_select_variants_SGP", "00_select_variants_SGP.sbatch")], jobs_id)
    return slurm_jobs_id



def merge_with_1KGP(jobs_id):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "04_merge_1KGP_SGP")):
        print "WARNING: 04_merge_1KGP_SGP already present, assuming this step has been completed with success."
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, "04_merge_1KGP_SGP"))
    
    sbatch_command = sbatch_header(sbatch_name="merge_1KGP_SGP", cwd=os.path.join(cwd, "04_merge_1KGP_SGP"))

    sbatch_command +=  "java -Xmx250g -jar {} -T CombineVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "--variant:SGP {} \\\n".format(os.path.join(cwd, "03_select_variants_SGP", "SRG_joincalled.filtered.snp.indels.vcf.gz"))
    sbatch_command += "--variant:1KGP {} \\\n".format(os.path.join(cwd, "02_select_variants_1KGP", "1KGP_selected.snp.indels.vcf.gz"))
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "04_merge_1KGP_SGP", "1KGP_SGP.vcf.gz")
    sbatch_command += "-o {} \n".format(output)

    sbatch_command += "\n"

    sbatch_command +=  "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    sbatch_command += "-V {} \\\n".format(os.path.join(cwd, "04_merge_1KGP_SGP", "1KGP_SGP.vcf.gz"))
    sbatch_command += "-select \'set == \"Intersection\"\' \\\n"
    output = os.path.join(cwd, "04_merge_1KGP_SGP", "1KGP_SGP.intersection.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "\n"
    
    with open(os.path.join(cwd, "04_merge_1KGP_SGP", "00_merge_1KGP_SGP.sbatch"), "w") as INTERSECT:
        INTERSECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, "04_merge_1KGP_SGP", "00_merge_1KGP_SGP.sbatch")], jobs_id)
    return slurm_jobs_id


def select_EU_samples():
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "08_EU_1KGP_SGP")):
        print "WARNING: 08_EU_1KGP_SGP already present, assuming this step has been completed with success."
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, "08_EU_1KGP_SGP"))



    sbatch_command = sbatch_header(sbatch_name="EU_1KGP_SGP", cwd=os.path.join(cwd, "08_EU_1KGP_SGP"))
    sbatch_command +=  "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-V {} \\\n".format(os.path.join(cwd, "04_merge_1KGP_SGP", "1KGP_SGP.intersection.vcf.gz"))
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "08_EU_1KGP_SGP", "EU_1KGP_SGP.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "-sf  {}\\\n".format(CONFIG["EU_samples"])

    sbatch_command += "\n"
    with open(os.path.join(cwd, "08_EU_1KGP_SGP", "00_select_EU_1KGP_SGP.sbatch"), "w") as SELECT:
        SELECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, "08_EU_1KGP_SGP", "00_select_EU_1KGP_SGP.sbatch")])
    return slurm_jobs_id



def runPCA(folder,output, VCF, populations):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, folder)):
        print "WARNING: {} already present, assuming this step has been completed with success.".format(folder)
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, folder))
    
    
    sbatch_command = sbatch_header(sbatch_name=folder, cwd=os.path.join(cwd, folder))
    #create tbed
    output_folder = os.path.join(cwd, folder)
    sbatch_command +=  "vcftools --gzvcf {} --plink-tped --out {}/{} \n".format(VCF, output_folder, output)
    #run plink on this set
    sbatch_command +=  "{} -tfile {}/{}  --pca --out {}/{}_PCA \n".format(CONFIG["PLINK"], output_folder, output, output_folder, output)
    #create PCA table with population
    sbatch_command += "python {} --pca {}/{}_PCA.eigenvec --populations ".format(CONFIG["PCA_to_plink"], output_folder, output)
    for population in populations:
        sbatch_command += " {} ".format(population)
    sbatch_command += " > {}/{}_PCA.pop.eigenvec \n".format(output_folder, output)
    sbatch_command += "\n"
    with open(os.path.join(cwd, folder, "00_runPCA.sbatch"), "w") as PCA:
        PCA.write(sbatch_command)
    slurm_jobs_id = submit_jobs([os.path.join(cwd, folder, "00_runPCA.sbatch")])



def main(args):
    config = conf.load_yaml_config(args.configuration)
    #merge SNPSs and INDELS
    jobs_id_merge = merge_snps_and_indels()
    #run QC on this set
    jobs_id_select_SGP = select_SGP(jobs_id_merge)
    #select same inteval for 1KGP
    jobs_id_select = select_1KGP()
    #merge 1KGP and SGP
    jobs_dependencies = []
    if jobs_id_select_SGP is not None:
        jobs_dependencies.append(jobs_id_select_SGP)
    if jobs_id_select is not None:
        jobs_dependencies.append(jobs_id_select)
    if len(jobs_dependencies) == 0:
        jobs_dependencies = None
    merge_with_1KGP(jobs_dependencies)
    #now run PCA, specify which vcf and which population to work on
    cwd = os.getcwd()
    SGPvcf = os.path.join(cwd,"03_select_variants_SGP", "SRG_joincalled.filtered.snp.indels.vcf.gz")
    runPCA("05_PCA_SGP_only", "SGP",  SGPvcf, [CONFIG["SGP_population"]])
    #now only for 1KGP
    OneKGP_vcf = os.path.join(cwd, "02_select_variants_1KGP", "1KGP_selected.snp.indels.vcf.gz")
    runPCA("06_PCA_1KGP_only", "1KGP",  OneKGP_vcf, [CONFIG["1KGP_superpopulation"]])
    #now for the mix one
    SGP_1KGP = os.path.join(cwd, "04_merge_1KGP_SGP", "1KGP_SGP.intersection.vcf.gz")
    runPCA("07_PCA_SGP_1KGP", "SGP_1KGP",  SGP_1KGP, [CONFIG["1KGP_superpopulation"], CONFIG["SGP_superpopulation"]])
    #select only EU samples from merges vcf file
    select_EU_samples()
    #run PCA on this set
    EU_SGP_1KGP = os.path.join(cwd, "08_EU_1KGP_SGP", "EU_1KGP_SGP.vcf.gz")
    runPCA("09_PCA_EU_SGP_1KGP", "EU_SGP_1KGP",  EU_SGP_1KGP, [CONFIG["1KGP_superpopulation"], CONFIG["SGP_superpopulation"]])

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Scripts performs standard validation steps on joint calling results""")
    parser.add_argument('--configuration', help="configuration file, give a look to the example one to see what to do", type=str)
    args = parser.parse_args()
    main(args)



