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


def merge_snps_and_indels(step, jobs_id=None):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, step)):
        print "WARNING: {} already present, assuming this step has been completed with success.".format(step)
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, step))

    sbatch_command = sbatch_header(sbatch_name=step, cwd=os.path.join(cwd, step))
    sbatch_command += "java -Xmx250g -jar {} -T CombineVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "--variant:snps {} \\\n".format(CONFIG["popSNPs"])
    sbatch_command += "--variant:indels {} \\\n".format(CONFIG["popINDELs"])
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, "01_merge_snp_indels", "SGP_joincalled.snp.indels.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "-genotypeMergeOptions PRIORITIZE \\\n"
    sbatch_command += "-priority snps,indels \\\n"
    sbatch_command += "\n"
    with open(os.path.join(cwd, step, "{}.sbatch".format(step)), "w") as MERGE:
        MERGE.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, step, "{}.sbatch".format(step))], jobs_id)
    return slurm_jobs_id


def select(step , vcf_in, vcf_out, options, jobs_id=None):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, step)):
        print "WARNING: {} already present, assuming this step has been completed with success.".format(step)
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, step))

    sbatch_command = sbatch_header(sbatch_name=step, cwd=os.path.join(cwd, step))
    sbatch_command += "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-V {} \\\n".format(vcf_in)
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, step, vcf_out)
    sbatch_command += "-o {} \\\n".format(output)
    for option in options:
        sbatch_command += "{} \\\n".format(option)
    sbatch_command += "\n"
    with open(os.path.join(cwd, step, "{}.sbatch".format(step)), "w") as SELECT:
        SELECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, step, "{}.sbatch".format(step))], jobs_id)
    return slurm_jobs_id


def merge_with_1KGP(step, vcf_one, vcf_two, jobs_id):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, step)):
        print "WARNING: {} already present, assuming this step has been completed with success.".format(step)
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, step))
    
    sbatch_command = sbatch_header(sbatch_name=step, cwd=os.path.join(cwd, step))

    sbatch_command +=  "java -Xmx250g -jar {} -T CombineVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "--variant:SGP {} \\\n".format(vcf_one)
    sbatch_command += "--variant:1KGP {} \\\n".format(vcf_two)
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, step, "1KGP_SGP.vcf.gz")
    sbatch_command += "-o {} \n".format(output)

    sbatch_command += "\n"

    sbatch_command +=  "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    sbatch_command += "-V {} \\\n".format(os.path.join(cwd, step, "1KGP_SGP.vcf.gz"))
    sbatch_command += "-select \'set == \"Intersection\"\' \\\n"
    output = os.path.join(cwd, step, "1KGP_SGP.intersection.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "\n"
    
    with open(os.path.join(cwd, step, "{}.sbatch".format(step)), "w") as INTERSECT:
        INTERSECT.write(sbatch_command)

    slurm_jobs_id = submit_jobs([os.path.join(cwd, step, "{}.sbatch".format(step))], jobs_id)
    return slurm_jobs_id


def select_EU_samples(step, vcf_all, jobs_id):
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, step)):
        print "WARNING: {} already present, assuming this step has been completed with success.".format(step)
        return
    #create the folder structure
    os.mkdir(os.path.join(cwd, step))

    sbatch_command = sbatch_header(sbatch_name=step, cwd=os.path.join(cwd, step))
    sbatch_command +=  "java -Xmx250g -jar {} -T SelectVariants  \\\n".format(CONFIG["GATK"])
    sbatch_command += "-R {} \\\n".format(CONFIG["reference"])
    sbatch_command += "-nt 16  \\\n"
    sbatch_command += "-V {} \\\n".format(vcf_all)
    sbatch_command += "-L {} \\\n".format(CONFIG["intervals"])
    output = os.path.join(cwd, step, "EU_1KGP_SGP.vcf.gz")
    sbatch_command += "-o {} \\\n".format(output)
    sbatch_command += "-sf  {}\\\n".format(CONFIG["EU_samples"])

    sbatch_command += "\n"
    with open(os.path.join(cwd, step, "{}.sbatch".format(step)), "w") as SELECT:
        SELECT.write(sbatch_command)
    slurm_jobs_id = submit_jobs([os.path.join(cwd, step, "{}.sbatch".format(step))], jobs_id)
    return slurm_jobs_id



def runPCA(folder, output, VCF, populations, jobs_id=None):
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
    slurm_jobs_id = submit_jobs([os.path.join(cwd, folder, "00_runPCA.sbatch")], jobs_id)




def main(args):
    config = conf.load_yaml_config(args.configuration)
    #01 merge SNPSs and INDELS for SweGene
    jobs_id_merge = merge_snps_and_indels(step="01_merge_snp_indels")
    SGP_raw_vcf   = os.path.join(os.getcwd(), "01_merge_snp_indels", "SGP_joincalled.snp.indels.vcf.gz")
    #02 select subset of vcf file merged at step 01
    jobs_id_select_SGP = select( step    = "02_select_variants_SGP" ,
                                 vcf_in  = SGP_raw_vcf ,
                                 vcf_out = "SGP_joincalled.filtered.snp.indels.vcf.gz" ,
                                 options = ["--restrictAllelesTo BIALLELIC", "--excludeFiltered"],
                                 jobs_id = jobs_id_merge
                                 )
    SGP_vcf = os.path.join(os.getcwd(), "02_select_variants_SGP", "SGP_joincalled.filtered.snp.indels.vcf.gz")
    #03 select subset of 1KGP
    jobs_id_select_1KGP = select(step    = "03_select_variants_1KGP" ,
                                 vcf_in  =  CONFIG["1KGP_VCF"],
                                 vcf_out = "1KGP_selected.snp.indels.vcf.gz" ,
                                 options = [],
                                 )
    OneKGP_vcf = os.path.join(os.getcwd(), "03_select_variants_1KGP", "1KGP_selected.snp.indels.vcf.gz")
    #merge 1KGP and SGP
    jobs_dependencies = []
    if jobs_id_select_SGP is not None:
        jobs_dependencies.extend(jobs_id_select_SGP)
    if jobs_id_select_1KGP is not None:
        jobs_dependencies.extend(jobs_id_select_1KGP)
    if len(jobs_dependencies) == 0:
        jobs_dependencies = None

    jobs_id_merge_1KGP_SGP = merge_with_1KGP(step="04_merge_1KGP_SGP", vcf_one=SGP_vcf, vcf_two=OneKGP_vcf, jobs_id=jobs_dependencies)
    SGP_1KGP_vcf = os.path.join(os.getcwd(), "04_merge_1KGP_SGP", "1KGP_SGP.intersection.vcf.gz")

    #select only EU samples from merged vcf file
    jobs_id_select_EU = select_EU_samples(step="05_EU_1KGP_SGP", vcf_all=SGP_1KGP_vcf, jobs_id=jobs_id_merge_1KGP_SGP)
    EU_SGP_1KGP_vcf   = os.path.join(os.getcwd(), "05_EU_1KGP_SGP", "EU_1KGP_SGP.vcf.gz")

    #now run PCA, specify: FOLDER, OUTPUT_NAME, VCF, POPLUATION, DEPENDECY_IDs
    cwd = os.getcwd()
    runPCA("06_PCA_SGP_only", "SGP",  SGP_vcf, [CONFIG["SGP_population"]], jobs_id_select_SGP)

    #now only for 1KGP
    runPCA("07_PCA_1KGP_only", "1KGP",  OneKGP_vcf, [CONFIG["1KGP_superpopulation"]], jobs_id_select_1KGP)

    #now for the mix one
    runPCA("08_PCA_SGP_1KGP", "SGP_1KGP",  SGP_1KGP_vcf, [CONFIG["1KGP_superpopulation"], CONFIG["SGP_superpopulation"]], jobs_id_merge_1KGP_SGP)

    #run PCA EU population only
    runPCA("09_PCA_EU_SGP_1KGP", "EU_SGP_1KGP",  EU_SGP_1KGP_vcf, [CONFIG["1KGP_population"], CONFIG["SGP_superpopulation"]], jobs_id_select_EU)



#eigenvec_table <- read.table('SGP_1KGP_PCA.pop.eigenvec')
#plot(eigenvec_table[4:5], col=factor(eigenvec_table[,3]), main="PCA", xlab="first component", ylab="second component")
#legend("bottomright", legend=levels(factor(eigenvec_table[,3])), text.col=seq_along(levels(factor(eigenvec_table[,3]))), pch=19, col=seq_along(levels(factor(eigenvec_table[,3]))) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Scripts performs standard validation steps on joint calling results""")
    parser.add_argument('--configuration', help="configuration file, give a look to the example one to see what to do", type=str)
    args = parser.parse_args()
    main(args)



