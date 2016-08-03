#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J CatVariants
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/03_CatVariants/std_out/CatVariants.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/03_CatVariants/std_err/CatVariants.err
module load bioinfo-tools
module load GATK/3.5.0
java -cp /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-assumeSorted \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/VCF/SRG_joincalled_interval0001.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/VCF/SRG_joincalled_interval0002.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/VCF/SRG_joincalled_interval0003.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/VCF/SRG_joincalled_interval0004.g.vcf.gz \
 -out /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/03_CatVariants/VCF/SRG_joincalled.g.vcf.gz

