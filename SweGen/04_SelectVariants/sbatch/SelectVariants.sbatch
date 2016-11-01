#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J SelectVariants
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/std_out/SelectVariants.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/std_err/SelectVariants.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/03_CatVariants/VCF/SweGen_joincalled.g.vcf.gz \
-selectType SNP \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.snp.g.vcf.gz



java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/03_CatVariants/VCF/SweGen_joincalled.g.vcf.gz \
-selectType INDEL \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.indel.g.vcf.gz

java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantEval -nt 16 \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-comp:hapmap /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/hapmap_3.3.b37.vcf \
-D  /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf \
--eval /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.snp.g.vcf.gz \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.snp.eval 


java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantEval -nt 16 \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-comp:hapmap /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/hapmap_3.3.b37.vcf \
-D  /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf \
--eval /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.indel.g.vcf.gz \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_SelectVariants/VCF/SweGen_joincalled.indel.eval 
