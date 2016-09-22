#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J GenotypeGVCFs_interval0005
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/02_GenotypeGVCFs/std_out/GenotypeGVCFs_interval0005.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/02_GenotypeGVCFs/std_err/GenotypeGVCFs_interval0005.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T GenotypeGVCFs  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
--max_alternate_alleles 6 \
-D /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/01_CombineGVCFs/VCF/SweGen_batch1_interval0005.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/01_CombineGVCFs/VCF/SweGen_batch2_interval0005.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/01_CombineGVCFs/VCF/SweGen_batch3_interval0005.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/01_CombineGVCFs/VCF/SweGen_batch4_interval0005.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/01_CombineGVCFs/VCF/SweGen_batch5_interval0005.g.vcf.gz \
 -L /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/00_intervals/interval0005.intervals \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/02_GenotypeGVCFs/VCF/SweGen_joincalled_interval0005.g.vcf.gz

