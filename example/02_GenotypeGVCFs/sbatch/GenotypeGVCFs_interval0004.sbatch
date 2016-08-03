#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J GenotypeGVCFs_interval0004
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/std_out/GenotypeGVCFs_interval0004.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/std_err/GenotypeGVCFs_interval0004.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T GenotypeGVCFs  \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
--max_alternate_alleles 3 \
-D /lupus/ngi/resources//piper/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/01_CombineGVCFs/VCF/SRG_batch1_interval0004.g.vcf.gz \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/01_CombineGVCFs/VCF/SRG_batch2_interval0004.g.vcf.gz \
 -L /proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/00_intervals/interval0004.intervals \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/02_GenotypeGVCFs/VCF/SRG_joincalled_interval0004.g.vcf.gz

