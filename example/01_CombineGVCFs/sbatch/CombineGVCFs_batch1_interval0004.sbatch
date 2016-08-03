#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 10-00:00:00
#SBATCH -J CombineGVCFs_batch1_interval0004
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/01_CombineGVCFs/std_out/CombineGVCFs_batch1_interval0004.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/01_CombineGVCFs/std_err/CombineGVCFs_batch1_interval0004.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx200g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T CombineGVCFs \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-V /proj/ngi2016003/nobackup/NGI/ANALYSIS/P2655/piper_ngi/07_variant_calls/P2655_102.clean.dedup.recal.bam.genomic.vcf.gz \
-V /proj/ngi2016003/nobackup/NGI/ANALYSIS/P2652/piper_ngi/07_variant_calls/P2652_213.clean.dedup.recal.bam.genomic.vcf.gz \
-V /proj/ngi2016003/nobackup/NGI/ANALYSIS/P2652/piper_ngi/07_variant_calls/P2652_241.clean.dedup.recal.bam.genomic.vcf.gz \
-V /proj/ngi2016003/nobackup/NGI/ANALYSIS/P2652/piper_ngi/07_variant_calls/P2652_160.clean.dedup.recal.bam.genomic.vcf.gz \
 -L /proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/00_intervals/interval0004.intervals \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/01_CombineGVCFs/VCF/SRG_batch1_interval0004.g.vcf.gz

