#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J ApplyRecalibration_INDEL
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/06_ApplyRecalibration/std_out/ApplyRecalibration_INDEL.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/06_ApplyRecalibration/std_err/ApplyRecalibration_INDEL.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T ApplyRecalibration  \
-input /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/04_SelectVariants/VCF/SweGen_joincalled.indel.g.vcf.gz \
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/05_VariantRecalibrator/VCF/SweGen_joincalled.indel.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/05_VariantRecalibrator/VCF/SweGen_joincalled.indel.tranches \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
--ts_filter_level 99.0 \
-mode INDEL \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_20/06_ApplyRecalibration/VCF/SweGen_joincalled.indel.recalibrated.filtered.vcf.gz 

