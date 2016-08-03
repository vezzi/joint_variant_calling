#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J VQSR_INDEL
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/05_VariantRecalibrator/std_out/VQSR_INDEL.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/05_VariantRecalibrator/std_err/VQSR_INDEL.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantRecalibrator  \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
--maxGaussians 4 \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf \
-an QD \
-an DP \
-an FS \
-an SOR \
-an ReadPosRankSum \
-an MQRankSum \
-mode INDEL \
-input /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.indel.g.vcf.gz \
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/05_VariantRecalibrator/VCF/SRG_joincalled.indel.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/05_VariantRecalibrator/VCF/SRG_joincalled.indel.tranches 

