#!/bin/bash -l
#SBATCH -A a2014205
#SBATCH --qos seqver
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J VQSR
#SBATCH -o VQSR.out
#SBATCH -e VQSR.err
module load bioinfo-tools
module load GATK/3.5.0

#### following this
#### http://gatkforums.broadinstitute.org/gatk/discussion/2805/howto-recalibrate-variant-quality-scores-run-vqsr


#### Build the SNP recalibration model
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantRecalibrator  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/hapmap_3.3.b37.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/1000G_omni2.5.b37.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an DP \
-an InbreedingCoeff \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-input /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/03_CatVariants/VCF/SweGen_joincalled.g.vcf.gz
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.snp.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.snp.tranches 


####  Apply the desired level of recalibration to the SNPs in the call set
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T ApplyRecalibration  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-input SweGen_joincalled.g.vcf.gz \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.snp.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.snp.tranches \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.recal_snp_raw_indels.vcf.gz


#### Prepare recalibration parameters for Indels
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantRecalibrator  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
--maxGaussians 4 \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf \
-an QD \
-an DP \
-an FS \
-an SOR \
-an ReadPosRankSum \
-an MQRankSum \
-an InbreedingCoeff \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-input /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.recal_snp_raw_indels.vcf.gz \
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.indel.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.indel.tranches


###  Apply the desired level of recalibration to the Indels in the call set
java -Xmx64g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T ApplyRecalibration  \
-R /sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-input SweGen_joincalled.recal_snp_raw_indels.vcf.gz \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.indel.recal \
-tranchesFile /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.indel.tranches \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/SwedishGenomeReference_16_09_26/04_VQSR/SweGen_joincalled.recal_snp_recal_indels.vcf.gz
















