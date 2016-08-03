#!/bin/bash -l
#SBATCH -A ngi2016003
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 10-00:00:00
#SBATCH -J SelectVariants
#SBATCH -o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/std_out/SelectVariants.out
#SBATCH -e /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/std_err/SelectVariants.err
module load bioinfo-tools
module load GATK/3.5.0
java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants  \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/03_CatVariants/VCF/SRG_joincalled.g.vcf.gz \
-selectType SNP \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.snp.g.vcf.gz



java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants  \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-nt 16 \
-V /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/03_CatVariants/VCF/SRG_joincalled.g.vcf.gz \
-selectType INDEL \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.indel.g.vcf.gz

java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantEval -nt 16 \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-comp:hapmap /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/hapmap_3.3.b37.vcf \
-D  /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf \
--eval /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.snp.g.vcf.gz \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.snp.eval 


java -Xmx250g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantEval -nt 16 \
-R /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/human_g1k_v37.fasta \
-comp:hapmap /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/hapmap_3.3.b37.vcf \
-D  /lupus/ngi/resources/piper/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf \
--eval /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.indel.g.vcf.gz \
-o /lupus/proj/ngi2016003/nobackup/vezzi/SwedishReferenceGenomeProjects/develop/test/04_SelectVariants/VCF/SRG_joincalled.indel.eval 
