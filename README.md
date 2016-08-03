# joint_variant_calling
Python package to run Join Calling on population at NGI (National Genomics Infrastructure) Sweden.
The script `joint_variant_calling.py` implements the GATK-Workflow described in 
https://www.broadinstitute.org/gatk/guide/article?id=3893

Samples to be be join called can be specified in three ways:
 - in the sample field of the config.yaml file to be provided as input. Complete path to the gvcf file needs to be provided
 - in the project field of the config.yaml file to be provided as input. In this case NGI-specific folder structure is assumed.
 - a file named 00_samples.txt in the cwd. If present samples specified in this file will be used instead of those speccified in the config

if run like
 ``` python joint_variant_calling.py --config config.yaml ```
it creates the following folder structure 
 - `00_intervals`: optional see Intervals section
 - `00_samples.txt`: 
 - `01_CombineGVCFs`
 - `02_GenotypeGVCFs`
 - `03_CatVariants`
 - `04_SelectVariants`
 - `05_VariantRecalibrator`
 - `006_ApplyRecalibration`


## Intervals
