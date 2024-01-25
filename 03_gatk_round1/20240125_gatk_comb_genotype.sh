#!/bin/bash


##----------------------------------------------------------------------------------------
## LOAD MODULES

module load samtools/1.15.1
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.12.1 
module load picard-tools/2.17.10



#----------------------------------------------------------------------------------------

#################### CONFIGURATION INFORMATION ####################

# Location of the file containing the list of sample names
file="/home/awijeratne/class/20240127_NRT/subsetAra/data/sample.name.txt"

# Main project directory
PROJECT_DIR=/home/awijeratne/class/20240127_NRT/subsetAra

# Directory containing raw data
RAW_DATA=/home/awijeratne/class/20240127_NRT/subsetAra/data

# Directory to store quality trimming data 

ADP_REM_BB="$PROJECT_DIR/01_preprocessing/qualitytrim"


##REFERENCE 
REF="/home/awijeratne/class/20240127_NRT/subsetAra/02_alginment/genome/Arabidopsis_halleri.Ahal2.2.dna.toplevel.fa"


## ALIGNMENT RESULT
RESULTS_DIR=$PROJECT_DIR/02_alginment

GATK_R1=$PROJECT_DIR/03_gatk_round1
GVCF_FILES=$GATK_R1/vcf/



##----------------------------------------------------------------------------------------


#gatk CombineGVCFs --reference $REF --output $RESULTS_DIR/vcf/combined_avr.g.vcf --variant ERR1760144_final_hapcall.g.vcf --variant ERR1760146_final_hapcall.g.vcf --variant ERR1760145_final_hapcall.g.vcf  --variant ERR1760147_final_hapcall.g.vcf

gatk CombineGVCFs --reference $REF --output $GVCF_FILES/combined_avr.g.vcf --variant $GVCF_FILES/ERR1760144_final_hapcall.g.vcf --variant $GVCF_FILES/ERR1760146_final_hapcall.g.vcf --variant $GVCF_FILES/ERR1760145_final_hapcall.g.vcf  --variant $GVCF_FILES/ERR1760147_final_hapcall.g.vcf