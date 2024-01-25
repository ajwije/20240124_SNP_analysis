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

#PREPARE GENOME FOR GATK 
gatk CreateSequenceDictionary -R $REF
samtools faidx $REF


## ALIGNMENT RESULT
RESULTS_DIR=$PROJECT_DIR/02_alginment

GATK_R1=$PROJECT_DIR/03_gatk_round1
mkdir $GATK_R1/vcf



##----------------------------------------------------------------------------------------
## Loop through samples
## mark duplicates, index bams, call variants

while IFS=" " read -r value1
do {


# mark duplicates in aligned_sorted.bam files
# this marks duplicates so that GATK tools will ignore them. They are not removed from the dataset unless requested.
#java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar MarkDuplicates INPUT=${RESULTS_DIR}/bam/${value1}_align_sort.bam OUTPUT=${RESULTS_DIR}/bam/${value1}_final.bam METRICS_FILE=${RESULTS_DIR}/bam/${value1}_final_metrics.txt

# index *final.bam files (i.e., aligned, sorted, duplMarked)
samtools index $RESULTS_DIR/bam/${value1}_final.bam

####



# variant calling
# run HaplotypeCaller for each sample individually, in GVCF mode (for jointly calling downstream, that is the --emit-ref-confidence GVCF flag))
gatk HaplotypeCaller --reference $REF --input $RESULTS_DIR/bam/${value1}_final.bam --output $GATK_R1/vcf/${value1}_final_hapcall.g.vcf --emit-ref-confidence GVCF

} done <"$file"

