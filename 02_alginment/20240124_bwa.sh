#!/bin/bash


##----------------------------------------------------------------------------------------
## LOAD MODULES

module load bwa/0.7.17
module load samtools/1.15.1
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load java/sunjdk_1.8.0 fastqc/0.12.1 picard-tools/2.17.10
module load python/anaconda-3.12
conda activate multiqc-3.12

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



REF=/home/awijeratne/class/20240127_NRT/subsetAra/02_alginment/genome/Arabidopsis_halleri.Ahal2.2.dna.toplevel.fa




## ALIGNMENT RESULT DIR PER SPECIES
RESULTS_DIR=$PROJECT_DIR/02_alginment

## MAKE SUBDIRS IN RESULTS_DIR
mkdir ${RESULTS_DIR}/bam
mkdir ${RESULTS_DIR}/sam
mkdir ${RESULTS_DIR}/flagstat
mkdir ${RESULTS_DIR}/fastqc



##----------------------------------------------------------------------------------------
## Index reference genome 
## only have to do this once per reference
#bwa index $REF


##----------------------------------------------------------------------------------------
## LOOP THROUGH SAMPLES

while IFS=" " read -r value1
do {

# Paths to the first and second sample files for the current value1
    FIRST_SAMPLE_LOC=${RAW_DATA}/${value1}_1.fastq
    SECCOND_SAMPLE_LOC=${RAW_DATA}/${value1}_2.fastq

    # Print information about the current sample

    echo "First Sample Location: $FIRST_SAMPLE_LOC"
    echo "Second Sample Location: $SECOND_SAMPLE_LOC"
    
    
# align fastq reads to reference genome using BWA, output is a .sam file
# -M: tells bwa to consider split reads as secondary
# -R: read group info
bwa mem -t 32 -M -R "@RG\tID:${value1}\tSM:${value1}\tPL:ILLUMINA" $REF "$FIRST_SAMPLE_LOC" $SECCOND_SAMPLE_LOC > ${RESULTS_DIR}/sam/${value1}_align.sam

# convert .sam to .bam file (bam is a compressed binary version of sam) using samtools
samtools view --threads 30 -S -b ${RESULTS_DIR}/sam/${value1}_align.sam > ${RESULTS_DIR}/bam/${value1}_align.bam

# sort the bam files using samtools
samtools sort --threads 30 ${RESULTS_DIR}/bam/${value1}_align.bam -o ${RESULTS_DIR}/bam/${value1}_align_sort.bam 

# index the align_sort.bam files
samtools index ${RESULTS_DIR}/bam/${value1}_align_sort.bam

# mapping stats
samtools flagstat --threads 30 -O tsv ${RESULTS_DIR}/bam/${value1}_align_sort.bam > ${RESULTS_DIR}/flagstat/${value1}_align_sort_flagstat_out.txt

fastqc --threads 32 ${RESULTS_DIR}/bam/${value1}_align_sort.bam -o ${RESULTS_DIR}/fastqc
} done <"$file"

cd ${RESULTS_DIR}/fastqc

multiqc ./
