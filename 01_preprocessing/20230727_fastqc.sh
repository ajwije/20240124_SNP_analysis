#!/bin/bash

############cluster specific parameters modules#############################
# Load required modules (Java and FastQC)


module load java/sunjdk_11.0.10

module load gcc-11.2.1/SKYLAKEX/fastqc/0.11.9

module load python/anaconda-3.12
conda activate multiqc-3.12




# Set the number of threads for multi-threaded tasks
export OMP_NUM_THREADS=32

####################################CONFIGERATION INFORMATION#########################################################
###LOCATIONS AND SAMPLE NAMES#########################################################################################

# Location of the file containing the list of sample names
file="/home/awijeratne/class/20240127_NRT/subsetAra/data/sample.name.txt"

# Main project directory
PROJECT_DIR=/home/awijeratne/class/20240127_NRT/subsetAra

# Directory containing raw data
RAW_DATA=/home/awijeratne/class/20240127_NRT/subsetAra/data

# Directory to store FASTQC results after quality trimming
FASTQC_B_QTRIM=$PROJECT_DIR/01_preprocessing/fastqc_beforeQtrim


###################################SOFTWARE SETTINGS##################################################################
###Q_TRIM AND ADAPTER REMOVING#####
# (Add the relevant commands or tools for quality trimming and adapter removal here.)

#######################################################################################################################
# Loop through each line in the file_names.txt file and process each sample
while IFS=" " read -r value1

do {

    ###################################Raw data#####################################################################

    # Paths to the first and second sample files for the current value1
    FIRST_SAMPLE_LOC=${RAW_DATA}/${value1}_1.fastq.gz
    SECCOND_SAMPLE_LOC=${RAW_DATA}/${value1}_2.fastq.gz

    # Print the first sample location (for demonstration purposes)
    echo $FIRST_SAMPLE_LOC

    ####################################Fastqc####################################################################
    # Change to the directory where FASTQC results will be stored
    cd $FASTQC_B_QTRIM

    # Run FASTQC on the first and second sample files, storing results in $FASTQC_B_QTRIM
    fastqc $FIRST_SAMPLE_LOC --threads 32 --outdir $FASTQC_B_QTRIM
    fastqc $SECCOND_SAMPLE_LOC --threads 32 --outdir $FASTQC_B_QTRIM

} done <"$file"

##run multiQC on the fastQC output

cd $FASTQC_B_QTRIM

multiqc ./ 

conda deactivate

