#!/bin/bash

#SBATCH --job-name=20240124_QTRIM
#SBATCH --partition=comp01
#SBATCH --output=prepro_qtrim_%j.txt
#SBATCH --error=prepro_qtrim_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=john.doe@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00



############ Cluster-specific parameters and modules ##############

# Set the number of threads for multi-threaded tasks
export OMP_NUM_THREADS=24

#Modules
module load bbmap/38.82
module load java/sunjdk_11.0.10
module load gcc-11.2.1/SKYLAKEX/fastqc/0.11.9
module load python/anaconda-3.12
conda activate multiqc-3.12





#################### CONFIGURATION INFORMATION ####################

# Location of the file containing the list of sample names
file="/home/awijeratne/class/20240127_NRT/subsetAra/data/sample.name.txt"

# Main project directory
PROJECT_DIR=/home/awijeratne/class/20240127_NRT/subsetAra

# Directory containing raw data
RAW_DATA=/home/awijeratne/class/20240127_NRT/subsetAra/data

# Directory to store FASTQC results before quality trimming
FASTQC_B_QTRIM=$PROJECT_DIR/01_preprocessing/fastqc_beforeQtrim

# Directory to store quality trimming data 

ADP_REM_BB="$PROJECT_DIR/01_preprocessing/qualitytrim"

# Directory to store FASTQC results after quality trimming
A_Q_T="$PROJECT_DIR/01_preprocessing/fastqc_afterQtrim"


##################### SOFTWARE SETTINGS ############################

FASTQ_TYPE="sanger"
# Minimum quality score
QUALITY_THRESHOLD=20
# Minimum read length after trimming
MIN_LENGTH=50

ADAPTER_SEQ="/scrfs/storage/awijeratne/home/Projects/20210504_sequencing/01_Preproprocessing/adapter.fa"

# conda activate cutadaptenv

####################################################################

while IFS=" " read -r value1
do

    # Paths to the first and second sample files for the current value1
    FIRST_SAMPLE_LOC=${RAW_DATA}/${value1}_1.fastq.gz
    SECCOND_SAMPLE_LOC=${RAW_DATA}/${value1}_2.fastq.gz

    # Print information about the current sample

    echo "First Sample Location: $FIRST_SAMPLE_LOC"
    echo "Second Sample Location: $SECOND_SAMPLE_LOC"

    ############################## Fastqc ##############################

    # Change to the directory where FASTQC results will be stored
    cd "$FASTQC_B_QTRIM"

    # Run FASTQC on the first and second sample files, storing results in $FASTQC_B_QTRIM
  
    fastqc --threads 32 "$FIRST_SAMPLE_LOC" "$SECOND_SAMPLE_LOC" -o "$FASTQC_B_QTRIM"


    ######################## Adapter Removing ########################

    cd "$ADP_REM_BB"

    #bbduk.sh in1=$FIRST_SAMPLE_LOC in2=$SECOND_SAMPLE_LOC out1=${value1}_bb_trim_R1.fastq out2=${value1}_bb_trim_R2.fastq ref=$ADAPTER_SEQ ktrim=r k=23 mink=11 hdist=1 ftr=100 tpe tbo qtrim=rl trimq=$QUALIY_THRESHOLD minlen=$MIN_LENGTH
    
bbduk.sh in1=$FIRST_SAMPLE_LOC in2=$SECOND_SAMPLE_LOC out1=${value1}_bb_trim_R1.fastq out2=${value1}_bb_trim_R2.fastq ref=$ADAPTER_SEQ ktrim=r k=23 mink=11 hdist=1 ftr=100 ftl=0 tpe tbo qtrim=rl trimq=$QUALITY_THRESHOLD minlen=$MIN_LENGTH
# "ktrim=r” is for right-trimming (3′ adapters)
# “k” specifies the kmer size to use (23)
# “mink” allows it to use shorter kmers at the ends of the read (for example, k=11 for the last 11 bases)
# “hdist” means “hamming distance”; this allows 1 mismatch
# "frt" means trim the right end, retaining the value and lower (ftr=75 means keep positions 75 and lower)
# For normal paired-end fragment libraries, I (bbduk folks) recommend adding the flags “tbo”, which specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences), and “tpe”, which specifies to trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).
# “qtrim=rl” quality-trim to Q10 using the Phred algorithm on both sides

    
    
    ############################ Fastqcaftertrimming ####################

    cd "$A_Q_T"

    fastqc --threads 32 "$ADP_REM_BB/${value1}_bb_trim_R1.fastq" "$ADP_REM_BB/${value1}_bb_trim_R2.fastq" -o "$A_Q_T"

done < "$file"


##run multiQC on the fastQC output

cd $FASTQC_B_QTRIM

multiqc ./ 


cd $A_Q_T

multiqc ./ 
