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

GATK_R2=$PROJECT_DIR/04_gatk_round2




##----------------------------------------------------------------------------------------


#gatk CombineGVCFs --reference $REF --output $RESULTS_DIR/vcf/combined_avr.g.vcf --variant ERR1760144_final_hapcall.g.vcf --variant ERR1760146_final_hapcall.g.vcf --variant ERR1760145_final_hapcall.g.vcf  --variant ERR1760147_final_hapcall.g.vcf

gatk CombineGVCFs --reference $REF --output $GVCF_FILES/combined_avr.g.vcf --variant $GVCF_FILES/ERR1760144_final_hapcall.g.vcf --variant $GVCF_FILES/ERR1760146_final_hapcall.g.vcf --variant $GVCF_FILES/ERR1760145_final_hapcall.g.vcf  --variant $GVCF_FILES/ERR1760147_final_hapcall.g.vcf

##----------------------------------------------------------------------------------------
## GENOTYPE THE COMBINED GVCFs

gatk GenotypeGVCFs --reference $REF --output $GATK_R2/ERR_comb_genotyped.vcf --variant $GVCF_FILES/ERR1760144_final_hapcall.g.vcf



##----------------------------------------------------------------------------------------
## SELECT VARIANTS (using output from GenotypeGVCFs, i.e., genotyped.vcf)
# select variants first, then can filter them separately

# SNPs only
gatk SelectVariants --reference $REF --output $GATK_R2/ERR_raw_comb_SNP.vcf --select-type-to-include SNP --variant $GATK_R2/ERR_comb_genotyped.vcf



##----------------------------------------------------------------------------------------
## FILTER VARIANTS (use output from SelectVariants, separate for SNP and INDEL)

# GATK recommendations (meant to be "very lenient") for filters here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# and more recs here: https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html


# SNPs only
gatk VariantFiltration --reference $REF --output $GATK_R2/ERR_filtered_comb_SNP.vcf --variant $GATK_R2/ERR_raw_comb_SNP.vcf --filter-name "QD" --filter-expression "QD < 2.0" --filter-name "FS" --filter-expression "FS > 60.0" --filter-name "SOR" --filter-expression "SOR > 3.0" --filter-name "MQ" --filter-expression "MQ < 40.0" --filter-name "MQRankSum" --filter-expression " MQRankSum < -12.5" --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0"


##----------------------------------------------------------------------------------------
## FROM HERE, YOU CAN... 
# - use SelectVariants on these filtered files to only select SNPs/INDELs that passed all filters (--exclude-filtered true)
# - filter INDEL
# - use MergeVcfs to merge the SNPs and INDELs back into a single vcf file
# - probably lots of other things