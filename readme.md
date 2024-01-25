# Before you start:
##Download data:
- You can get the data from the following Dropbox folder: https://www.dropbox.com/scl/fi/ot9cg0udzqv9cunqufjzq/data.tar.gz?rlkey=5jya8v3qyxcnukm90tbg4gxx7&dl=0

##Create a folder in your working directory:

\```bash
mkdir foldername
\```

Then run the following script with folders you want to create:

\```bash
bash folder_structure.sh 01_preprocessing 01_preprocessing/fastqc_beforeQtrim 01_preprocessing/fastqc_afterQtrim 01_preprocessing/qualitytrim 02_alginment 03_gatk_round1 04_gatk_round2
\```

## Preprocessing of reads

### 01_preprocessing/20240125_quality_bbtrim.sh

Run the script:

Interactive mode:

\```bash
bash 20240125_quality_bbtrim.sh
\```

Non-interactive Slurm script (need to change information related to node, etc.):

\```bash
sbatch 20240125_quality_bbtrim.slurm
\```

This script is designed to perform quality control and preprocessing on sequencing data:

- **Cluster-specific parameters and modules:** The script sets the number of threads for multi-threaded tasks and loads necessary modules.
- **Configuration Information:** Sets up directories for the project, raw data, and results. Specifies the location of the file containing the list of sample names.
- **Software Settings:** Sets parameters for the quality control process.
- **Main Loop:** Reads the file containing the list of sample names line by line and performs FastQC before trimming, adapter removal, FastQC after trimming.
- **MultiQC:** Aggregates FastQC results into a single report for both pre-trimmed and post-trimmed data.

## Alignment step

### 02_alginment/20240124_bwa.sh

- **Load Modules:** Loads necessary modules and activates the MultiQC environment.
- **Configuration Information:** Sets up directories for the project, raw data, and results. Specifies the location of the file containing the list of sample names.
- **Index Reference Genome:** Indexes the reference genome using BWA. This step is commented out and only needs to be done once per reference.
- **Main Loop:** Reads the file containing the list of sample names line by line and performs alignment, SAM to BAM conversion, sorting, indexing, mapping stats, and FastQC. Adds read group header to the SAM/BAM file.
- **MultiQC:** Aggregates FastQC results into a single report.

## Genotyping

### Initial genotyping with each file:

#### 03_gatk_round1/20240125_gatk_r1.sh

\```bash
bash 20240125_gatk_r1.sh
\```

- **Configuration Information:** Defines paths to various files and directories. Prepares the reference genome for use with GATK. Defines the directory for alignment results and creates a directory for storing the output of GATK’s first round of variant calling.
- **Loop Through Samples:** Marks duplicates in the aligned and sorted BAM file using Picard’s MarkDuplicates tool (commented out). Indexes the BAM file using Samtools’ index tool. Calls variants on the BAM file using GATK’s HaplotypeCaller tool, outputting a GVCF file.

### Combine genotyping:

#### 03_gatk_round1/20240125_gatk_comb_genotype.sh

- **Configuration Information:** Defines paths to various files and directories. Defines the directory for alignment results and creates a directory for storing the output of GATK’s first round of variant calling.
- **Combine GVCFs:** Runs GATK’s CombineGVCFs tool to combine GVCF files into a single GVCF file.

### SNP filtering

#### 04_gatk_round2/20240125_gatk_filter_genotype.sh

- **Configuration Information:** Defines paths to various files and directories. Defines the directory for alignment results and creates a directory for storing the output of GATK’s first round of variant calling.
- **Steps:** Combine GVCFs: Runs GATK’s CombineGVCFs tool to combine GVCF files into a single GVCF file. Genotype the Combined GVCFs: Runs GATK’s GenotypeGVCFs tool to perform joint genotyping on the combined GVCF file. Select Variants: Runs GATK’s SelectVariants tool to select only SNP variants from the genotyped VCF file. Filter Variants: Runs GATK’s VariantFiltration tool to apply quality filters to the SNP variants.

    - The `--reference` option specifies the reference genome.
    - The `--output` option specifies the output file.
    - The `--variant` option specifies the input VCF file.
    - The `--filter-name` and `--filter-expression` options specify the filters to be applied.
