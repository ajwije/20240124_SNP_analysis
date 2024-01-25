#Preprocessing of reads

## Interactive mode:
bash 20240125_quality_bbtrim.sh

## None interactive slurm script (need to change information related to node etc.)
sbatch 20240125_quality_bbtrim.slurm


This script is designed to perform quality control and preprocessing on sequencing data. It uses several modules and tools including BBMap, Java, GCC, FastQC, Python, and MultiQC.

Cluster-specific parameters and modules: The script sets the number of threads for multi-threaded tasks and loads necessary modules.

Configuration Information: The script sets up directories for the project, raw data, and results. It also specifies the location of the file containing the list of sample names.

Software Settings: The script sets parameters for the quality control process, such as the type of FASTQ files, the minimum quality score, the minimum read length after trimming, and the sequence of the adapter to be removed.

Main Loop: The script reads the file containing the list of sample names line by line. For each sample, it performs the following steps:

FastQC Before Trimming: Runs FastQC on the raw data files to assess their quality before trimming.
Adapter Removal: Uses BBduk (from the BBMap suite) to trim adapters and low-quality ends from the reads.
FastQC After Trimming: Runs FastQC again on the trimmed data files to assess their quality after trimming.
MultiQC: Finally, the script runs MultiQC to aggregate the FastQC results into a single report. This is done for both the pre-trimmed and post-trimmed data.

Please note that this script assumes a paired-end sequencing data format (two FASTQ files per sample). The script is designed to be run on a cluster that uses environment modules to manage software. Make sure to replace the paths and parameters with those that match your specific setup.