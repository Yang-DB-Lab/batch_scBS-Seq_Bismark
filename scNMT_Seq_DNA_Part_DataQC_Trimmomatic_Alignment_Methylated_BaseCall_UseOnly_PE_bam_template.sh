#!/bin/bash

# folder with all folders of fastq files
top_level_folder="/data/guang/Experiment_20211208/raw_data/scNMT_Seq_DNA"
######### Parameters related to Trim-Galore ##########
# How to run Trimmomatic (inlcude path of Trimmomatic)
execute_Trimmomatic="java -jar /home/guang/bio_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"

## Files with adapter sequences
# Paired-end Trimmomatic TruSeq3 adapter file
TruSeq3_PE__adapter="/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
# Single-end Trimmomatic TruSeq3 adapter file
TruSeq3_SE_adapter="/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
# Paired-end Trimmomatic NexteraPE-PE adapter file
NexteraPE_PE_adapter="/home/guang/bio_softwares/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
######################################################

######### Parameters related to Bismark ###############
genome_Bismark_index_Dir="/home/guang/mouse_genome_index/Bismark_GRCm38.101_with_Lambda_DNA"
######### Prepare genome index for Bismark if it is not available
## Check whether Bismark genome index folder(-d option) available
# if [ ! -d "$genome_Bismark_index_Dir" ]
# then 
# 	echo "Bismark genome index does not exist. Build it from beginning."
	## Build Bismark index
	
# 	mkdir -p $genome_Bismark_index_Dir
	# wget -P $genome_Bismark_index_Dir http://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	# cd $genome_Bismark_index_Dir
	# gunzip  Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	# ls $genome_Bismark_index_Dir
	## Mus_musculus.GRCm38.dna.primary_assembly.fa  Mus_musculus.GRCm38.dna.primary_assembly.fa.gz lambda_DNA.fasta
	# cat lambda_DNA.fasta >> Mus_musculus.GRCm38.dna.primary_assembly.fa
	# mv Mus_musculus.GRCm38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly_with_lambda_DNA.fa
	# rm lambda_DNA.fasta 
	# ls $genome_Bismark_index_Dir
	## Mus_musculus.GRCm38.dna.primary_assembly.fa.gz  Mus_musculus.GRCm38.dna.primary_assembly_with_lambda_DNA.fa

	# Convert C to T using bismark_genome_preparation
# 	bismark_genome_preparation --path_to_aligner /home/guang/bio_softwares/bowtie2-2.3.5.1-linux-x86_64/ --verbose $genome_Bismark_index_Dir
# fi

##

##******************* Parameter Definition Finished Here *********************##

##*************** The actual alignment script starts here. ********************##
################################################################################

## Global variables used by all script pieces.
## The global variables will be re-defined in each part. Although this re-definition
## is not necessary for execution of this merged script, re-definition of these
## global variables will make each part still a complete script and can be copied
## out to run independently.
## Re-define of the global variables to make  this part independent.
## top_level_folder="/data/guang/Experiment_20211208/raw_data/scNMT_Seq_DNA"
## Initialize the top_level folder. Must use full path!!
## This folder should contain the sample folders with single-end fatq.gz files.
## This folder will be used by all script pieces below.
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"
## Get names of sample folders:
## The command on the right site will first use 'ls -l' to check all files and folders
## then 'grep "^d" ' will select the ones will a 'd' property at the beginning of 
## 'ls -l' command which are folders
## finally 'awk '{print $NF}' will print out the last column of 'ls -l' which are 
## the actual folder names. 
## NOTE: There should be no space in the folder names. If space exists in folder names
##       only the last word of folder name will be print out. This is not what we want!
## sample_folder_names will be used by all script pieces below.

################## Part 1 Quality Control and fastq data trimming.#####################
############# Part 1.1 Check fastq data quality using FastQC ##########################

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	# Get into the sample_folder with fastq.gz file(s)
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
	# Make a new folder in sample_folder to store FastQC result; Full path used here.
	fastqc -t 8 *.fq.gz -O ./$sample_folder\_fastqc_results/
	# -t 8: use 8 threads
	# relative path is used. The full path is 
	# $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
	# Two files generate after calling fastqc: 
	#                (fastq.gz filename)_fastqc.html and (fastq.gz filename)_fastqc.zip
done



########################################################################################
################# Trimming use Trimmomatic
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	# Get into the sample_folder with fastq.gz file
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_Trimmomatic_trimmed
	# Make a new folder in sample_folder to store trimmed result; Full path used here.
	fastq_gz_files="$(ls *.fq.gz)"
	
	$execute_Trimmomatic PE -threads 8 -phred33 \
	$fastq_gz_files \
	./$sample_folder\_Trimmomatic_trimmed/trimmed.$sample_folder\_paired.R1.fq.gz \
	./$sample_folder\_Trimmomatic_trimmed/trimmed.$sample_folder\_unpaired.R1.fq.gz \
	./$sample_folder\_Trimmomatic_trimmed/trimmed.$sample_folder\_paired.R2.fq.gz \
	./$sample_folder\_Trimmomatic_trimmed/trimmed.$sample_folder\_unpaired.R2.fq.gz \
	ILLUMINACLIP:$NexteraPE_PE_adapter:2:30:10 LEADING:3 \
	TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	
	# remove old fastq files to save space
	# rm *.gz
done




######### Parameters related to Bismark ###############
# bismark genome index is generated in above parameter part
genome_Bismark_index_Dir="/home/guang/mouse_genome_index/Bismark_GRCm38.101_with_Lambda_DNA/"
####################################################

# top_level_folder="/data/guang/Experiment_20211208/raw_data/scNMT_Seq_DNA"

cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

#### 2.1 align the reads using bismark
for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	#get into sample_folder
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_bismark_results
	
	trimmed_fastq_gz_files="$(ls -d $PWD/$sample_folder\_*\_trimmed/*.fq.gz)"
	#$PWD is current path
	#This command get the full path of trimmed fast.gz files
	
	paired_trimmed_read1="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_paired.R1.fq.gz)"
	# paired_trimmed_read1_basename="$(basename $paired_trimmed_read1)"
	paired_trimmed_read2="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_paired.R2.fq.gz)"
	# paired_trimmed_read2_basename="$(basename $paired_trimmed_read2)"
	unpaired_trimmed_reads="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_unpaired*.fq.gz)"
	# use unpaired_reads to store the filenames of the 2 unpaired read files trimmed out by trimmomatic
	
	# process trim galore trimmed paired-end read files
	# Bowtie2 is used as the aligner which is default
	cd $top_level_folder/$sample_folder/$sample_folder\_bismark_results
	# cd: get into the folder where the results will be put their
	
	# turn on "--non_directional" parameter when do alignment
	# add "--dovetail" to make overlapped pair-mates legal
	bismark --non_directional --dovetail --maxins 5000 --genome $genome_Bismark_index_Dir \
			-1 $paired_trimmed_read1 -2 $paired_trimmed_read2

#### 2.2 bismark deduplicate

	deduplicate_bismark --paired *bt2_pe.bam
	
	
#### 2.3 Bismark methylation extractor
	# bismark_methylation_extractor --gzip --bedGraph *.deduplicated.bam
	bismark_methylation_extractor --gzip --bedGraph \
	--buffer_size 10G --cytosine_report \
	--genome_folder $genome_Bismark_index_Dir \
	*.deduplicated.bam
	
#### 2.4 Use coverage2cytosine module for scNMT-Seq
	coverage2cytosine --nome-seq
	# turn on the  "--nome-seq" paramter for scNMT-seq
	
	# The option --nome-seq:
	# (i) filters the genome-wide CpG-report to only output cytosines in ACG and TCG context
	# (ii) filters the GC context output to only report cytosines in GCA, GCC and GCT context
	
	# PLEASE NOTE that NOMe-seq data requires a .cov.gz file as input which has been generated 
	# in non-CG mode!! (--CX), else the GpC output file will be empty...
#### 2.5 The Bismark HTML Processing Report
	bismark2report
#### 2.6 The Bismark Summary Report
	bismark2summary
#### 2.7 Bismark Nucleotide Coverage report 
	bam2nuc --genome_folder $genome_Bismark_index_Dir 
	# genome_folder contains the .fa genome sequence file
#### 2.8 Filtering out non-bisulfite converted reads
	filter_non_conversion
done
##****************************Part2 Finished here *******************************##





