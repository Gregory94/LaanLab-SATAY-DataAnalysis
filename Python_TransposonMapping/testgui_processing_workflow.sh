#!/bin/bash

### This is test version of the processing_workflow.sh that is used for developing a gui using zenity.



#filepath=`zenity --file-selection --file-filter="*.fq" --file-filter="*.fastq" --file-directory="/home/laanlab/Documents/satay"`
#zenity --warning --title="file selection" --text="selected file is "$filepath
#
#settings=`zenity --entry --text="Input settings"`
#zenity --info --title="Input settings" --text=$settings
#
#settings=`zenity --forms --title="Processing settings" --text="Data type" --add-entry="Using paired end data?" --add-entry="Enter settings for BBDuk"`
#echo 'settings are '$settings
#
fileselection1=`yad --width=1000 --height=400 --title="Select fastq file" --file-selection="Please select datafile" --file-filter="*.fq" --file-filter="*.fastq"`
#fileselection2=`yad --width=1000 --height=400 --title="Select fastq file" --file-selection="Please select datafile" --file-filter="*.fq" --file-filter="*fastq" --file-directory="/home/laanlab/Documents/satay"`
settings=`yad --width=1000 --height=800 --title="Processing settings" --text="Settings" --form --field="Selected file primary reads":RO --field="Enter data type":CB --field="Enter trimming settings" --field="Enter alignment settings":CHK $fileselection1 "single-read!paired-end" "" ""`
echo $settings

####################### USER SETTINGS ######################
# Define whether data is paired-end ('t' for paired-end, 'f' for single end)
paired=T

# Define filename (can also be a zipped file ending with .gz). Use filename2 for paired end or leave empty for single end or interleaved paired end (i.e. paired end reads are in single file).
filepath1=/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/wt1_enzo_dataset_demultiplexed_interleaved_sample1/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_pairs.fq
#filepath1=/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/wt1_enzo_dataset_demultiplexed_interleaved_sample2/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_pairs.fq
filepath2=''


###### Set options for trimming software ######
# set 'b' for bbduk, set 't' for trimmomatic
trimming_software='b'

###    bbduk    ###
trimming_settings_bbduk='ktrim=l k=15 mink=10 hdist=1 tpe tbo qtrim=r trimq=10 minlen=30'
#trimming_settings_bbduk='k=20 mink=8 ktrim=l restrictleft=50 hdist=3 hdist2=1 qtrim=r trimq=10 minlen=25 tpe=t tbo=t'
## Set adapter sequences
## Open file using xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa
###################

### trimmomatic ###
trimmomatic_initialization='-phred33'
trimming_settings_trimmomatic='ILLUMINACLIP:adapters.fa:0:30:10 SLIDINGWINDOW:10:4 MINLEN:30'
## Set adapter sequences
## Open file using xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa
###################
###############################################


# Set options for alignment software (bwa mem) (note that for paired end data the parameter -p does not need to be set as long as paired=T)
alignment_settings='-M -B 3 -O 3,3 -S -v 2'

# Trim the reads with the options set in trimming software section above ('T' for yes, 'F' for no)?
trimming=T

# Create sorted and indexed bam file ('T' for yes, 'F' for no)?
sort_and_index=T


# Apply transposon mapping (requires sort_and_index=T)
mapping=F


# Delete sam file ('T' for yes, 'F' for no)? This file is always converted to its binary equivalent (.bam ) and the sam file is rarely used but takes up relatively a lot of memory.
delete_sam=F


# Create a quality report of the alignment based on the sam file (this also works when the sam file is being deleted, i.e delete_sam=T)
flagstat_report=T


# Open adapters.fa file after the first quality check in order to change the adapters for trimming.
#open_adapters=F


# Create quality report of raw data (before trimming)?
quality_check_raw=F


# Create quality report of trimmed data (after trimming)?
quality_check_trim=F


# Determine whether the script should automatically continue after creating the first quality report. Set to True if you might want to make changes depending on the quality report of the raw data.
qualitycheck_interrupt=F

############################################################



echo 'Processing finished'
