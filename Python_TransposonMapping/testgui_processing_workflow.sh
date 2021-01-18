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

cachefile="/home/laanlab/Documents/satay/software/processing_workflow_cache.txt"

if [ ! -f $cachefile ];
then
	fileselections=`yad --width=1000 --height=400 --title="Select fastq file" --center --on-top --buttons-layout=spread --multiple --file-selection="Please select datafile (or two files in case of paired-end noninterleaved fastq files)" --file-filter="*.fq" --file-filter="*.fastq" --file-filter="*.fq.gz" --file-filter="*.fastq.gz"`
	filepath1=$(echo $fileselections | awk 'BEGIN {FS="|" } { print $1 }')
	filepath2=$(echo $fileselections | awk 'BEGIN {FS="|" } { print $2 }')
	[ -z $filepath1 ] && filepath1='none'
	[ -z $filepath2 ] && filepath2='none'

	settings=`yad --width=1000 --height=500 --title="Processing settings" --text="Settings" --center --on-top --buttons-layout=spread --form \
	--field="Selected file primary reads":RO \
	--field="Selected file secondary reads":RO \
	--field="Data type":CB \
	--field="Which trimming to use":CB \
	--field="Enter trimming settings" \
	--field="Enter alignment settings" \
	--field="Quality checking raw data":CHK \
	--field="Quality checking trimmed data":CHK \
	--field="Quality check interrupt\n (allows for changing trimming and alignment settings after quality report raw data)":CHK \
	--field="Delete sam file":CHK \
	--field="Sort and index bam files":CHK \
	--field="Transposon mapping (NOTE: requires sorting and indexing)":CHK \
	--field="Create flagstat report":CHK \
	--field="Open adapters file":FBTN \
	$filepath1 \
	$filepath2 \
	"Paired-end!Single-end" \
	"bbduk!trimmomatic!Do not trim" \
	"ktrim=l k=15 mink=10 hdist=1 tpe tbo qtrim=r trimq=10 minlen=30" \
	" -p -M -S -P -v 2" \
	"False" \
	"TRUE" \
	"FALSE" \
	"TRUE" \
	"TRUE" \
	"TRUE" \
	"TRUE" \
	"bash -c 'xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa'"`

	if [ ! -z "$settings" ] && [ $filepath1 != "none" ]
	then
		echo $settings >> $cachefile
	fi

elif [ -f $cachefile ];
then

	previoussettings=`head -n 1 $cachefile`
	echo $previoussettings

	settings=`yad --width=1000 --height=500 --title="Processing settings" --text="Settings" --center --on-top --buttons-layout=spread --form \
	--field="Selected file primary reads":RO \
	--field="Selected file secondary reads":RO \
	--field="Data type":RO \
	--field="Which trimming to use":RO \
	--field="Enter trimming settings" \
	--field="Enter alignment settings" \
	--field="Quality checking raw data":CHK \
	--field="Quality checking trimmed data":CHK \
	--field="Quality check interrupt\n (allows for changing trimming and alignment settings after quality report raw data)":CHK \
	--field="Delete sam file":CHK \
	--field="Sort and index bam files":CHK \
	--field="Transposon mapping (NOTE: requires sorting and indexing)":CHK \
	--field="Create flagstat report":CHK \
	--field="Open adapters file":FBTN \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $1 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $2 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $3 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $4 }') \
	"$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $5 }')" \
	"$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $6 }')" \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $7 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $8 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $9 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $10 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $11 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $12 }') \
	$(echo $previoussettings | awk 'BEGIN {FS="|" } { print $13 }') \
	"bash -c 'xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa'"`

	rm $cachefile
fi


####################### USER SETTINGS ######################
# Define whether data is paired-end ('t' for paired-end, 'f' for single end)
paired=$(echo $settings | awk 'BEGIN {FS="|" } { print $3 }')
echo 'paired ' $paired

###### Set options for trimming software ######
# set 'b' for bbduk, set 't' for trimmomatic
trimming_software=$(echo $settings | awk 'BEGIN {FS="|" } { print $4 }')
echo 'trimming_software ' $trimming_software

###    bbduk    ###
trimming_settings_bbduk=$(echo $settings | awk 'BEGIN {FS="|" } { print $5 }')
echo 'trimming_settings_bbduk ' $trimming_settings_bbduk
#trimming_settings_bbduk='k=20 mink=8 ktrim=l restrictleft=50 hdist=3 hdist2=1 qtrim=r trimq=10 minlen=25 tpe=t tbo=t'
## Set adapter sequences
## Open file using xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa
###################

### trimmomatic ###
trimmomatic_initialization='-phred33'
trimming_settings_trimmomatic=$(echo $settings | awk 'BEGIN {FS="|" } { print $5 }')
echo 'trimming_settings_trimmomatic ' $trimming_settings_trimmomatic
## Set adapter sequences
## Open file using xdg-open /home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa
###################
###############################################


# Set options for alignment software (bwa mem) (note that for paired end data the parameter -p does not need to be set as long as paired=T)
alignment_settings=$(echo $settings | awk 'BEGIN {FS="|" } { print $6 }')
echo 'alignment_settings ' $alignment_settings

# Trim the reads with the options set in trimming software section above ('T' for yes, 'F' for no)?
#trimming=T

# Create sorted and indexed bam file ('T' for yes, 'F' for no)?
sort_and_index=$(echo $settings | awk 'BEGIN {FS="|" } { print $11 }')
echo 'sort_and_index ' $sort_and_index


# Apply transposon mapping (requires sort_and_index=T)
mapping=$(echo $settings | awk 'BEGIN {FS="|" } { print $12 }')
echo 'mapping ' $mapping


# Delete sam file ('T' for yes, 'F' for no)? This file is always converted to its binary equivalent (.bam ) and the sam file is rarely used but takes up relatively a lot of memory.
delete_sam=$(echo $settings | awk 'BEGIN {FS="|" } { print $10 }')
echo 'delete_sam ' $delete_sam


# Create a quality report of the alignment based on the sam file (this also works when the sam file is being deleted, i.e delete_sam=T)
flagstat_report=$(echo $settings | awk 'BEGIN {FS="|" } { print $13 }')
echo 'flagstat_report ' $flagstat_report


# Open adapters.fa file after the first quality check in order to change the adapters for trimming.
#open_adapters=F


# Create quality report of raw data (before trimming)?
quality_check_raw=$(echo $settings | awk 'BEGIN {FS="|" } { print $7 }')
echo 'quality_check_raw ' $quality_check_raw


# Create quality report of trimmed data (after trimming)?
quality_check_trim=$(echo $settings | awk 'BEGIN {FS="|" } { print $8 }')
echo 'quality_check_trim ' $quality_check_trim


# Determine whether the script should automatically continue after creating the first quality report. Set to True if you might want to make changes depending on the quality report of the raw data.
qualitycheck_interrupt=$(echo $settings | awk 'BEGIN {FS="|" } { print $9 }')
echo 'qualitycheck_interrupt ' $qualitycheck_interrupt

############################################################

#After interrupt, rm $cachefile

echo 'Processing finished'
