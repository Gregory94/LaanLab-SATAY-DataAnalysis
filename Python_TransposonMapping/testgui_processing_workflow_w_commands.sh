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
	"bbduk!trimmomatic!Do_not_trim" \
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

#	if [ ! -z "$settings" ] && [ $filepath1 != "none" ]
#	then
#		echo $settings >> $cachefile
#	fi

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


# Set options for alignment software (bwa mem) (note that for paired end data the parameter -p does not need to be set as long as paired=Paired-end)
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

echo 'Preparing processing for' $(basename ${filepath1}) '...'
echo ''


pathdata=$(dirname ${filepath1})
filename1=$(basename ${filepath1})

if [[ ${filepath1} =~ 'none' ]]
then
	echo 'ERROR: No file selected. Process canceled.' && exit 1
elif [ ! -f ${filepath1} ]
then
	echo 'ERROR: File' ${filepath1} 'not found.' && exit 1
fi


if [[ ${paired} =~ 'Paired-end' ]] && ! [[ ${filepath2} =~ 'none' ]]
then
	if [ ! -f ${filepath2} ]
	then
		echo 'ERROR: File' ${filepath2} 'not found.' && exit 1
	else
		filename2=$(basename ${filepath2})
	fi
elif [[ ${paired} =~ 'Single-end' ]] && ! [[ ${filepath2} =~ 'none' ]]
then
	echo 'WARNING: A secondary reads file was specified but paired was set to single-end. Therefore the secondary reads file will be ignored.'
fi

if [[ -z "${settings}" ]]
then
	echo 'Process canceled.' && exit 1
fi



# Get extension of the file
extension='.'$(echo $filename1 | rev | cut -d. -f1 | rev)
if [[ '$extension' == '.gz' ]]
then
	extension='.'$(echo $filename1 | rev | cut -d. -f2 | rev)
fi



# Define filename for trimming and alignment results
if ! [[ ${trimming_software} =~ 'Do not trim' ]]
then
	filename_trimmed1=${filename1%$extension*}'_trimmed.fastq'
	if ! [[ ${filepath2} =~ 'none' ]] #if not filename2 is none
	then
		filename_trimmed2=${filename2%$extension*}'_trimmed.fastq'
	fi
	filename_sam=${filename1%$extension*}'_trimmed.sam'
	filename_bam=${filename1%$extension*}'_trimmed.bam'
	filename_sort=${filename1%$extension*}'_trimmed.sorted.bam'
else
	filename_trimmed1=${filename1}
	if ! [[ ${filepath2} =~ 'none' ]] #if not filename2 is none
	then
		filename_trimmed2=${filename2}
	fi
	filename_sam=${filename1%$extension*}'_notrimmed.sam'
	filename_bam=${filename1%$extension*}'_notrimmed.bam'
	filename_sort=${filename1%$extension*}'_notrimmed.sorted.bam'
fi



# Define path output directory fastqc
if [[ ${quality_check_raw} =~ TRUE ]] || [[ ${quality_check_trim} =~ TRUE ]]
then
	path_fastqc_out=${pathdata}/fastqc_out
	[ ! -d ${path_fastqc_out} ] && echo 'Creating fastqc output folder ...' && mkdir ${path_fastqc_out} || echo 'Folder for fastqc output exists with path:' ${path_fastqc_out}
fi

exit 1
echo 'This you should not read:('

# Define path output directory trimming
if [[ ${trimming} =~ ^[tT]$ ]]
then
	path_trimm_out=${pathdata}/trimm_out
	[ ! -d ${path_trimm_out} ] && echo 'Creating trimming output folder ...' && mkdir ${path_trimm_out} || echo 'Folder for trimming output exists with name:' $(basename ${path_trimm_out})
else
	path_trimm_out=${pathdata}
fi

# Define path output directory alignment
path_align_out=${pathdata}/align_out
[ ! -d ${path_align_out} ] && echo 'Creating alignment output folder ...' && mkdir ${path_align_out} || echo 'Folder for alignment output exists with name:' $(basename ${path_align_out})

# Define paths to reference genomes (both S288C and W303)
path_refgenome=/home/laanlab/Documents/satay/reference_sequences/Cerevisiae_S288C_reference/S288C_reference_sequence_R64-2-1_20150113.fsa
name_refgenome='S288C'
if [ ! -f ${path_refgenome} ] #if path to reference genome does not exist
then
	echo 'ERROR: Reference genome not found at location:' ${path_refgenome} && exit 1
else
	echo 'Reference genome:' ${name_refgenome}
fi

# Define path bbduk software
path_bbduk_software=/home/laanlab/Documents/satay/software/bbmap/
path_bbduk_adapters=${path_bbduk_software}/resources/adapters.fa
[ ! -d ${path_ddbuk_software} ] && echo 'WARNING: Path to bbduk software does not exists.'

# Define path trimmomatic software
path_trimm_software=/home/laanlab/Documents/satay/software/Trimmomatic-0.39/
[ ! -d ${path_trimm_software} ] && echo 'WARNING: Path to trimmomatic software does not exists.'

# Define path to python script
path_python_codes=/home/laanlab/Documents/satay/software/python_codes/
[ ! -d ${path_python_codes} ] && echo 'WARNING: Path to python codes does not exists.'






### Start processing workflow
echo ''
echo 'Start processing ...'
echo ''


# Quality checking raw data
if [[ ${quality_check_raw} =~ ^[tT]$ ]]
then
	if [[ ! -e ${path_fastqc_out}/${filename1%$extension*}'_fastqc.html' ]]
	then
		echo 'Quality checking raw data ...'
		fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename1}
		echo 'Quality checking raw data completed. Results are stored at' ${path_fastqc_out}
		echo ''

		if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
		then
			fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename2}
			echo 'Quality checking raw data paired end reads completed. Results are stored at' ${path_fastqc_out}
			echo ''
		fi
	else
		echo 'Quality report raw data already exists. Skipping fastqc'
	fi


	if [[ ${qualitycheck_interrupt} =~ ^[tT]$ ]]
	then
		read -p 'Continue processing? (press "y" if yes, press "n" if no): ' -n 1 -r
		echo
		if [[ ! $REPLY =~ ^[yY]$ ]]
		then
			exit 1
		fi
	fi
fi


#if [[ ${open_adapters} =~ ^[tT]$ ]]
#then
#	echo "Adapter.fa file is being opened..."
#	xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
#	read -s -p "Press enter to continue"
#fi


# Trimming
if [[ ${trimming} =~ ^[tT]$ ]]
then
	if [[ ${trimming_software} =~ ^[bB]$ ]]
	then
		if [[ ${paired} =~ ^[fF]$ ]]
		then
			echo 'Data trimming using bbduk single end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
			echo ''
		elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
		then
			echo 'Data trimming using bbduk paired end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g in1=${pathdata}/${filename1} out1=${path_trimm_out}/${filename_trimmed1} in2=${pathdata}/${filename2} out2=${path_trimm_out}/${filename_trimmed2} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1} 'and for the paired end reads in' ${path_trimm_out}/${filename_trimmed2}
			echo ''

		elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filepath2} ]]
		then
			echo 'Data trimming using bbduk paired end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g interleaved=t in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
			echo ''
		fi

	elif [[ ${trimming_software} =~ ^[tT]$ ]]
	then
		if [[ ${paired} =~ ^[fF]$ ]]
		then
			echo 'Data trimming using trimmomatic ...'
			currentpath=$(pwd)
			cd ${path_bbduk_software}/resources/
			java -jar ${path_trimm_software}trimmomatic-0.39.jar SE ${trimmomatic_initialization} ${pathdata}/${filename1} ${path_trimm_out}/${filename_trimmed1} ${trimming_settings_trimmomatic}
			cd ${currentpath}

		elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
		then
			echo 'Data trimming using trimmomatic ...'
			currentpath=$(pwd)
			cd ${path_bbduk_software}/resources/
			java -jar ${path_trimm_software}trimmomatic-0.39.jar PE ${trimmomatic_initialization} ${pathdata}/${filename1} ${pathdata}/${filename2} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${path_trimm_out}/${filename_trimmed2} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${trimming_settings_trimmomatic}
			cd ${currentpath}

		elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filepath2} ]]
		then
			echo 'Enter two input files for using paired end reads with Trimmomatic.'
			exit 1
		fi

	else
		echo 'Trimming software not recognized, please check settings' && exit 1
	fi
fi

# Quality report trimmed data
if [[ ${quality_check_trim} =~ ^[tT]$ ]] && [[ ${trimming}  =~ ^[tT]$ ]]
then
	echo 'Quality checking trimmed data ...'
	fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed1}
	echo 'Quality checking trimmed data completed. Results are stored at' ${path_fastqc_out}
	echo ''
	if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
	then
		echo 'Quality checking trimmed data paired end reads ...'
		fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed2}
		echo 'Quality checking trimmed data paired end reads completed. Results are stored at' ${path_fastqc_out}
		echo ''
	fi
fi


# Sequence alignment
if [[ ${paired} =~ ^[fF]$ ]]
then
	echo 'Sequence alignment ...'
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} > ${path_align_out}/${filename_sam}
elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filepath2} ]]
then
	echo 'Sequence alignment paired end interleaved ...'
	bwa mem -p ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} > ${path_align_out}/${filename_sam}
elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
then
	echo 'Sequence alignment paired end ...'
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed2} > ${path_align_out}/${filename_sam}
fi
echo 'Sequence alignment is completed. Results are stored in' ${path_align_out}/${filename_sam}
echo ''


# Creating alignment quality report
if [[ ${flagstat_report} =~ ^[tT]$ ]]
then
	echo 'Creating flagstat report sam file ...'
	samtools flagstat ${path_align_out}/${filename_sam} > ${path_align_out}/${filename1%$extension*}'_trimmed_flagstatreport.txt'
	echo 'Flagstat report sam file completed. Results are stored in' ${path_align_out}/${filename1%$extension*}'_trimmed_flagstatreport.txt'
	echo ''
fi


# Converting sam file to its binary equivalent
echo 'Converting sam to bam ...'
samtools view -b ${path_align_out}/${filename_sam} > ${path_align_out}/${filename_bam}
echo 'Converting sam to bam completed. Results are stored in' ${path_align_out}/${filename_bam}
echo ''


# Quickcheck bam file
echo 'Checking bam file ...'
samtools quickcheck ${path_align_out}/${filename_bam}
echo ''


# Indexing and sorting bam file
if [[ ${sort_and_index} =~ ^[tT]$ ]]
then
	echo 'Indexing bam file ...'
	sambamba-0.7.1-linux-static sort -m 500MB ${path_align_out}/${filename_bam}
	echo 'Indexing completed. Results are stored in' ${path_align_out}
	echo ''
fi


# Transposon mapping
if [[ ${mapping} =~ ^[tT]$ ]]
then
	echo 'Transposon mapping ...'
	cd /home/laanlab/Documents/satay/software/python_codes
	python3 ${path_python_codes}transposonmapping_satay.py ${path_align_out}/${filename_sort}
	cd ~/
	echo ''
	echo 'Transposon mapping complete. Results are stored in' ${path_align_out}
	echo ''
fi



if [[ ${delete_sam} =~ ^[tT]$ ]]
then
	echo 'Removing .sam file ...'
	rm ${path_align_out}/${filename_sam}
	echo 'sam file removed.'
fi






### Creating log file
echo ''
echo 'Creating log file ...'
echo ${filename1}	$(date +%F_%T) > ${pathdata}/${filename1%$extension*}'_log.txt'
if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filepath2} ]]
then
	echo 'Paired end reads with paired file:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo ${filename2} >> ${pathdata}/${filename1%$extension*}'_log.txt'
elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filepath2} ]]
then
	echo 'Interleaved paired end reads:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

if [[ ${trimming} =~ ^[tT]$ ]]
then
	echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo 'Trimming options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	if [[ ${trimming_software} =~ ^[bB]$ ]]
	then
		echo 'BBDuk' >> ${pathdata}/${filename1%$extension*}'_log.txt'
		echo ${trimming_settings_bbduk} >> ${pathdata}/${filename1%$extension*}'_log.txt'
	elif [[ ${trimming_software} =~ ^[tT]$ ]]
	then
		echo 'Trimmomatic' >> ${pathdata}/${filename1%$extension*}'_log.txt'
		echo ${trimmomatic_initialization} ${trimming_settings_trimmomatic} >> ${pathdata}/${filename1%$extension*}'_log.txt'
	fi
else
	echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo 'No trimming performed on reads.' >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Alignment options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
if [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filepath2} ]]
then
	echo -p ${alignment_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
else
	echo ${alignment_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Reference genome used:' ${name_refgenome} >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Adapter sequences from adapters.fa:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
cat ${path_bbduk_adapters} >> ${pathdata}/${filename1%$extension*}'_log.txt'





echo 'Processing finished.'
