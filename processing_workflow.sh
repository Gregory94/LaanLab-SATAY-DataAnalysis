#!/bin/bash

#######################
#__Author__ = Gregory van Beek. LaanLab, department of Bionanoscience, Delft University of Technology
#__version__ = 1.1
#__Date last update__ = 2021-02-02
#
#Version history:
#	1.0; First working version of the workflow where the user had to open the script to change the parameters and file paths [2020-07-27]
#	1.1; Integration of a Graphical User Interface and help text (accessible via --help or -h) [2021-02-01]
#######################




if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]
then
	echo
	echo "This is a processing workflow designed for SAturated Transposon Analysis in Yeast (SATAY)."
	echo
	echo "The program can be used by running the following command: bash [path]/processing_workflow.sh, where [path] is the path to the processing_workflow.sh file."
	echo
	echo "The program can trim sequencing reads, create quality reports of raw data and of the trimmed data, align the reads to a reference genome (S288C yeast genome, downloaded from the SGD) in .sam- and .bam-format, sort and index the bam file and perform transposon-mapping."
	echo "The following tools are used:"
	echo "- quality report: FASTQC"
	echo "- trimming: BBDuk or Trimmomatic"
	echo "- alignment: BWA MEM"
	echo "- create flagstat report after alignment: SAMTools"
	echo "- sort and index .bam-files: SAMBamba"
	echo "- transposon-mapping: Python3 together with custom python software (https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/dev_Gregory/Python_TransposonMapping/transposonmapping_satay.py)"
	echo
	echo "The program does not take inputs from the command line, except for the help text which can be accessed with the arguments --help or -h. It can handle both single-end data, paired-end data and paired-end interleaved data."
	echo
	echo "When the program is started, a window appears where the datafile(s) can be selected. The datafiles should be in fastq format and must have the extension .fastq or .fq. They can be either unpacked or gzipped (i.e. having the extension .gz)."
	echo "Select the right extension in the bottom right corner and navigate to the datafile(s)."
	echo "In case of single-end data or paired-end interleaved data select only one file. In case of paired-end data where the pairs are stored in two seperate files, select two files by holding ctrl-button and clicking the two files."
	echo
	echo "After pressing 'ok', a new window appears where some options and parameters can be selected."
	echo "- 'file primary reads' and 'file secondary reads': Show the selected file(s). If only one file is chosen, the 'file secondary reads' will show 'none'."
	echo "- 'Data type': Select whether the reads are paired-end or single-end. If two data files were chosen but this setting is set to 'Single-end', the secondary reads file will be ignored."
	echo "- 'which trimming to use': Select whether to use bbduk or trimmomatic for trimming the reads or select 'Do_not_trim' to prevent trimming of the reads."
	echo "- 'trimming settings': Input trimming settings. See the documentation of the selected trimming software which settings can be applied. When 'which trimming to use' is set to 'Do_not_trim', this field will be ignored. Sequences that need to trimmed (e.g. adapter or primer sequences) have to be entered in the adapters file which can be accessed using the 'Open adapters file' button on the bottom of the window. NOTE 1: For bbduk do not input 'interleaved=t' when using interleaved data. For trimmomatic do not input 'SE' or 'PE' to indicate single-end or paired-end data. This will all be automatically set depending on your selection in the 'Data type'-field. NOTE 2: For trimmomatic, when using ILLUMINACLIP, do not specify the path to the adapters file as this is inserted automatically (see example settings below)."
	echo "EXAMPLE SETTINGS:"
	echo "	- Trimmmomatic: ILLUMINACLIPPING:1:30:10 TRAILING:20 SLIDINGWINDOW:5:10 MINLEN:15 [http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf]"
	echo "	- bbduk: ktrim=l k=15 mink=10 hdist=1 qtrim=r trimq=10 minlen=30 [https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/]"
	echo "- 'alignment settings': Input alignment settings. See the documentation of BWA MEM which settings can be applied. NOTE: Do not set -p for smart pairing (i.e. interleaved paired-end data). This will be automatically set depending on your selection in the 'Data type'-field. [http://bio-bwa.sourceforge.net/bwa.shtml]"
	echo "- 'Quality checking raw data': Perform a fastqc quality check on the raw reads."
	echo "- 'Quality checking trimmed data': Perform a fastqc quality check on the trimmed reads. This setting is ignored if 'which trimming to use' it set to 'Do_not_trim'."
	echo "- 'Quality check interrupt': This quits the program after performing the quality report on the raw dataset and creating a temporary file with your settings. This allows you to check the quality report before continuing. To continue the program, restart the program (using bash processing_workflow.sh). It will automatically set the options you have chosen the first time, but these can be changed if this is necessary depending on the outcome of the quality report. This option can be useful if you have no idea how the dataset looks. This requires 'Quality checking raw data'."
	echo "- 'Delete sam file': After alignment the .sam file is converted to its binary equivalent and only this .bam file is used for downstream processing. Since the .sam file typically requires a lot of memory, this is can be deleted. It is recommended to keep the .sam file only for manual checking te alignment results."
	echo "- 'Sort and index bam files.': This is needed for transposon-mapping and for many other downstream processes. It is recommended to always leave this on."
	echo "- 'Transposon mapping': This custom python script requires sorting and indexing of the .bam file and creates the following files:"
	echo "	- .bed file: Creates list of all insertion locations with the number of reads in each location in bed-format."
	echo "	- .wig file: Creates list of all insertion locations with the number of reads in each location in wig-format. Small difference with the bed-file is that here reads from insertions at the same location but with different orientation are added up. In the bed-file these are regarded two separate insertions."
	echo "	- 4 .txt-files: List of all genes with the number of insertions and reads in each gene. The files are different in whether they show all genes or all annotated essential genes and whether they also show the distribution of insertions within the genes."
	echo "- 'Create flagstat report': Creates a flagstat report based on the .bam file."
	echo "- 'Open adapters file': Opens the text file where the adapter and primer sequences can be entered that will be trimmed. Enter the sequences in fasta format."
	echo
	echo "Questions, recommendations and issues can be noted at https://www.github.com/Gregory94/LaanLab-SATAY-DataAnalysis/issues/33"
	echo
	echo
	echo "Dependencies:"
	echo "- Zenity"
	echo "- YAD"
	echo "- fastqc"
	echo "- bbmap"
	echo "- trimmomatic"
	echo "- bwa"
	echo "- samtools"
	echo "- sambamba"
	echo "- python3"
	echo "	- transposonmapping_satay.py [https://www.github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Python_TransposonMapping/transposonmapping_satay.py]"
	echo "	- numpy"
	echo "	- pysam"
	echo "	- python_modules [https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Python_TransposonMapping/python_modules]"
	echo "	- data_files [https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Data_Files]"
	exit 0
fi


################### Define paths ###########################

#CACHEFILE (this is a temporary file that is created to store user settings when 'Quality check interrupt is set to true).
#cachefile="/home/gregoryvanbeek/Desktop/processing_workflow_cache.txt"
cachefile="/home/laanlab/Documents/satay/software/processing_workflow_cache.txt"

#ADAPTERFILE (this refers to the file with adapter sequences that are used for trimming).
#adapterfile="/home/gregoryvanbeek/Documents/Software/BBMap/bbmap/resources/adapters.fa"
adapterfile="/home/laanlab/Documents/satay/software/bbmap/resources/adapters.fa"

#REFERENCE GENOME (path to the fasta file of the reference genome).
#path_refgenome='/home/gregoryvanbeek/Documents/Reference_Sequences/Reference_Sequence_S288C/S288C_reference_sequence_R64-2-1_20150113.fsa'
path_refgenome=/home/laanlab/Documents/satay/reference_sequences/Cerevisiae_S288C_reference/S288C_reference_sequence_R64-2-1_20150113.fsa

#DDBUK SOFTWARE (path to bbduk for trimming).
#path_bbduk_software=/home/gregoryvanbeek/Documents/Software/BBMap/bbmap/
path_bbduk_software=/home/laanlab/Documents/satay/software/bbmap/

#TRIMMOMATIC (path to trimmomatic for trimming).
#path_trimm_software=/home/gregoryvanbeek/Documents/Software/Trimmomatic-0.39/
path_trimm_software=/home/laanlab/Documents/satay/software/Trimmomatic-0.39/

#PYTHON CODES (path to python code for transposon_mapping).
#path_python_codes=/home/gregoryvanbeek/Documents/Software/python_codes/
path_python_codes=/home/laanlab/Documents/satay/software/python_codes/

############################################################




if [ ! -f $cachefile ];
then
	fileselections=`yad --width=1000 --height=400 --title="Select fastq file" --center --on-top --buttons-layout=spread --multiple --file-selection="Please select datafile (or two files in case of paired-end noninterleaved fastq files)" --file-filter="*.fq" --file-filter="*.fastq" --file-filter="*.fq.gz" --file-filter="*.fastq.gz"`
	filepath1=$(echo $fileselections | awk 'BEGIN {FS="|" } { print $1 }')
	filepath2=$(echo $fileselections | awk 'BEGIN {FS="|" } { print $2 }')
	
	#if no file selected, set filepath to none. If no filepath1, exit.
	[ -z $filepath1 ] && filepath1='none' && exit 0
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
	"Single-end!Paired-end" \
	"bbduk!trimmomatic!Do_not_trim" \
	"ktrim=l k=15 mink=10 hdist=1 qtrim=r trimq=10 minlen=30" \
	" -v 2" \
	"False" \
	"TRUE" \
	"FALSE" \
	"TRUE" \
	"TRUE" \
	"TRUE" \
	"TRUE" \
	"bash -c 'xdg-open ${adapterfile}'"`

	if [ ! -z "$settings" ] && [ $filepath1 != "none" ] && [ $(echo $settings | awk 'BEGIN {FS="|" } { print $9 }') == TRUE ] && [ $(echo $settings | awk 'BEGIN {FS="|" } { print $7 }') == TRUE ] #Create cachefile only if settings or filepath1 is not empty and Qualitycheck interrupt is set to True and Quality check raw files is set to True.
	then
		echo $settings >> $cachefile
		echo 'Cache file created.'
	fi

elif [ -f $cachefile ];
then

	previoussettings=`head -n 1 $cachefile`
	echo $previoussettings

	settings=`yad --width=1100 --height=500 --title="Processing settings" --text="Settings" --center --on-top --buttons-layout=spread --form \
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
	"bash -c 'xdg-open ${adapterfile}'"`

	filepath1=$(echo $settings | awk 'BEGIN {FS="|" } { print $1 }')
	filepath2=$(echo $settings | awk 'BEGIN {FS="|" } { print $2 }')

	rm ${cachefile}
	if [[ -z ${settings} ]]
	then
		echo 'Process canceled.'
		exit 1
	fi
fi


############################################################

# Define whether data is paired-end
paired=$(echo $settings | awk 'BEGIN {FS="|" } { print $3 }')
echo 'paired ' $paired

# Define which trimming tool to use
trimming_software=$(echo $settings | awk 'BEGIN {FS="|" } { print $4 }')
echo 'trimming_software ' $trimming_software

###    bbduk    ###
trimming_settings=$(echo $settings | awk 'BEGIN {FS="|" } { print $5 }')
echo 'trimming_settings bbduk ' $trimming_settings

### trimmomatic ###
trimmomatic_initialization='-phred33'
trimming_settings=$(echo $settings | awk 'BEGIN {FS="|" } { print $5 }')
echo 'trimming_settings trimmomatic ' $trimming_settings

# Set options for alignment software (bwa mem)
alignment_settings=$(echo $settings | awk 'BEGIN {FS="|" } { print $6 }')
echo 'alignment_settings ' $alignment_settings

# Create sorted and indexed bam file ('TRUE' or 'FALSE')
sort_and_index=$(echo $settings | awk 'BEGIN {FS="|" } { print $11 }')
echo 'sort_and_index ' $sort_and_index

# Apply transposon mapping ('TRUE' or 'FALSE') (requires sort_and_index=TRUE)
mapping=$(echo $settings | awk 'BEGIN {FS="|" } { print $12 }')
echo 'mapping ' $mapping

# Delete sam file ('TRUE' or 'FALSE')? This file is always converted to its binary equivalent (.bam ) and the sam file is rarely used but takes up relatively a lot of memory.
delete_sam=$(echo $settings | awk 'BEGIN {FS="|" } { print $10 }')
echo 'delete_sam ' $delete_sam

# Create a quality report of the alignment based on the sam file ('TRUE' or 'FALSE')
flagstat_report=$(echo $settings | awk 'BEGIN {FS="|" } { print $13 }')
echo 'flagstat_report ' $flagstat_report

# Create quality report of raw data (before trimming) ('TRUE' or 'FALSE')?
quality_check_raw=$(echo $settings | awk 'BEGIN {FS="|" } { print $7 }')
echo 'quality_check_raw ' $quality_check_raw

# Create quality report of trimmed data (after trimming) ('TRUE' or 'FALSE')?
quality_check_trim=$(echo $settings | awk 'BEGIN {FS="|" } { print $8 }')
echo 'quality_check_trim ' $quality_check_trim

# Determine whether the script should automatically continue after creating the first quality report. Set to True if you might want to make changes depending on the quality report of the raw data ('TRUE' or 'FALSE').
qualitycheck_interrupt=$(echo $settings | awk 'BEGIN {FS="|" } { print $9 }')
echo 'qualitycheck_interrupt ' $qualitycheck_interrupt

############################################################




echo 'Preparing processing for' $(basename ${filepath1}) '...'
echo ''


# Check paths selected datafiles.
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
	echo ''
fi

if [[ -z "${settings}" ]]
then
	echo 'Process canceled.' && exit 1
fi



# Get extension of the file
extension='.'$(echo $filename1 | rev | cut -d. -f1 | rev)
if [[ ${extension} == '.gz' ]]
then
	extension='.'$(echo $filename1 | rev | cut -d. -f2 | rev)
fi


# Define filename for trimming and alignment results
if ! [[ ${trimming_software} == 'Do_not_trim' ]] && [[ ${trimming_settings} == '' ]]
then
	echo 'trimming settings:' ${trimming_settings}
	echo 'WARNING: No trimming setting were specified. Trimming set to Do_not_trim'
	trimming_software='Do_not_trim'
fi

if ! [[ ${trimming_software} =~ 'Do_not_trim' ]]
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


# Define path output directory trimming
if ! [[ ${trimming_software} == 'Do_not_trim' ]]
then
	path_trimm_out=${pathdata}/trimm_out
	[ ! -d ${path_trimm_out} ] && echo 'Creating trimming output folder ...' && mkdir ${path_trimm_out} || echo 'Folder for trimming output exists with name:' $(basename ${path_trimm_out})
else
	path_trimm_out=${pathdata}
fi

# Define path output directory alignment
path_align_out=${pathdata}/align_out
[ ! -d ${path_align_out} ] && echo 'Creating alignment output folder ...' && mkdir ${path_align_out} || echo 'Folder for alignment output exists with name:' $(basename ${path_align_out})

# Define paths to reference genomes S288C
name_refgenome='S288C'
if [ ! -f ${path_refgenome} ] #if path to reference genome does not exist
then
	echo 'ERROR: Reference genome not found at location:' ${path_refgenome} && exit 1
else
	echo 'Reference genome:' ${name_refgenome}
fi

# Define path bbduk software
[ ! -d ${path_ddbuk_software} ] && echo 'WARNING: Path to bbduk software does not exists.'

# Define path trimmomatic software
[ ! -d ${path_trimm_software} ] && echo 'WARNING: Path to trimmomatic software does not exists.'

# Define path to python script
[ ! -d ${path_python_codes} ] && echo 'WARNING: Path to python codes does not exists.'




### Start processing workflow
echo ''
echo 'Start processing ...'
echo ''


# Quality checking raw data
if [[ ${quality_check_raw} == TRUE ]]
then
	if [[ ! -e ${path_fastqc_out}/${filename1%$extension*}'_fastqc.html' ]]
	then
		echo 'Quality checking raw data ...'
		fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename1}
		echo 'Quality checking raw data completed. Results are stored at' ${path_fastqc_out}
		echo ''

		if [[ ${paired} =~ 'Paired-end' ]] && ! [[ ${filepath2} =~ 'none' ]]
		then
			fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename2}
			echo 'Quality checking raw data paired end reads completed. Results are stored at' ${path_fastqc_out}
			echo ''
		fi
	else
		echo 'Quality report raw data already exists. Skipping fastqc'
	fi

	if [[ ${qualitycheck_interrupt} == TRUE ]]
	then
		read -p 'Continue processing? (press "y" if yes, press "n" if no): ' -n 1 -r
		echo
		if [[ ! $REPLY =~ ^[yY]$ ]]
		then
			exit 0
		else
			rm ${cachefile}
		fi
	fi
fi


# Trimming
if ! [[ ${trimming_software} == 'Do_not_trim' ]]
then
	if [[ ${trimming_software} == 'bbduk' ]]
	then
		if [[ ${paired} == 'Single-end' ]]
		then
			echo 'Data trimming using bbduk single end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${adapterfile} ${trimming_settings}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
			echo ''
		elif [[ ${paired} == 'Paired-end' ]] && ! [[ ${filepath2} == 'none' ]]
		then
			echo 'Data trimming using bbduk paired end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g in1=${pathdata}/${filename1} out1=${path_trimm_out}/${filename_trimmed1} in2=${pathdata}/${filename2} out2=${path_trimm_out}/${filename_trimmed2} ref=${adapterfile} ${trimming_settings}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1} 'and for the paired end reads in' ${path_trimm_out}/${filename_trimmed2}
			echo ''

		elif [[ ${paired} == 'Paired-end' ]] && [[ ${filepath2} == 'none' ]]
		then
			echo 'Data trimming using bbduk paired end reads...'
			${path_bbduk_software}bbduk.sh -Xmx2g interleaved=t in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${adapterfile} ${trimming_settings}
			echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
			echo ''
		fi

	elif [[ ${trimming_software} == 'trimmomatic' ]]
	then
		if [[ ${trimming_settings} == *'ILLUMINACLIP'* ]]
		then
			clip_startlocation=$(echo ${trimming_settings} | grep -b -o ILLUMINACLIP | awk 'BEGIN {FS=":"}{print $1}')
			clip_endlocation=$(echo ${clip_startlocation}+12 | bc)
			trimming_settings_trimmomatic=$(echo ${trimming_settings:0:${clip_startlocation}}${trimming_settings:${clip_startlocation}:12}:${adapterfile}${trimming_settings:${clip_endlocation}})
		else
			trimming_settings_trimmomatic=${trimming_settings}
		fi

		echo 'Trimming settings trimmomatic: '${trimming_settings_trimmomatic}

		if [[ ${paired} == 'Single-end' ]]
		then
			echo 'Data trimming using trimmomatic ...'
			currentpath=$(pwd)
			cd ${path_bbduk_software}/resources/
			java -jar ${path_trimm_software}trimmomatic-0.39.jar SE ${trimmomatic_initialization} ${pathdata}/${filename1} ${path_trimm_out}/${filename_trimmed1} ${trimming_settings_trimmomatic}
			cd ${currentpath}

		elif [[ ${paired} == 'Paired-end' ]] && ! [[ ${filepath2} == 'none' ]]
		then
			echo 'Data trimming using trimmomatic ...'
			currentpath=$(pwd)
			cd ${path_bbduk_software}/resources/
			java -jar ${path_trimm_software}trimmomatic-0.39.jar PE ${trimmomatic_initialization} ${pathdata}/${filename1} ${pathdata}/${filename2} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${path_trimm_out}/${filename_trimmed2} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${trimming_settings_trimmomatic}
			cd ${currentpath}

		elif [[ ${paired} = 'Paired-end' ]] && [[ ${filepath2} == 'none' ]]
		then
			echo 'Enter two input files for using paired end reads with Trimmomatic.'
			exit 1
		fi

	else
		echo 'Trimming software not recognized, please check settings' && exit 1
	fi
fi

# Quality report trimmed data
if [[ ${quality_check_trim} == TRUE ]] && ! [[ ${trimming_software} == 'Do_not_trim' ]]
then
	echo 'Quality checking trimmed data ...'
	fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed1}
	echo 'Quality checking trimmed data completed. Results are stored at' ${path_fastqc_out}
	echo ''
	if [[ ${paired} == 'Paired-end' ]] && ! [[ ${filepath2} == 'none' ]]
	then
		echo 'Quality checking trimmed data paired end reads ...'
		fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed2}
		echo 'Quality checking trimmed data paired end reads completed. Results are stored at' ${path_fastqc_out}
		echo ''
	fi
fi



# Sequence alignment
if [[ ${paired} == 'Single-end' ]]
then
	echo 'Sequence alignment ...'
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} > ${path_align_out}/${filename_sam}
elif [[ ${paired} == 'Paired-end' ]] && [[ ${filepath2} == 'none' ]]
then
	echo 'Sequence alignment paired end interleaved ...'
	bwa mem -p ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} > ${path_align_out}/${filename_sam}
elif [[ ${paired} == 'Paired-end' ]] && ! [[ ${filepath2} == 'none' ]]
then
	echo 'Sequence alignment paired end ...'
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed2} > ${path_align_out}/${filename_sam}
fi
echo 'Sequence alignment is completed. Results are stored in' ${path_align_out}/${filename_sam}
echo ''



# Creating alignment quality report
if [[ ${flagstat_report} == TRUE ]]
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
if [[ ${sort_and_index} == TRUE ]]
then
	echo 'Indexing bam file ...'
	sambamba-0.7.1-linux-static sort -m 500MB ${path_align_out}/${filename_bam}
	echo 'Indexing completed. Results are stored in' ${path_align_out}
	echo ''
fi


# Transposon mapping
if [[ ${mapping} == TRUE ]]
then
	echo 'Transposon mapping ...'
	cd ${path_python_codes}
	python3 ${path_python_codes}transposonmapping_satay.py ${path_align_out}/${filename_sort}
	cd ~/
	echo ''
	echo 'Transposon mapping complete. Results are stored in' ${path_align_out}
	echo ''
fi



if [[ ${delete_sam} == TRUE ]]
then
	echo 'Removing .sam file ...'
	rm ${path_align_out}/${filename_sam}
	echo 'sam file removed.'
fi




### Creating log file
echo ''
echo 'Creating log file ...'
echo ${filename1}	$(date +%F_%T) > ${pathdata}/${filename1%$extension*}'_log.txt'
if [[ ${paired} == 'Paired-end' ]] && ! [[ ${filepath2} == 'none' ]]
then
	echo 'Paired end reads with paired file:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo ${filename2} >> ${pathdata}/${filename1%$extension*}'_log.txt'
elif [[ ${paired} == 'Paired-end' ]] && [[ ${filepath2} == 'none' ]]
then
	echo 'Interleaved paired end reads:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

if ! [[ ${trimming_software} == 'Do_not_trim' ]]
then
	echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo 'Trimming options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	if [[ ${trimming_software} == 'bbduk' ]]
	then
		echo 'BBDuk' >> ${pathdata}/${filename1%$extension*}'_log.txt'
		echo ${trimming_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
	elif [[ ${trimming_software} == 'trimmomatic' ]]
	then
		echo 'Trimmomatic' >> ${pathdata}/${filename1%$extension*}'_log.txt'
		echo ${trimmomatic_initialization} ${trimming_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
	fi
else
	echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo 'No trimming performed on reads.' >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Alignment options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
if [[ ${paired} == 'Paired-end' ]] && [[ ${filepath2} == 'none' ]]
then
	echo ${alignment_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
else
	echo ${alignment_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Reference genome used:' ${name_refgenome} >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Adapter sequences from adapters.fa:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
cat ${adapterfile} >> ${pathdata}/${filename1%$extension*}'_log.txt'





echo 'Processing finished.'
