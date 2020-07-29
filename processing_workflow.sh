#!/bin/bash

### This workflow is developed for automatically preprocessing sequencing data.
### The variables stored in the USER SETTINGS block should be checked and potentially changed for each dataset.
### The other lines can be left as they are.
### A folder with name given in the 'foldername' variable is created in '~/Documents/data_processing'.
### Within this folder, three more folders are generated for the output of the quality report, trimming and alignment.
### Settings for trimming and alignment should be set in the 'trimming_settings' and 'alignment_settings' variables respectively.
### Adapter sequences for trimming should be stored in '~/Documents/Software/BBMap/bbmap/resources/adapters.fa' respecting fasta convention.
### The software is called in this order: 1.Quality checking (fastqc) 2.trimming (bbduk) 3.Quality checking trimmed data (fastqc) 4.alignment (bwa) 5.converting sam to bam (samtools) 6.indexing bam file (sambamba) 7.transposon mapping (python).
### Finally, the data file and all the created files are moved to the shared folder.
### There are two trimming software packages, 'bbduk' and 'trimmomatic'.
### You can use either one of them by settings the appropriate option in the USER SETTINGS.
### When selecting one trimming software, the options of the other program will be ignored.
### The workflow will ask after the first quality report if you want to continue.
### Pressing 'n' will abort the workflow and allows you to make changes according to the quality report.
### Restarting the workflow cause it to skip over the initial quality report (unless you deleted it) and continues with the trimming and alignment steps.
### In case of emergency, press ctrl-c (possibly multiple times) in the terminal.





####################### USER SETTINGS ######################

# Define filename (can also be a zipped file ending with .gz)
filename='SRR062634.filt.fastq.gz'

# Define foldername where the analysis results are stored
foldername='SRR062634_test'


###### Set options for trimming software ######
# set 'b' for bbduk, set 't' for trimmomatic
trimming_software='b'

###    bbduk    ###
trimming_settings_bbduk='ktrim=r k=4 hdist=0 qtrim=r trimq=4'
## Set adapter sequences
## Open file using xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
###################

### trimmomatic ###
trimmomatic_initialization='SE -phred33'
trimming_settings_trimmomatic='ILLUMINACLIP:adapters.fa:0:30:10 SLIDINGWINDOW:10:4 MINLEN:30'
## Set adapter sequences
## Open file using xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
###################
###############################################


# Set options for alignment software (bwa mem)
alignment_settings='-B 2 -O 3'

# Set which reference genome you want to use (type 's' for S288C or 'w' for W303)
refgenome='s'


# Create sorted and indexed bam file ('y' for yes, 'n' for no)?
sort_and_index='y'


# Apply transposon mapping
mapping='y'


# Save sam file ('y' for yes, 'n' for no)? This file is always converted to its binary equivalent (.bam ) and the sam file is rarely used but takes up relatively a lot of memory.
delete_sam='n'

############################################################





# Ask for confirmation to continue after quality report raw data (t for True or f for False).
# When False, the program continues automatically.
ask_user=T




echo 'Preparing processing for' ${filename} '...'
echo ''

# Define filename for trimming and alignment results
filename_trimmed=${filename%.fastq*}'_trimmed.fastq'
filename_sam=${filename%.fastq*}'_trimmed.sam'
filename_bam=${filename%.fastq*}'_trimmed.bam'
filename_sort=${filename%.fastq*}'_trimmed.sorted.bam'

# Define full path to data folder and create it if it doesn't exists
pathdata=~/Documents/data_processing/${foldername}
[ ! -d ${pathdata} ] && echo 'Creating datafolder ...' && mkdir ${pathdata}

# Define path output directory fastqc
path_fastqc_out=${pathdata}/fastqc_out
[ ! -d ${path_fastqc_out} ] && echo 'Creating fastqc output folder ...' && mkdir ${path_fastqc_out} || echo 'Folder for fastqc output exists with name:' $(basename ${path_fastqc_out})

# Define path output directory trimming
path_trimm_out=${pathdata}/trimm_out
[ ! -d ${path_trimm_out} ] && echo 'Creating trimming output folder ...' && mkdir ${path_trimm_out} || echo 'Folder for trimming output exists with name:' $(basename ${path_trimm_out})

# Define path output directory alignment
path_align_out=${pathdata}/align_out
[ ! -d ${path_align_out} ] && echo 'Creating alignment output folder ...' && mkdir ${path_align_out} || echo 'Folder for alignment output exists with name:' $(basename ${path_align_out})

# Define path shared folder
path_sf=/media/sf_VMSharedFolder_Ubuntu64_1/

# Define paths to reference genomes (both S288C and W303)
if [[ ${refgenome} =~ ^[sS]$ ]]
then
	path_refgenome=/home/gregoryvanbeek/Documents/Reference_Sequences/Reference_Sequence_S288C/S288C_reference_sequence_R64-2-1_20150113.fsa
	name_refgenome='S288C'
	echo 'Reference genome:' ${name_refgenome}
elif [[ ${refgenome} =~ ^[wW]$ ]]
then
	path_refgenome=/home/gregoryvanbeek/Documents/Reference_Sequences/Reference_Sequence_W303/W303_SGD_2015_JRIU00000000.fsa
	name_refgenome='W303'
	echo 'Reference genome:' ${name_refgenome}
else
	echo 'ERROR: Reference genome not defined. Please check settings.' && exit 1
fi

# Define path bbduk software
path_bbduk_software=~/Documents/Software/BBMap/bbmap/
path_bbduk_adapters=${path_bbduk_software}/resources/adapters.fa
[ ! -d ${path_ddbuk_software} ] && echo 'WARNING: Path to bbduk software does not exists.'

# Define path trimmomatic software
path_trimm_software=~/Documents/Software/Trimmomatic-0.39/
[ ! -d ${path_trimm_software} ] && echo 'WARNING: Path to trimmomatic software does not exists.'

# Define path to python script
path_python_codes=~/Documents/Software/python_codes/
[ ! -d ${path_python_codes} ] && echo 'WARNING: Path to python codes does not exists.'



# Check if datafile is already in the datafolder. If not, move it to the datafolder
[ -e ${path_sf}${filename} ] && echo 'Moving' ${filename} 'to' ${foldername} '...' && mv ${path_sf}${filename} ${pathdata} && echo 'Moving complete.' && sleep 1s

[ ! -e ${pathdata}/${filename} ] && echo 'ERROR:' ${filename} 'does not exists in' $(basename ${pathdata}) '. Cannot proceeed with processing.' && exit 1





### Creating log file
echo ''
echo 'Creating log file ...'
echo ${filename}	$(date +%F_%T) > ${pathdata}/${filename%.fastq*}'_log.txt'
echo '' >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo 'Trimming options:' >> ${pathdata}/${filename%.fastq*}'_log.txt'

if [[ ${trimming_software} =~ ^[bB]$ ]]
then
	echo 'BBDuk' >> ${pathdata}/${filename%.fastq*}'_log.txt'
	echo ${trimming_settings_bbduk} >> ${pathdata}/${filename%.fastq*}'_log.txt'
elif [[ ${trimming_software} =~ ^[tT]$ ]]
then
	echo 'Trimmomatic' >> ${pathdata}/${filename%.fastq*}'_log.txt'
	echo ${trimmomatic_initialization} ${trimming_settings_trimmomatic} >> ${pathdata}/${filename%.fastq*}'_log.txt'
fi

echo '' >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo 'Alignment options:' >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo ${alignment_settings} >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo '' >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo 'Reference genome used:' ${name_refgenome} >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo '' >> ${pathdata}/${filename%.fastq*}'_log.txt'
echo 'Adapter sequences from adapters.fa:' >> ${pathdata}/${filename%.fastq*}'_log.txt'
cat ${path_bbduk_adapters} >> ${pathdata}/${filename%.fastq*}'_log.txt'





### Start processing workflow
echo ''
echo 'Start processing ...'
echo ''


# Quality checking raw data
if [[ ! -e ${path_fastqc_out}/${filename%.fastq*}'_fastqc.html' ]]
then
	echo 'Quality checking raw data ...'
	fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename}
	echo 'Quality checking raw data completed. Results are stored at' ${path_fastqc_out}
	echo ''
else
	echo 'Quality report raw data already exists. Skipping fastqc'
fi


if [[ ${ask_user} =~ ^[tT]$ ]]
then
	read -p 'Continue processing? (press "y" if yes, press "n" if no): ' -n 1 -r
	echo
	if [[ ! $REPLY =~ ^[yY]$ ]]
	then
		exit 1
	fi
fi

# Trimming
if [[ ${trimming_software} =~ ^[bB]$ ]]
then
	echo 'Data trimming using bbduk ...'
	${path_bbduk_software}bbduk.sh -Xmx1g in=${pathdata}/${filename} out=${path_trimm_out}/${filename_trimmed} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
	echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed}
	echo ''
elif [[ ${trimming_software} =~ ^[tT]$ ]]
then
	echo 'Data trimming using trimmomatic ...'
	currentpath=$(pwd)
	cd ${path_bbduk_software}/resources/
	java -jar ${path_trimm_software}trimmomatic-0.39.jar ${trimmomatic_initialization} ${pathdata}/${filename} ${path_trimm_out}/${filename_trimmed} ${trimming_settings_trimmomatic}
	cd ${currentpath}
else
	echo 'Trimming software not recognized, please check settings'
	exit 1
fi

# Quality report trimmed data
echo 'Quality checking trimmed data ...'
fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed}
echo 'Quality checking trimmed data completed. Results are stored at' ${path_fastqc_out}
echo ''

# Sequence alignment
echo 'Sequence alignment ...'
bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed} > ${path_align_out}/${filename_sam}
echo 'Sequence alignment is completed. Results are stored in' ${path_align_out}/${filename_sam}
echo ''

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
if [[ ${sort_and_index} =~ ^[yY]$ ]]
then
	echo 'Indexing bam file ...'
	sambamba-0.7.1-linux-static sort -m 500MB ${path_align_out}/${filename_bam}
	echo 'Indexing completed. Results are stored in' ${path_align_out}
	echo ''
fi


# Transposon mapping
if [[ ${mapping} =~ ^[yY]$ ]]
then
	echo 'Transposon mapping ...'
	cd ~/Documents/Software/python_codes
	python3 ${path_python_codes}transposonmapping_satay.py ${path_align_out}/${filename_sort}
	cd ~/
	echo ''
	echo 'Transposon mapping complete. Results are stored in' ${path_align_out}
	echo ''
fi



# Moving results to shared folder.
echo 'Processing completed.'

if [[ ${delete_sam} =~ ^[yY]$ ]]
then
	echo 'Removing .sam file ...'
	rm ${path_align_out}/${filename_sam}
	echo 'sam file removed.'
fi

echo 'Moving results to shared folder ...'
mv ${pathdata} ${path_sf}
[ -d ${path_sf}$(basename ${pathdata}) ] && echo 'Files sucessfully moved to shared folder.' || 'WARNING: Files not moved to shared folder.'





