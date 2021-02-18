#!/bin/bash

#inputfile1="/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
#inputfile2="/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"
inputfile1="C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
inputfile2="C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"
sample_per_reads_file="C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\reads_dataframe.csv"



if [ ! -f $inputfile1 ]
then
	echo 'ERROR:' $inputfile1 'not found.'
fi
if [ ! -f $inputfile2 ]
then
	echo 'ERROR:' $inputfile2 'not found.'
fi



echo "Input files:" $inputfile1 $inputfile2

extension='.'$(echo $inputfile1 | rev | cut -d. -f1 | rev)

outputfile_sample1=${inputfile1%$extension*}'_sample1_interleaved.fq'
outputfile_sample2=${inputfile1%$extension*}'_sample2_interleaved.fq'



echo "Demultiplexing input files ..."

awk -F" " '/^@/ && FNR==NR {lines[$1]; next} $1 in lines {print $2}' $inputfile1 $inputfile2

# awk -F" " '/^@/ && FNR==NR lines[$1]; next {
# if($1 in lines)
# {
# {print $1}
# }' $inputfile1 $inputfile2







