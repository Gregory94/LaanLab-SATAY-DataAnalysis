Welcome to the Ubuntu Virtual Machine that is created for processing of SATAY data.
This file explains the layout of the virtual machine (VM) and how to perform the processing of raw sequencing data, specifically for SATAY data.
All the required software tools are already preinstalled and ready to use.

If the VM is setup according to the installation guide (https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/docs/Installation_Guide_SATAY_Analysis_Software.md), there should appear a shared folder on the Desktop with the name: 'sf_VMSharedFolder_Ubuntu64_1'.
This shared folder can be used for easy sharing of files between the VM and the host system (i.e. Windows).
In the host system the files are located in the folder that was selected during the setup of the VM.

The structure of the most important folders and files in the VM is shown next, where each indent indicates a subfolder and names that are defined by the user are indicated between ${}:
~/Documents
	data_processing
		${datafolder}
			
			fastqc_out
			trimm_out
			align_out

The main operations are performed in the terminal app (located in the left bar on the screen).
When you open this, the default location is ~/.


The structure of the most important folders and files in the VM is shown next, where each arrow indicates a subfolder (note that the subfolders in the data_processing folder are created when the workflows is started):

~/Documents
*processing_workflow.sh*

> data_processing

> > datafolder

> > >*datafile.fastq*

> > > fastqc_out

> > > trimm_out

> > > align_out

> Reference_sequences

> > Reference_sequence_S288C

> > >*S288C_reference_sequence_R64-2-1_20190113.fsa*

> > Reference_sequence_W303

> > >*W303_SGD_2015_JRIU00000000.fsa*

> Software

> > BBMap

> > > bbmap

> > > > resouces

> > > > >*adapters.fa*

> > BWA

> > cmpfastq

> > FastQC

> > java_download

> > Sambamba

> > SAMTools_bcftools

> > Trimmomatic_0.39