#!/bin/bash
#run a FastQC loop for all gzipped fastq's

cd /hdd5/eDNA_Data/Raw/12S #change folder here per marker
#unzip if you want to
# for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done

mkdir fastqc_output    
fastqc -o fastqc_output/ *.fastq.gz    #first -o part is the output then the input files
multiqc -o fastqc_output/ fastqc_output    
