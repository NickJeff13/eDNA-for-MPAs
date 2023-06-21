#!/bin/bash
#run a FastQC loop for all gzipped fastq's

cd ~/eDNA/Musquash/Data/Raw/Musquash-18S #change folder here per marker

#unzip if you want to
# for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done

mkdir fastqc_output    
fastqc -t 60 -o fastqc_output/ *.fastq.gz    #first -o part is the output then the input files, -t sets number of cores
multiqc -o fastqc_output/ fastqc_output    
