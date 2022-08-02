###Using MJOLNIR to analyze eDNA metabarcoding data###
##Date: July 29, 2022
#First going to try analyzing the 6 Gully samples as a test

#load libraries
library(mjolnir)
library(ggplot2)
setwd("Raw/16S/Gully/")
#we start with demultiplexed paired end data - mjolnir can handle both de- and multiplexed data
#I first created the Gully metadata tsv file which needs 3 columns: mjolnir_agnomens, original_samples, and fastq_name_R1
#Info at https://github.com/uit-metabarcoding/MJOLNIR/tree/main/example_MJOLNIR_demultiplexed_data 

#16S forward primer: AGCGYAATCACTTGTCTYTTAA
#16S reverse primer: CRBGGTCGCCCCAACCRAA

# FREYJA will do the paired-end alignment, demultiplexing & length filtering. This will be done for each sample file separately.
lib<-"GULL"
mjolnir2_FREYJA("", cores=32, Lmin=150, Lmax=159,lib,demultiplexed=T,
                primer_F="AGCGYAATCACTTGTCTYTTAA",
                primer_R="CRBGGTCGCCCCAACCRAA",
                R1_motif="_R1.",R2_motif="_R2.")

