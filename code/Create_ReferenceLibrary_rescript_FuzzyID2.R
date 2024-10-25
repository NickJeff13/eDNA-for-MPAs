#use rescript_createReferenceDB.sh to create the preliminary library, which includes dereplicating and trimming for sequence length

#Follow https://github.com/dfo-mar-mpas/can_marinefish_ref from Step 5 to polish the DB

library(haplotypes)
#setwd
setwd("~/eDNA/Musquash/DBs/12S")

#load data
refFish <- read.fas("12S_reference_library_Actinopterygii_FuzzyID2.fasta") #1474 seqs

data <- read.fas("fish-12S-ref-seqs-FINAL_FuzzyID2.fasta") #46063 seqs
test <- data[1:400,as.matrix=F]


#infer haplotypes
haplos <- haplotype(data, indels="s") #Code gaps using simple indel coding method for the distance calculations.
