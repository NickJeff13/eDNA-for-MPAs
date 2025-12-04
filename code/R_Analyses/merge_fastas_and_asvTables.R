#Filter fasta files to only include sequences from the filtered ASV files

library(dplyr)
library(tidyverse)


# 2024 Perley -------------------------------------------------------------
# This includes some AZMP samples from the Gully and Fundian Channel

## 12S
#Redoing this filtration step from filters_ASVs.R but keeping OTU.ID
sab24.12s.filt <- filter_low_reads(sab24.12s.merge %>% 
                                     drop_na() %>%
                                     filter(V3>97.99 & V7 %in% c("bony fishes","whales & dolphins","sharks & rays")) %>%
                                     select(-c(V2, V4, V5, V8, V7))) %>%
  relocate(V6) %>%
  relocate(V3, .after=V6)

#Load the fasta file which is a single column of ASV names followed by their DNA sequence.
# We'll need to split this into 2 columns to make it easy to merge
sab24.12s.fasta <- readLines("data/2024Perley/12S/denoised/dna-sequences.fasta")


# extract IDs (remove ">")
ids <- sab24.12s.fasta[grepl("^>", sab24.12s.fasta)] |> gsub("^>", "",x=_)

# extract sequences (lines not starting with ">")
seqs <- sab24.12s.fasta[!grepl("^>", sab24.12s.fasta)]

# combine into a data frame
fas <- data.frame(ASV = ids, sequence = seqs)


#Now merge with filtered ASV file

sab24.with.seqs <- left_join(sab24.12s.filt, fas, by=c("OTU.ID"="ASV"))

write.csv(x = sab24.with.seqs, file = "data/2024Perley/12S/SAB2024_ASVs_WithSEQS.csv", row.names = F, quote=F)
