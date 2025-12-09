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


## COI 
#Redoing this filtration step from filters_ASVs.R but keeping OTU.ID
sab24.coi.filt <- filter_low_reads(sab24.coi.merge %>% 
                                     select(-c(V2:V11,V13,V22, V25, V28)) %>%
                                     filter(V29>0.979, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta",
                                                                  "Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria",
                                                                  "Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada",
                                                                  "Rotifera", "Chaetognatha","Kinorhyncha","Acanthocephala_thorny-headed_worms")) %>%
                                     rename(Phylum=V12, Class=V15, Species=V27) %>% 
                                     as.data.frame())

#Load the fasta file which is a single column of ASV names followed by their DNA sequence.
# We'll need to split this into 2 columns to make it easy to merge
sab24.coi.fasta <- readLines("data/2024Perley/COI/dna-sequences.fasta")


# extract IDs (remove ">")
ids <- sab24.coi.fasta[grepl("^>", sab24.coi.fasta)] |> gsub("^>", "",x=_)

# extract sequences (lines not starting with ">")
seqs <- sab24.coi.fasta[!grepl("^>", sab24.coi.fasta)]

# combine into a data frame
fas <- data.frame(ASV = ids, sequence = seqs)


#Now merge with filtered ASV file

sab24.COIwith.seqs <- left_join(sab24.coi.filt, fas, by=c("OTU.ID"="ASV"))

write.csv(x = sab24.COIwith.seqs, file = "data/2024Perley/COI/SAB2024_COIASVs_WithSEQS.csv", row.names = F, quote=F)
