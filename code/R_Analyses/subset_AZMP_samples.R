
# AZMP samples ------------------------------------------------------------

# filter out AZMP samples from 2024 Perley data - this includes the Gully and Browns Bank Line

# Load libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(tidyverse)

#simple function to filter reads in an ASV file - set to 0.01% now but can be changed to 0.1% etc
filter_low_reads <- function(df) {
  numeric_df <-df[sapply(df,is.numeric)]
  total_reads <- sum(numeric_df, na.rm = TRUE)  # Total reads across species and sites
  species_sums <- rowSums(numeric_df, na.rm = TRUE)  # Sum of reads per species
  df_filtered <- df[species_sums >= 0.0001 * total_reads, ]  # Keep species with at least 0.01% of total reads
  return(df_filtered)
}


#load metadata
# AZMP2024_EDNA_001 through _021 are for Fundian Channel (BBL and NEC)
# AZMP2024_EDNA_037 through _042 are the Gully (stations GUL02 and GUL03)
# AZMP2024_EDNA_052 through _063 are St. Anns Bank (STAB line)

azmp.metadat <- read.csv("data/2024Perley/GOTeDNA_SAB_Sampling_2024_Metadata.csv", header = T) %>% glimpse()

#Load eDNA datasets

azmp.fish <- read.csv("data/2024Perley/12S/GOTeDNA_SAB2024_Perley_12S_filtered.csv", header = T) %>% glimpse()

azmp.inverts <- read.csv("data/2024Perley/COI/GOTeDNA_SAB2024_COI_formatted.csv", header = T) %>% glimpse()


# Filter fish data first 
azmp.fish.filt <- azmp.fish %>% 
  dplyr::select(starts_with("V"),starts_with("AZMP")) %>%
  filter(rowSums(across(-c(1, 2))) > 1)
                
# Filter invertebrates
azmp.inverts.filt <- azmp.inverts %>% 
  dplyr::select(Species,starts_with("V"),starts_with("AZMP"))


# Get just the Gully data

azmp.fish.gully <- azmp.fish.filt %>%
  dplyr::select(V6,V3,AZMP2024.EDNA.037:AZMP2024.EDNA.042) %>%
  rename(Species=V6, percentMatch = V3) %>%
  filter(rowSums(across(-c(1, 2))) > 1) #surprise, there's barely anything in here! A few northern bottlenose and veiled anglemouth reads


azmp.inverts.gully <- azmp.inverts.filt %>%
  dplyr::select(Species,V29,AZMP2024.EDNA.037:AZMP2024.EDNA.042) %>%
  rename(confidence = V29) %>%
  filter(rowSums(across(-c(1, 2))) > 1) 

write.csv(azmp.inverts.gully, file = "data/2024Perley/COI/GullyInvertebrates.csv", row.names = F, quote=F)
