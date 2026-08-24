
# Prep data for GOTeDNA App -----------------------------------------------

#This document is in two parts, the first organizing 12S and COI data for the GOTeDNA templates, and the second part using those templates for OBIS 


#Open packages and load libraries
#install.packages("openxlsx")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("worrms")

library(readxl)
library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)
library(worrms)
library(stringr)
library(purrr)

#Read in excel files
mifish <- read_csv("data/2025RVsurvey/12S/RV2025_12S_filtered.csv")
coi <- read_csv("data/2025RVsurvey/COI/GOTeDNA-RV2025-COI-ASVsFiltered.csv")
SAB2024_meta <- read_csv("data/2025RVsurvey/2025RVsurvey_metadata_forGOTeDNA.csv") %>% glimpse() #Just keep this object as sab2024_meta to make downstream steps easier


#View(SAB_12S2024)
#View(SAB_COI2024)
#View(SAB2024_meta)

#Start with organizing GOTeDNA metadata sheet
#Create df with GOTeDNA columns and add rows for species names from bioinformatic output file (xlsx/csv)

df_Metadata_SAB2024_12S <- data.frame(
  Project_ID = character(length(SAB2024_meta$SampleID)),
  GOTeDNA_ID = character(length(SAB2024_meta$SampleID)),
  protocol_version = character(length(SAB2024_meta$SampleID)),
  protocolVersion_rationale = character(length(SAB2024_meta$SampleID)),
  ownerContact = character(length(SAB2024_meta$SampleID)),
  bibliographicCitation = character(length(SAB2024_meta$SampleID)),
  basisOfRecord = character(length(SAB2024_meta$SampleID)),
  materialSampleID = character(length(SAB2024_meta$SampleID)),
  eventID = character(length(SAB2024_meta$SampleID)),
  eventDate = character(length(SAB2024_meta$SampleID)),
  eventTime = character(length(SAB2024_meta$SampleID)),
  decimalLatitude = character(length(SAB2024_meta$SampleID)),
  decimalLongitude = character(length(SAB2024_meta$SampleID)),
  recordedBy = character(length(SAB2024_meta$SampleID)),
  sampleDepth = character(length(SAB2024_meta$SampleID)),
  waterTemp = character(length(SAB2024_meta$SampleID)),
  volumeFiltered = character(length(SAB2024_meta$SampleID)),
  dateFiltration = character(length(SAB2024_meta$SampleID)),
  timeFiltration = character(length(SAB2024_meta$SampleID)),
  controlType = character(length(SAB2024_meta$SampleID)),
  LClabel = character(length(SAB2024_meta$SampleID)),
  samplingStation = character(length(SAB2024_meta$SampleID)),
  waterColumnDepth = character(length(SAB2024_meta$SampleID)),
  tide = character(length(SAB2024_meta$SampleID)),
  depthWaterTemp = character(length(SAB2024_meta$SampleID)),
  totalDNAconc = character(length(SAB2024_meta$SampleID)),
  unitsDNAconc = character(length(SAB2024_meta$SampleID)),
  stringsAsFactors = FALSE
)


#View dataframe
#View(df_Metadata_SAB2024_12S)

#Dataframe created! Next, we need to fill it with the ESI Coastal 2024 metadata file

df_Metadata_SAB2024_12S <- as.data.frame(SAB2024_meta)  %>%
    mutate(
    Project_ID = "Ecosystem Survey and AZMP 2025",
    GOTeDNA_ID = "22",
    protocol_version = "1",
    protocolVersion_rationale = NA_character_,
    ownerContact = "nick.jeffery@dfo-mpo.gc.ca",
    bibliographicCitation = NA_character_,
    basisOfRecord = "materialSample",
    materialSampleID = SampleID,
    eventID =  paste(GOTeDNA_ID, SampleID, sep = "-"),
    eventDate = eventDate,
    eventTime = eventTime,
    decimalLatitude = latitude,
    decimalLongitude = longitude,
    recordedBy = recordedBy,
    sampleDepth = sampleDepth,
    waterTemp = waterTemp,
    volumeFiltered = volumeFiltered,
    dateFiltration = NA_character_,
    timeFiltration = NA_character_,
    controlType = controlType,
    LClabel = "We acknowledge that this work has taken place in Mi’kma’ki, the ancestral and unceded territory of the Mi’kmaq People, who are part of the Wabanaki (Dawnland Confederacy). Our relationship to this place is one of a visiting scientist, learning from its waters. Our commitment is to move beyond words by ensuring scientific findings are shared publicly and accompanied by this acknowledgement to inform both science and a respectful relationship to the Mi’kmaq communities who have stewarded these lands and waters since time immemorial. We acknowledge their enduring sovereignty, knowledge, and contributions, past, present, and future.",
    samplingStation = samplingStation,
    waterColumnDepth = NA_character_,
    tide = tide,
    depthWaterTemp = depthWaterTemp,
    totalDNAconc = totalDNAconc,
    unitsDNAconc = unitsDNAconc,
  ) %>%
  #filter(MPA == "St. Anns Bank") %>%
  select(
    Project_ID, GOTeDNA_ID, protocol_version, protocolVersion_rationale, ownerContact,
    bibliographicCitation, basisOfRecord, materialSampleID, eventID, eventDate, eventTime,
    decimalLatitude, decimalLongitude, recordedBy, sampleDepth, waterTemp, volumeFiltered,
    dateFiltration, timeFiltration, controlType, LClabel, samplingStation, waterColumnDepth,
    tide, depthWaterTemp, totalDNAconc, unitsDNAconc
  )


df_Metadata_Gully2024_12S <- SAB2024_meta  %>%
  
  mutate(
    Project_ID = "Gully 2024",
    GOTeDNA_ID = "22",
    protocol_version = "1",
    protocolVersion_rationale = NA_character_,
    ownerContact = "nick.jeffery@dfo-mpo.gc.ca",
    bibliographicCitation = NA_character_,
    basisOfRecord = "materialSample",
    materialSampleID = SampleID,
    eventID =  paste(GOTeDNA_ID, SampleID, sep = "-"),
    eventDate = eventDate,
    eventTime = eventTime,
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    recordedBy = recordedBy,
    sampleDepth = sampleDepth,
    waterTemp = waterTemp,
    volumeFiltered = volumeFiltered,
    dateFiltration = NA_character_,
    timeFiltration = NA_character_,
    controlType = controlType,
    LClabel = NA_character_,
    samplingStation = samplingStation,
    waterColumnDepth = NA_character_,
    tide = tide,
    depthWaterTemp = depthWaterTemp,
    totalDNAconc = totalDNAconc,
    unitsDNAconc = unitsDNAconc,
  ) %>%
  filter(MPA == "Gully") %>%
  select(
    Project_ID, GOTeDNA_ID, protocol_version, protocolVersion_rationale, ownerContact,
    bibliographicCitation, basisOfRecord, materialSampleID, eventID, eventDate, eventTime,
    decimalLatitude, decimalLongitude, recordedBy, sampleDepth, waterTemp, volumeFiltered,
    dateFiltration, timeFiltration, controlType, LClabel, samplingStation, waterColumnDepth,
    tide, depthWaterTemp, totalDNAconc, unitsDNAconc
  )


df_Metadata_Fundian2024_12S <- SAB2024_meta  %>%
  
  mutate(
    Project_ID = "Fundian Channel - Browns Bank 2024",
    GOTeDNA_ID = "22",
    protocol_version = "1",
    protocolVersion_rationale = NA_character_,
    ownerContact = "nick.jeffery@dfo-mpo.gc.ca",
    bibliographicCitation = NA_character_,
    basisOfRecord = "materialSample",
    materialSampleID = SampleID,
    eventID =  paste(GOTeDNA_ID, SampleID, sep = "-"),
    eventDate = eventDate,
    eventTime = eventTime,
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    recordedBy = recordedBy,
    sampleDepth = sampleDepth,
    waterTemp = waterTemp,
    volumeFiltered = volumeFiltered,
    dateFiltration = NA_character_,
    timeFiltration = NA_character_,
    controlType = controlType,
    LClabel = NA_character_,
    samplingStation = samplingStation,
    waterColumnDepth = NA_character_,
    tide = tide,
    depthWaterTemp = depthWaterTemp,
    totalDNAconc = totalDNAconc,
    unitsDNAconc = unitsDNAconc,
  ) %>%
  filter(MPA == "Fundian Channel - Browns Bank") %>%
  select(
    Project_ID, GOTeDNA_ID, protocol_version, protocolVersion_rationale, ownerContact,
    bibliographicCitation, basisOfRecord, materialSampleID, eventID, eventDate, eventTime,
    decimalLatitude, decimalLongitude, recordedBy, sampleDepth, waterTemp, volumeFiltered,
    dateFiltration, timeFiltration, controlType, LClabel, samplingStation, waterColumnDepth,
    tide, depthWaterTemp, totalDNAconc, unitsDNAconc
  )


# View dataframes
#View(df_Metadata_SAB2024_12S)
#View(df_Metadata_Gully2024_12S)
#View(df_Metadata_Fundian2024_12S)

#Combine all dataframes
GOTeDNA_22_metadata <- df_Metadata_SAB2024_12S
  #rbind(df_Metadata_SAB2024_12S, df_Metadata_Gully2024_12S, df_Metadata_Fundian2024_12S)
#View(GOTeDNA_22_metadata)

# Change commas in recordedBy column to "|"
GOTeDNA_22_metadata$recordedBy <- gsub(",", " |", GOTeDNA_22_metadata$recordedBy)

#change underscores to periods in materialSampleID
GOTeDNA_22_metadata$materialSampleID <- gsub("-", ".", GOTeDNA_22_metadata$materialSampleID)

#View(GOTeDNA_22_metadata)


#Write excel file and input into Google sheets GOTeDNA-23

writexl::write_xlsx(
  GOTeDNA_22_metadata,
  "data/2025RVsurvey/GOTeDNA-22_Sample-Metadata_12SCOI_RVsurvey2025_Final.xlsx"
)


#Now organize the sample_metabarcoding data
#Create df with GOTeDNA columns and add rows for species names from bioinformatic output file (xlsx/csv)

GOTeDNA_22_12S<- data.frame(
  Project_ID = character(length(mifish$X6)),
  GOTeDNA_ID = character(length(mifish$X6)),
  GOTeDNA_version = character(length(mifish$X6)),
  basisOfRecord = character(length(mifish$X6)),
  materialSampleID = character(length(mifish$X6)),
  eventID = character(length(mifish$X6)),
  recordedBy = character(length(mifish$X6)),
  target_gene = character(length(mifish$X6)),
  target_subfragment = character(length(mifish$X6)),
  scientificName = mifish$X6,
  taxonID = character(length(mifish$X6)),
  kingdom = character(length(mifish$X6)),
  phylum = character(length(mifish$X6)),
  class = character(length(mifish$X6)),
  order = character(length(mifish$X6)),
  family = character(length(mifish$X6)),
  genus = character(length(mifish$X6)),
  occurrenceID = character(length(mifish$X6)),
  organismQuantity = numeric(length(mifish$X6)),
  organismQuantityType = character(length(mifish$X6)),
  sampleSizeValue = character(length(mifish$X6)),
  samplingSizeUnit = character(length(mifish$X6)),
  samplingProtocol = character(length(mifish$X6)),
  associatedSequences = character(length(mifish$X6)),
  identificationRemarks = character(length(mifish$X6)),
  stringsAsFactors = FALSE
)

#View dataframe
#View(GOTeDNA_22_12S)

# Pivot from wide (samples as columns) to long (sample, value)
# Linking bioinformatic output to GOTeDNA file format
GOTeDNA_22_12S <- mifish %>%
  pivot_longer(
    cols = -c(X6, X3),   # keep Species fixed, everything else pivots
    names_to = "materialSampleID",
    values_to = "organismQuantity"
  ) %>%
  mutate(
    Project_ID = NA_character_,
    GOTeDNA_ID = "22",
    GOTeDNA_version = "1",
    basisOfRecord = "MaterialSample",
    eventID = paste(GOTeDNA_ID, materialSampleID, sep = "-"),
    recordedBy = NA_character_,
    target_gene = "12S",
    target_subfragment = "MiFish_12S",
    scientificName = X6,  #Species name
    taxonID = NA_character_,
    kingdom = NA_character_,
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    occurrenceID = NA_character_,
    organismQuantityType = "DNA sequence reads",
    sampleSizeValue = NA_character_,
    sampleSizeUnit = NA_character_,
    samplingProtocol = NA_character_,
    associatedSequences = NA_character_,
    identificationRemarks = X3,
  ) %>%
  select(
    Project_ID, GOTeDNA_ID, GOTeDNA_version, basisOfRecord, materialSampleID, eventID,
    recordedBy, target_gene, target_subfragment, scientificName, taxonID,
    kingdom, phylum, class, order, family, genus, occurrenceID,
    organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,
    samplingProtocol, associatedSequences, identificationRemarks
  )



# 1) Build a one-time lookup for unique species -> AphiaID (fast & robust)
species_lookup <- tibble::tibble(scientificName = unique(GOTeDNA_22_12S$scientificName)) %>%
  mutate(
    taxonID = sapply(scientificName, function(s) {
      id <- tryCatch(wm_name2id(s, marine_only = TRUE), error = function(e) NA)
      if (is.null(id) || length(id) == 0) NA_integer_ else as.integer(id[1])
    })
  )

# 2) Join the AphiaIDs into your df
GOTeDNA_22_12S <- GOTeDNA_22_12S %>%
  select(-any_of("taxonID")) %>%        # drop old taxonID column if present
  left_join(species_lookup, by = "scientificName") %>%
  filter(!is.na(taxonID))               # 3) remove rows with no AphiaID match

# Optional: see how many were dropped
dropped <- sum(is.na(species_lookup$taxonID))
message(sprintf("Dropped %d species with no AphiaID match.", dropped))


# Pick one AphiaID you have and inSpeciesect what wm_classification returns
one_id <- unique(GOTeDNA_22_12S$taxonID)[1]
wm_classification(one_id) %>% str()   # you'll see columns 'rank' and 'scientificname'


# Safe extractor for a given rank using the 'scientificname' column
.rank_val <- function(classif, rank_name) {
  if (is.null(classif)) return(NA_character_)
  if (!all(c("rank","scientificname") %in% names(classif))) return(NA_character_)
  v <- classif$scientificname[classif$rank == rank_name]
  if (length(v) == 0) NA_character_ else as.character(v[1])
}

get_taxonomy_for_id <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  tibble(
    taxonID = id,
    kingdom = .rank_val(classif, "Kingdom"),
    phylum  = .rank_val(classif, "Phylum"),
    class   = .rank_val(classif, "Class"),
    order   = .rank_val(classif, "Order"),
    family  = .rank_val(classif, "Family"),
    genus   = .rank_val(classif, "Genus")
  )
}

# Build a lookup for the unique AphiaIDs present in df
unique_ids <- unique(GOTeDNA_22_12S$taxonID)
taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)

# (Optional) see how many came back totally empty
# sum(rowSums(is.na(taxonomy_lookup[,c("kingdom","phylum","class","order","family","genus")])) == 6)

# Join back to your df
GOTeDNA_22_12S <- GOTeDNA_22_12S %>%
  select(-any_of(c("kingdom","phylum","class","order","family","genus"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")



#Fill out and edit other data columns with project details
GOTeDNA_22_12S <- GOTeDNA_22_12S %>%
  mutate(
    occurrenceID = paste(
      GOTeDNA_ID, materialSampleID, target_gene, substr(genus, 1, 4),
      sep = "-"
    )
  )

#Order the columns to be in GOTeDNA Order
GOTeDNA_22_12S <- GOTeDNA_22_12S %>%
  select(
    GOTeDNA_ID,              # Protocol ID
    GOTeDNA_version,         # protocol version
    basisOfRecord,           # basis of record
    materialSampleID,        # material sample
    eventID,                 # event ID
    recordedBy,              # recorded by
    target_gene,             # target gene
    target_subfragment,      # target subfragment
    scientificName,          # scientific name
    taxonID,                 # taxon ID
    kingdom,                 # kingdom
    phylum,                  # phylum
    class,                   # class
    order,                   # order
    family,                  # family
    genus,                   # genus
    occurrenceID,            # occurrence ID
    organismQuantity,        # organism quantity
    organismQuantityType,    # organism quantity type
    sampleSizeValue,         # sample size value
    sampleSizeUnit,          # sample size unit
    samplingProtocol,        # sampling protocol
    associatedSequences,     # associated sequences
    identificationRemarks    # identification remarks
  )

#View df
View(GOTeDNA_22_12S)


#Do the same for COI

GOTeDNA_22_COI <- data.frame(
  Project_ID = character(length(coi$X6)),
  GOTeDNA_ID = character(length(coi$X6)),
  GOTeDNA_version = character(length(coi$X6)),
  basisOfRecord = character(length(coi$X6)),
  materialSampleID = character(length(coi$X6)),
  eventID = character(length(coi$X6)),
  recordedBy = character(length(coi$X6)),
  target_gene = character(length(coi$X6)),
  target_subfragment = character(length(coi$X6)),
  scientificName = coi$X6,
  taxonID = character(length(coi$X6)),
  kingdom = character(length(coi$X6)),
  phylum = character(length(coi$X6)),
  class = character(length(coi$X6)),
  order = character(length(coi$X6)),
  family = character(length(coi$X6)),
  genus = character(length(coi$X6)),
  occurrenceID = character(length(coi$X6)),
  organismQuantity = numeric(length(coi$X6)),
  organismQuantityType = character(length(coi$X6)),
  sampleSizeValue = character(length(coi$X6)),
  samplingSizeUnit = character(length(coi$X6)),
  samplingProtocol = character(length(coi$X6)),
  associatedSequences = character(length(coi$X6)),
  identificationRemarks = character(length(coi$X6)),
  stringsAsFactors = FALSE
)

#View(GOTeDNA_22_COI)


# Pivot from wide (samples as columns) to long (sample, value)
# Linking bioinformatic output to GOTeDNA file format
GOTeDNA_22_COI <- coi %>%
  pivot_longer(
    cols = -c(X6, X3),   # keep Species fixed, everything else pivots
    names_to = "materialSampleID",
    values_to = "organismQuantity"
  ) %>%
  mutate(
    Project_ID = NA_character_,
    GOTeDNA_ID = "22",
    GOTeDNA_version = "1",
    basisOfRecord = "MaterialSample",
    eventID = paste(GOTeDNA_ID, materialSampleID, sep = "-"),
    recordedBy = NA_character_,
    target_gene = "COI",
    target_subfragment = "Leray_COI",
    scientificName = X6,  #Species name
    taxonID = NA_character_,
    kingdom = NA_character_,
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    occurrenceID = NA_character_,
    organismQuantityType = "DNA sequence reads",
    sampleSizeValue = NA_character_,
    sampleSizeUnit = NA_character_,
    samplingProtocol = NA_character_,
    associatedSequences = NA_character_,
    identificationRemarks = X3,
  ) %>%
  select(
    Project_ID, GOTeDNA_ID, GOTeDNA_version, basisOfRecord, materialSampleID, eventID,
    recordedBy, target_gene, target_subfragment, scientificName, taxonID,
    kingdom, phylum, class, order, family, genus, occurrenceID,
    organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,
    samplingProtocol, associatedSequences, identificationRemarks
  )

# View dataframe
#View(GOTeDNA_22_COI)

# 1) Build a one-time lookup for unique species -> AphiaID (fast & robust)
species_lookup <- tibble::tibble(scientificName = unique(GOTeDNA_22_COI$scientificName)) %>%
  mutate(
    taxonID = sapply(scientificName, function(s) {
      id <- tryCatch(wm_name2id(s, marine_only = TRUE), error = function(e) NA)
      if (is.null(id) || length(id) == 0) NA_integer_ else as.integer(id[1])
    })
  )

# 2) Join the AphiaIDs into your df
GOTeDNA_22_COI <- GOTeDNA_22_COI %>%
  select(-any_of("taxonID")) %>%        # drop old taxonID column if present
  left_join(species_lookup, by = "scientificName") %>%
  filter(!is.na(taxonID))               # 3) remove rows with no AphiaID match

# Optional: see how many were dropped
dropped <- sum(is.na(species_lookup$taxonID))
message(sprintf("Dropped %d species with no AphiaID match.", dropped))


# Pick one AphiaID you have and inSpeciesect what wm_classification returns
one_id <- unique(GOTeDNA_22_COI$taxonID)[1]
wm_classification(one_id) %>% str()   # you'll see columns 'rank' and 'scientificname'


# Safe extractor for a given rank using the 'scientificname' column
.rank_val <- function(classif, rank_name) {
  if (is.null(classif)) return(NA_character_)
  if (!all(c("rank","scientificname") %in% names(classif))) return(NA_character_)
  v <- classif$scientificname[classif$rank == rank_name]
  if (length(v) == 0) NA_character_ else as.character(v[1])
}

get_taxonomy_for_id <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  tibble(
    taxonID = id,
    kingdom = .rank_val(classif, "Kingdom"),
    phylum  = .rank_val(classif, "Phylum"),
    class   = .rank_val(classif, "Class"),
    order   = .rank_val(classif, "Order"),
    family  = .rank_val(classif, "Family"),
    genus   = .rank_val(classif, "Genus")
  )
}

# Build a lookup for the unique AphiaIDs present in df
unique_ids <- unique(GOTeDNA_22_COI$taxonID)
taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)

# (Optional) see how many came back totally empty
# sum(rowSums(is.na(taxonomy_lookup[,c("kingdom","phylum","class","order","family","genus")])) == 6)

# Join back to your df
GOTeDNA_22_COI <- GOTeDNA_22_COI %>%
  select(-any_of(c("kingdom","phylum","class","order","family","genus"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")


#Fill out and edit other data columns with project details
GOTeDNA_22_COI <- GOTeDNA_22_COI %>%
  mutate(
    occurrenceID = paste(
      GOTeDNA_ID, materialSampleID, target_gene, substr(genus, 1, 4),
      sep = "-"
    )
  )

#Order the columns to be in GOTeDNA Order
GOTeDNA_22_COI <- GOTeDNA_22_COI %>%
  select(
    GOTeDNA_ID,              # Protocol ID
    GOTeDNA_version,         # protocol version
    basisOfRecord,           # basis of record
    materialSampleID,        # material sample
    eventID,                 # event ID
    recordedBy,              # recorded by
    target_gene,             # target gene
    target_subfragment,      # target subfragment
    scientificName,          # scientific name
    taxonID,                 # taxon ID
    kingdom,                 # kingdom
    phylum,                  # phylum
    class,                   # class
    order,                   # order
    family,                  # family
    genus,                   # genus
    occurrenceID,            # occurrence ID
    organismQuantity,        # organism quantity
    organismQuantityType,    # organism quantity type
    sampleSizeValue,         # sample size value
    sampleSizeUnit,          # sample size unit
    samplingProtocol,        # sampling protocol
    associatedSequences,     # associated sequences
    identificationRemarks    # identification remarks
  )

#View df
#View(GOTeDNA_22_COI)


#Combine 12S and COI df
GOTeDNA_22 <- rbind(GOTeDNA_22_12S, GOTeDNA_22_COI)
#View(GOTeDNA_22)


#Write excel file and input into Google sheets GOTeDNA-22

writexl::write_xlsx(
  GOTeDNA_22,
  "data/2025RVsurvey/GOTeDNA-22_Sample-Metabarcoding_12SCOIRV2025_Final.xlsx"
)


##COMPLETE

## OBIS MCT St Anns Bank, Gully, Fundian Channel-Browns Bank 12S

#Need to have the precursor R file because the data is too large to store in Google Drive 

#Install packages and Load libraries
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(worrms)
library(openxlsx)
library(lubridate)
library(tidyr)

##Obtain filled GOTeDNA template and link data to OBIS specific column headers

setwd("C:/Users/HEADK/Desktop/OBIS_prep")

GOT22 <- read.xlsx("GOTeDNA-22_data.xlsx", sheet = "Sample_Metadata")
GOT22_metabar <- GOTeDNA_22

#Standardize materialSampleID values

GOT22 <- GOT22 %>%
  mutate(
    materialSampleID = as.character(materialSampleID),
    materialSampleID = sub("\\.0$", "", materialSampleID),
    materialSampleID = gsub("_", ".", materialSampleID)
  )

GOT22_metabar <- GOT22_metabar %>%
  mutate(
    materialSampleID = as.character(materialSampleID),
    materialSampleID = sub("\\.0$", "", materialSampleID),
    materialSampleID = gsub("_", ".", materialSampleID)
  )

#Filter for data contact and primer
# Filter for data contact and sample IDs #Find polygons and fix coordinates here
GOT22_f <- GOT22 %>%
  filter(
    ownerContact == "nick.jeffery@dfo-mpo.gc.ca",
    #decimalLatitude  >= 44.80 & decimalLatitude <= 45.27,   #Filter for coordinates specific to SAB (to exclude Gully and Fundian data)
    #decimalLongitude >= -68.78 & decimalLongitude <= -66.07
  )



# Filter the metabarcoding dataset
GOT22_metabar_f <- GOT22_metabar %>%
  filter(
    target_gene == "12S",
  )


#Combine files
GOT22_joined <- GOT22_f %>%
  left_join(
    GOT22_metabar_f,
    by = "materialSampleID",
    relationship = "many-to-many"
  )

#Fix the values in eventDate

GOT22_joined <- GOT22_joined %>%
  mutate(
    # Always work from a character version
    eventDate_chr = as.character(eventDate),
    
    # Try to interpret as numeric (Excel serial)
    eventDate_num = suppressWarnings(as.numeric(eventDate_chr)),
    
    # Fix 2-digit year dates like 21/8/19 -> 21/8/2019
    eventDate_chr_2y = if_else(
      str_detect(eventDate_chr, "^\\d{1,2}/\\d{1,2}/\\d{2}$"),
      sub("(\\d{1,2}/\\d{1,2}/)(\\d{2})$", "\\120\\2", eventDate_chr),
      NA_character_
    ),
    
    eventDate = case_when(
      # 1) Excel numeric serials (42293 etc.)
      !is.na(eventDate_num) ~ as.Date(eventDate_num, origin = "1899-12-30"),
      
      # 2) dd/mm/yy we just expanded to dd/mm/20yy
      !is.na(eventDate_chr_2y) ~ as.Date(eventDate_chr_2y, format = "%d/%m/%Y"),
      
      # 3) Normal dd/mm/yyyy strings
      TRUE ~ as.Date(eventDate_chr, format = "%d/%m/%Y")
    )
  ) %>%
  select(-eventDate_chr, -eventDate_num, -eventDate_chr_2y)

#Create unique event and occurrence IDs
# Create unique event and occurrence IDs for non-controls with coordinates
GOT22_joined <- GOT22_joined %>%
  # 1) Remove any rows that are controls (controlType not NA/empty)
  filter(is.na(controlType) | controlType == "") %>%
  
  # 2) Drop any rows missing coordinates
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  ) %>%
  
  # 3) Work within each original eventID.x
  group_by(eventID.x) %>%
  mutate(
    row_suffix = row_number(),
    
    # eventID for every remaining row (all have coords now)
    eventID = paste0(eventID.x, "-", row_suffix),
    
    # occurrenceID for every remaining row (present + absent)
    occurrenceID = if_else(
      is.na(occurrenceID),
      paste0("DFO-MCT-PerleyOffshore_2024", materialSampleID, "-", row_suffix),
      paste0(occurrenceID, "-", row_suffix)
    )
  ) %>%
  ungroup() %>%
  select(-row_suffix)


#Set up occurrence dataframe

map_raw <- read.csv("OBIS_Perley_2024_12S.csv", stringsAsFactors = FALSE)

map_long <- bind_rows(
  map_raw %>%
    select(Seq_Run_ID, lib_id = Lib_ID) %>%
    mutate(lib_id = gsub("_", ".", lib_id)),
  
  map_raw %>%
    select(Seq_Run_ID = Seq_Run_ID.1,
           lib_id    = Lib_ID.1) %>%
    mutate(lib_id = gsub("_", ".", lib_id))
) %>%
  filter(!is.na(lib_id) & lib_id != "")

map_clean <- map_long %>%
  filter(!str_detect(lib_id, "ENEG")) %>%
  mutate(
    # Extract SAB24.EDNA.014 from SAB24.EDNA.014.S1.MiSeq_250131
    materialSampleID = str_extract(lib_id, "^[A-Za-z0-9]+\\.EDNA\\.\\d{3}")
  ) %>%
  distinct(materialSampleID, .keep_all = TRUE)

map_clean <- map_clean %>%
  # treat blanks as NA
  mutate(Seq_Run_ID = na_if(Seq_Run_ID, "")) %>%
  # fill seq_run_id down within the whole table
  fill(Seq_Run_ID)

occurrence <- GOT22_joined %>%
  mutate(
    GOTeDNAprotocol_ID= GOTeDNA_ID.x,
    GOTeDNAprotocolVersion = protocolVersion,
    project_contact = ownerContact,
    bibliographicCitation = bibliographicCitation,
    basisOfRecord = basisOfRecord.x,
    locationID = samplingStation,
    materialSampleID = materialSampleID,
    samp_name = materialSampleID,
    technical_rep_id = "1",
    recordedBy = recordedBy.x,
    geo_loc_name = "Canada: Nova Scotia, Offshore",
    nameAccordingTo = "WoRMS",
    country = "Canada",
    datasetID = "DFO-MCT-SAB-GUL-FUN-2024",
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    minimumDepthInMeters = sampleDepth,
    maximumDepthInMeters = sampleDepth,
    totalDNAconc = totalDNAconc,
    unitsDNAconc = unitsDNAconc,
    samp_size = "3",
    samp_size_unit = "L",
    filtrationType = "peristaltic",
    size_frac = "1.2",
    filter_name = "Smith-Root",
    eventDate = eventDate,
    occurrenceID = occurrenceID,
    occurrenceStatus = case_when(
      organismQuantity > 0 ~ "present",
      organismQuantity == 0 ~ "absent",
      TRUE ~ NA_character_
    ),
    language = "en",
    filter_material = "Polyethersulfone (PES)",
    materialSampleID = materialSampleID,
    month = month(eventDate),
    bibliographicCitation = bibliographicCitation,
    year = year(eventDate),
    occurrenceID = occurrenceID,
    samp_name = materialSampleID,
    project_name = "Perley Offshore 2024",
    scientificName = scientificName,
    env_broad_scale = "marine biome (ENVO:00000447)",
    env_local_scale = "oceanic benthopelagic zone biome (ENVO:01000040)",
    env_medium = "sea water (ENVO:00002149)",
    targetTaxonomicAssay = "vertebrate",
    samp_store_sol = "100% ethanol",
    nucl_acid_ext_kit = "DNeasy® Blood & Tissue Kit (Qiagen)",
    nucl_acid_ext_modify = "For the extractions, all sample filters were cut in half and processed according to ABL’s eDNA extraction protocol. This particular protocol utilized the QIACube Connect (Qiagen) and a modified version of the Qiagen DNeasy Blood and Tissue Kit (240)",
    seq_kit = "MiSeq Reagent Kit v3 (Illumina)",
    samp_store_temp = "-80",
    samp_collect_device = "Niskin bottle",
    ref_db = "NCBI nt_euk",
    otu_seq_comp_appr = "BLAST v2.15.0",
    otu_db = "QIIME2 Amplicon v2024.10, DADA2 v2024.10.0 (in R, DADA2 v1.30.0)",
    assay_name = "12S (metabarcoding)",
    assay_type = "metabarcoding",
    target_gene = target_gene,
    target_subfragment = target_subfragment,
    instrument = "Illumina MiSeq [OBI_0002003]",
    pcr_primer_forward = "CGTGCCAGCCACCGCGGTT",
    pcr_primer_reverse = "CATAGTGGGGTATCTAATCCCAGTTTG",
    pcr_primer_name_forward = "12S_248F_RADS_For",
    pcr_primer_name_reverse = "Mifish_UR_Miya",
    pcr_primer_reference = "doi.org/10.1139/cjfas-2021-0215 | doi.org/10.1098/RSOS.150088",
    identificationRemarks = identificationRemarks,
    checkls_ver = "1.0.2",
    sampleSizeValue = sampleSizeValue,
    sampleSizeUnit = sampleSizeUnit,
    filter_passive_active_0_1 = "1",
    associatedSequences = "In progress for NCBI",
    samp_category = "sample",
    samp_collec_device = samp_collect_device,
    project_id = Project_ID,
    pcr_0_1 = "1",
    samp_mat_process = "filtration",
    platform = "ILLUMINA",
    tax_assign_cat = "sequence similarity"
  ) %>%
  
  #filter(
  #!is.na(minimumDepthInMeters),
  #!is.na(maximumDepthInMeters),
  #!is.na(target_gene),
  #!is.na(scientificName)
  #) %>%
  left_join(map_clean, by = "materialSampleID"
  ) %>%
  mutate (seq_run_id = Seq_Run_ID)


#Look up AphiaID (taxonID) from WoRMS by ScientificName

species_lookup <- occurrence %>%
  distinct(scientificName) %>%
  mutate(
    taxonID = map_int(
      scientificName,
      ~ {
        # Keep NAs as NAs
        if (is.na(.x) || .x == "") return(NA_integer_)
        
        id <- tryCatch(
          wm_name2id(.x, marine_only = FALSE),
          error = function(e) NA
        )
        
        if (is.null(id) || length(id) == 0 || is.na(id[1])) {
          NA_integer_
        } else {
          as.integer(id[1])
        }
      }
    )
  )

# Join taxonID into occurrence, but DO NOT drop rows with no match
occurrence <- occurrence %>%
  select(-any_of(c("taxonID", "scientificNameID",
                   "kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "aphiaID"))) %>%  # we'll re-create these
  left_join(species_lookup, by = "scientificName") %>%
  mutate(
    # optional helper flag: which rows had a WoRMS match?
    worms_match = !is.na(taxonID),
    # keep AphiaID explicitly if you want it separate from taxonID
    aphiaID = taxonID
  )


# Build a taxonomy lookup for each unique non-NA AphiaID


# Helper to pull a rank from wm_classification()
.rank_val <- function(classif, rank_name) {
  if (is.null(classif)) return(NA_character_)
  if (!all(c("rank", "scientificname") %in% names(classif))) return(NA_character_)
  v <- classif$scientificname[classif$rank == rank_name]
  if (length(v) == 0) NA_character_ else as.character(v[1])
}

get_taxonomy_for_id <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  rec     <- tryCatch(wm_record(id),        error = function(e) NULL)
  
  tibble(
    taxonID           = id,
    scientificNameID  = if (!is.null(rec) && "lsid" %in% names(rec))
      rec$lsid else NA_character_,
    kingdom           = .rank_val(classif, "Kingdom"),
    phylum            = .rank_val(classif, "Phylum"),
    class             = .rank_val(classif, "Class"),
    order             = .rank_val(classif, "Order"),
    family            = .rank_val(classif, "Family"),
    genus             = .rank_val(classif, "Genus")
  )
}

unique_ids <- occurrence$taxonID %>%
  unique() %>%
  discard(is.na)

taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)


#Join WoRMS taxonomy back to occurrence

occurrence <- occurrence %>%
  select(-any_of(c("kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "scientificNameID"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")


#This code is used to extract information from otus.database.fasta.gz files (barque) to obtain the OTUids and sequences for input into the DNA-derived data extension for GBIF/OBIS

#install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)

# 1. Read the OTU database File
otu_df <- read.csv("SAB2024_ASVs_WithSEQS_12S.csv", stringsAsFactors = FALSE)

#DNA derived data extension

#Link DNA sequences from otu_df

dna_df <- occurrence %>%
  left_join(
    otu_df %>%
      select(
        scientificName = Species,
        seq_id = OTU.ID,
        dna_sequence = sequence
      ),
    relationship = "many-to-many"
  )

###########################################
#Run Quality Control checks

#install.packages("devtools")
#devtools::install_github("iobis/obistools")

library(obistools)
library(Hmisc)

occur <- occurrence
dna <- dna_df

#Check that the taxa names match with WoRMS
#worms <- match_taxa(unique(occur$scientificName)) #choose option 1 each time

occur <- merge(occur, worms, by="scientificName")
colnames(occur)[colnames(occur) == "scientificNameID.y"] = "scientificNameID"

#Check that all required fields are present in the occurrence table
#Check the uniqueness of the occurrenceID field (Want to = TRUE)
length(occur$occurrenceID) == length(unique(occur$occurrenceID))

#clean files
cleaned_occ <- occur %>%
  select(
    occurrenceID, bibliographicCitation, materialSampleID, eventDate, decimalLatitude, decimalLongitude, scientificName,
    organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,  associatedSequences,
    basisOfRecord, locationID,  recordedBy, country, datasetID, occurrenceStatus, minimumDepthInMeters, maximumDepthInMeters,
    language,  month, year, taxonID, scientificNameID, kingdom, phylum, class, order, family, genus
  )

cleaned_dna <- dna %>%
  select(
    project_name, dna_sequence, target_gene, pcr_primer_forward, pcr_primer_reverse, samp_name,
    env_broad_scale, env_local_scale, env_medium, samp_mat_process, size_frac, samp_size,
    samp_size_unit, otu_db, seq_kit, otu_seq_comp_appr, pcr_primer_name_forward, pcr_primer_name_reverse,
    pcr_primer_reference,  occurrenceID
  )

##Write updated + cleaned files; When confident in the data, write these files and put into OBIS
#Write CSV files for input into OBIS
write.csv(cleaned_occ, "OBIS_MCT_SABGULFUN2024_12S_occurrence.csv")
write.csv(cleaned_dna, "OBIS_MCT_SABGULFUN2024_12S_dnaderiveddata.csv")


#Create eMOF df
#I need to create code to organize the measurement or fact df -> see OBIS Manual for details

emof <- dna_df %>%
  mutate(
    occurrenceID,
    measurementType = NA_character_, #This is what I put all of the FAIRe names in
    measurementValue = NA_character_,
    measurementUnit = NA_character_,
    measurementTypeID = NA_character_, #This is where I cite all of the FAIRe names
    measurementValueID = NA_character_,
    measurementUnitID = NA_character_,
    measurementRemarks = NA_character_
  ) %>%
  select(
    occurrenceID, measurementType, measurementValue, measurementUnit, measurementTypeID, measurementValueID, measurementUnitID, measurementUnitID, measurementRemarks
  )


#I need to create code that will provide information for each FAIRe name in the measurementType column per occurrenceID
#Insert proper FAIRe names from occurrence below

url_map <- c(
  LClabel                  = "https://github.com/GOTeDNA-OBON",
  ownerContact             = "https://github.com/GOTeDNA-OBON",
  totalDNAconc             = "https://github.com/GOTeDNA-OBON",
  unitsDNAconc             = "https://github.com/GOTeDNA-OBON",
  dateFiltration           = "https://github.com/GOTeDNA-OBON",
  timeFiltration           = "https://github.com/GOTeDNA-OBON",
  depthWaterTemp           = "https://github.com/GOTeDNA-OBON",
  occurrenceID             = "https://manual.obis.org/darwin_core.html"
)

emof <- dna_df %>%
  select(
    seq_id, samp_category, checkls_ver, assay_name, assay_type, targetTaxonomicAssay,
    geo_loc_name, technical_rep_id, project_contact, seq_run_id, lib_id, project_id,
    pcr_0_1, samp_store_sol, samp_store_temp, platform, instrument, tax_assign_cat,
    LClabel, occurrenceID, nucl_acid_ext_kit, filter_material
  ) %>%
  distinct() %>%
  mutate(across(-occurrenceID, as.character)) %>%
  pivot_longer(
    cols = -occurrenceID,
    names_to = "measurementType",
    values_to = "measurementValue"
  ) %>%
  mutate(
    measurementTypeID = dplyr::recode(
      measurementType,
      !!!url_map,
      .default = "https://github.com/FAIR-eDNA/FAIRe_checklist"
    ),
    measurementUnit    = NA_character_,
    measurementValueID = NA_character_,
    measurementUnitID  = NA_character_,
    measurementRemarks = NA_character_
  )

#write csv file and put into OBIS
write.csv(emof, "OBIS_MCT_SABGULFUN2024_12S_emof.csv")

#---------------------------------------------------------------------------------------------------------------------------------------

#Same for COI

## OBIS MCT St Anns Bank, Gully, Fundian Channel-Browns Bank 12S

#Need to have the precursor R file because the data is too large to store in Google Drive 

#Install packages and Load libraries
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(worrms)
library(openxlsx)
library(lubridate)
library(tidyr)

##Obtain filled GOTeDNA template and link data to OBIS specific column headers

setwd("C:/Users/HEADK/Desktop/OBIS_prep")

GOT22 <- read.xlsx("GOTeDNA-22_data.xlsx", sheet = "Sample_Metadata")
GOT22_metabar <- GOTeDNA_22

#Standardize materialSampleID values

GOT22 <- GOT22 %>%
  mutate(
    materialSampleID = as.character(materialSampleID),
    materialSampleID = sub("\\.0$", "", materialSampleID),
    materialSampleID = gsub("_", ".", materialSampleID)
  )

GOT22_metabar <- GOT22_metabar %>%
  mutate(
    materialSampleID = as.character(materialSampleID),
    materialSampleID = sub("\\.0$", "", materialSampleID),
    materialSampleID = gsub("_", ".", materialSampleID)
  )

#Filter for data contact and primer
# Filter for data contact and sample IDs #Find polygons and fix coordinates here
GOT22_f <- GOT22 %>%
  filter(
    ownerContact == "nick.jeffery@dfo-mpo.gc.ca",
    #decimalLatitude  >= 44.80 & decimalLatitude <= 45.27,   #Filter for coordinates specific to SAB (to exclude Gully and Fundian data)
    #decimalLongitude >= -68.78 & decimalLongitude <= -66.07
  )



# Filter the metabarcoding dataset
GOT22_metabar_f <- GOT22_metabar %>%
  filter(
    target_gene == "COI",
  )


#Combine files
GOT22_joined <- GOT22_f %>%
  left_join(
    GOT22_metabar_f,
    by = "materialSampleID",
    relationship = "many-to-many"
  )

#Fix the values in eventDate

GOT22_joined <- GOT22_joined %>%
  mutate(
    # Always work from a character version
    eventDate_chr = as.character(eventDate),
    
    # Try to interpret as numeric (Excel serial)
    eventDate_num = suppressWarnings(as.numeric(eventDate_chr)),
    
    # Fix 2-digit year dates like 21/8/19 -> 21/8/2019
    eventDate_chr_2y = if_else(
      str_detect(eventDate_chr, "^\\d{1,2}/\\d{1,2}/\\d{2}$"),
      sub("(\\d{1,2}/\\d{1,2}/)(\\d{2})$", "\\120\\2", eventDate_chr),
      NA_character_
    ),
    
    eventDate = case_when(
      # 1) Excel numeric serials (42293 etc.)
      !is.na(eventDate_num) ~ as.Date(eventDate_num, origin = "1899-12-30"),
      
      # 2) dd/mm/yy we just expanded to dd/mm/20yy
      !is.na(eventDate_chr_2y) ~ as.Date(eventDate_chr_2y, format = "%d/%m/%Y"),
      
      # 3) Normal dd/mm/yyyy strings
      TRUE ~ as.Date(eventDate_chr, format = "%d/%m/%Y")
    )
  ) %>%
  select(-eventDate_chr, -eventDate_num, -eventDate_chr_2y)

#Create unique event and occurrence IDs
# Create unique event and occurrence IDs for non-controls with coordinates
GOT22_joined <- GOT22_joined %>%
  # 1) Remove any rows that are controls (controlType not NA/empty)
  filter(is.na(controlType) | controlType == "") %>%
  
  # 2) Drop any rows missing coordinates
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  ) %>%
  
  # 3) Work within each original eventID.x
  group_by(eventID.x) %>%
  mutate(
    row_suffix = row_number(),
    
    # eventID for every remaining row (all have coords now)
    eventID = paste0(eventID.x, "-", row_suffix),
    
    # occurrenceID for every remaining row (present + absent)
    occurrenceID = if_else(
      is.na(occurrenceID),
      paste0("DFO-MCT-PerleyOffshore_2024", materialSampleID, "-", row_suffix),
      paste0(occurrenceID, "-", row_suffix)
    )
  ) %>%
  ungroup() %>%
  select(-row_suffix)


#Set up occurrence dataframe

map_raw <- read.csv("OBIS_Perley_2024_COI.csv", stringsAsFactors = FALSE)

map_long <- bind_rows(
  map_raw %>%
    select(Seq_Run_ID, lib_id = Lib_ID) %>%
    mutate(lib_id = gsub("_", ".", lib_id)),
  
  map_raw %>%
    select(Seq_Run_ID = Seq_Run_ID.1,
           lib_id    = Lib_ID.1) %>%
    mutate(lib_id = gsub("_", ".", lib_id))
) %>%
  filter(!is.na(lib_id) & lib_id != "")

map_clean <- map_long %>%
  filter(!str_detect(lib_id, "ENEG")) %>%
  mutate(
    # Extract SAB24.EDNA.014 from SAB24.EDNA.014.S1.MiSeq_250131
    materialSampleID = str_extract(lib_id, "^[A-Za-z0-9]+\\.EDNA\\.\\d{3}")
  ) %>%
  distinct(materialSampleID, .keep_all = TRUE)

map_clean <- map_clean %>%
  # treat blanks as NA
  mutate(Seq_Run_ID = na_if(Seq_Run_ID, "")) %>%
  # fill seq_run_id down within the whole table
  fill(Seq_Run_ID)

occurrence <- GOT22_joined %>%
  mutate(
    GOTeDNAprotocol_ID = GOTeDNA_ID.x,
    GOTeDNAprotocolVersion = protocolVersion,
    project_contact = ownerContact,
    bibliographicCitation = bibliographicCitation,
    basisOfRecord = basisOfRecord.x,
    locationID = samplingStation,
    materialSampleID = materialSampleID,
    samp_name = materialSampleID,
    technical_rep_id = "1",
    recordedBy = recordedBy.x,
    geo_loc_name = "Canada: Nova Scotia, Offshore",
    nameAccordingTo = "WoRMS",
    country = "Canada",
    datasetID = "DFO-MCT-SAB-GUL-FUN-2024",
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    minimumDepthInMeters = sampleDepth,
    maximumDepthInMeters = sampleDepth,
    totalDNAconc = totalDNAconc,
    unitsDNAconc = unitsDNAconc,
    samp_size = "3",
    samp_size_unit = "L",
    filtrationType = "peristaltic",
    size_frac = "1.2",
    filter_name = "Smith-Root",
    eventDate = eventDate,
    occurrenceID = occurrenceID,
    occurrenceStatus = case_when(
      organismQuantity > 0 ~ "present",
      organismQuantity == 0 ~ "absent",
      TRUE ~ NA_character_
    ),
    language = "en",
    filter_material = "Polyethersulfone (PES)",
    materialSampleID = materialSampleID,
    month = month(eventDate),
    bibliographicCitation = bibliographicCitation,
    year = year(eventDate),
    occurrenceID = occurrenceID,
    samp_name = materialSampleID,
    project_name = "Perley Offshore 2024",
    scientificName = scientificName,
    env_broad_scale = "marine biome (ENVO:00000447)",
    env_local_scale = "oceanic benthopelagic zone biome (ENVO:01000040)",
    env_medium = "sea water (ENVO:00002149)",
    targetTaxonomicAssay = "vertebrate",
    samp_store_sol = "100% ethanol",
    nucl_acid_ext_kit = "DNeasy® Blood & Tissue Kit (Qiagen)",
    nucl_acid_ext_modify = "For the extractions, all sample filters were cut in half and processed according to ABL’s eDNA extraction protocol. This particular protocol utilized the QIACube Connect (Qiagen) and a modified version of the Qiagen DNeasy Blood and Tissue Kit (240)",
    seq_kit = "MiSeq Reagent Kit v3 (Illumina)",
    samp_store_temp = "-80",
    samp_collect_device = "Niskin bottle",
    ref_db = "NCBI nt_euk",
    otu_seq_comp_appr = "BLAST v2.15.0",
    otu_db = "QIIME2 Amplicon v2024.10, DADA2 v2024.10.0 (in R, DADA2 v1.30.0)",
    assay_name = "COI-1 (metabarcoding)",
    assay_type = "metabarcoding",
    target_gene = target_gene,
    target_subfragment = target_subfragment,
    instrument = "Illumina MiSeq [OBI_0002003]",
    pcr_primer_forward = "GGWACWGGWTGAACWGTWTAYCCYCC",
    pcr_primer_reverse = "ACTTTCGTTCTTGATYRA",
    pcr_primer_name_forward = "mICOIintF",
    pcr_primer_name_reverse = "jgHCO2198",
    pcr_primer_reference = "https://doi.org/10.1002/ece3.4213 | https://doi.org/10.1186/1742-9994-10-34 | https://doi.org/10.1111/1755-0998.12138",
    identificationRemarks = identificationRemarks,
    checkls_ver = "1.0.2",
    sampleSizeValue = sampleSizeValue,
    sampleSizeUnit = sampleSizeUnit,
    filter_passive_active_0_1 = "1",
    associatedSequences = "In progress for NCBI",
    samp_category = "sample",
    samp_collec_device = samp_collect_device,
    project_id = Project_ID,
    pcr_0_1 = "1",
    samp_mat_process = "filtration",
    platform = "ILLUMINA",
    tax_assign_cat = "sequence similarity"
  ) %>%
  
  #filter(
  #!is.na(minimumDepthInMeters),
  #!is.na(maximumDepthInMeters),
  #!is.na(target_gene),
  #!is.na(scientificName)
  #) %>%
  left_join(map_clean, by = "materialSampleID"
  ) %>%
  mutate(seq_run_id = Seq_Run_ID)

#Look up AphiaID (taxonID) from WoRMS by ScientificName

species_lookup <- occurrence %>%
  distinct(scientificName) %>%
  mutate(
    taxonID = map_int(
      scientificName,
      ~ {
        # Keep NAs as NAs
        if (is.na(.x) || .x == "") return(NA_integer_)
        
        id <- tryCatch(
          wm_name2id(.x, marine_only = FALSE),
          error = function(e) NA
        )
        
        if (is.null(id) || length(id) == 0 || is.na(id[1])) {
          NA_integer_
        } else {
          as.integer(id[1])
        }
      }
    )
  )

# Join taxonID into occurrence, but DO NOT drop rows with no match
occurrence <- occurrence %>%
  select(-any_of(c("taxonID", "scientificNameID",
                   "kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "aphiaID"))) %>%  # we'll re-create these
  left_join(species_lookup, by = "scientificName") %>%
  mutate(
    # optional helper flag: which rows had a WoRMS match?
    worms_match = !is.na(taxonID),
    # keep AphiaID explicitly if you want it separate from taxonID
    aphiaID = taxonID
  )


# Build a taxonomy lookup for each unique non-NA AphiaID


# Helper to pull a rank from wm_classification()
.rank_val <- function(classif, rank_name) {
  if (is.null(classif)) return(NA_character_)
  if (!all(c("rank", "scientificname") %in% names(classif))) return(NA_character_)
  v <- classif$scientificname[classif$rank == rank_name]
  if (length(v) == 0) NA_character_ else as.character(v[1])
}

get_taxonomy_for_id <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  rec     <- tryCatch(wm_record(id),        error = function(e) NULL)
  
  tibble(
    taxonID           = id,
    scientificNameID  = if (!is.null(rec) && "lsid" %in% names(rec))
      rec$lsid else NA_character_,
    kingdom           = .rank_val(classif, "Kingdom"),
    phylum            = .rank_val(classif, "Phylum"),
    class             = .rank_val(classif, "Class"),
    order             = .rank_val(classif, "Order"),
    family            = .rank_val(classif, "Family"),
    genus             = .rank_val(classif, "Genus")
  )
}

unique_ids <- occurrence$taxonID %>%
  unique() %>%
  discard(is.na)

taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)


#Join WoRMS taxonomy back to occurrence

occurrence <- occurrence %>%
  select(-any_of(c("kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "scientificNameID"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")


#This code is used to extract information from otus.database.fasta.gz files (barque) to obtain the OTUids and sequences for input into the DNA-derived data extension for GBIF/OBIS

#install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)

# 1. Read the OTU database File
otu_df <- read.csv("SAB2024_COIASVs_WithSEQS.csv", stringsAsFactors = FALSE)

#DNA derived data extension

#Link DNA sequences from otu_df

dna_df <- occurrence %>%
  left_join(
    otu_df %>%
      select(
        scientificName = Species,
        seq_id = OTU.ID,
        dna_sequence = sequence
      ),
    relationship = "many-to-many"
  )

###########################################
#Run Quality Control checks

#install.packages("devtools")
#devtools::install_github("iobis/obistools")

library(obistools)
library(Hmisc)

occur <- occurrence
dna <- dna_df

#Check that the taxa names match with WoRMS
#worms <- match_taxa(unique(occur$scientificName)) #choose option 1 each time of match options

#occur <- merge(occur, worms, by="scientificName")
colnames(occur)[colnames(occur) == "scientificNameID.y"] = "scientificNameID"

#Check that all required fields are present in the occurrence table
#Check the uniqueness of the occurrenceID field (Want to = TRUE)
length(occur$occurrenceID) == length(unique(occur$occurrenceID))

#clean files
cleaned_occ <- occur %>%
  select(
    occurrenceID, bibliographicCitation, materialSampleID, eventDate, decimalLatitude, decimalLongitude, scientificName,
    organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,  associatedSequences,
    basisOfRecord, locationID,  recordedBy, country, datasetID, occurrenceStatus, minimumDepthInMeters, maximumDepthInMeters,
    language,  month, year, taxonID, scientificNameID, kingdom, phylum, class, order, family, genus
  )

cleaned_dna <- dna %>%
  select(
    project_name, dna_sequence, target_gene, pcr_primer_forward, pcr_primer_reverse, samp_name,
    env_broad_scale, env_local_scale, env_medium, samp_mat_process, size_frac, samp_size,
    samp_size_unit, otu_db, seq_kit, otu_seq_comp_appr, pcr_primer_name_forward, pcr_primer_name_reverse,
    pcr_primer_reference,  occurrenceID
  )

##Write updated + cleaned files; When confident in the data, write these files and put into OBIS
#Write CSV files for input into OBIS
write.csv(cleaned_occ, "OBIS_MCT_SABGULFUN2024_COI_occurrence.csv")
write.csv(cleaned_dna, "OBIS_MCT_SABGULFUN2024_COI_dnaderiveddata.csv")


#Create eMOF df
#I need to create code to organize the measurement or fact df -> see OBIS Manual for details

emof <- dna_df %>%
  mutate(
    occurrenceID,
    measurementType = NA_character_, #This is what I put all of the FAIRe names in
    measurementValue = NA_character_,
    measurementUnit = NA_character_,
    measurementTypeID = NA_character_, #This is where I cite all of the FAIRe names
    measurementValueID = NA_character_,
    measurementUnitID = NA_character_,
    measurementRemarks = NA_character_
  ) %>%
  select(
    occurrenceID, measurementType, measurementValue, measurementUnit, measurementTypeID, measurementValueID, measurementUnitID, measurementUnitID, measurementRemarks
  )


#I need to create code that will provide information for each FAIRe name in the measurementType column per occurrenceID
#Insert proper FAIRe names from occurrence below

url_map <- c(
  GOTeDNAprotocol_ID              = "https://github.com/GOTeDNA-OBON",
  GOTeDNAprotocolVersion          = "https://github.com/GOTeDNA-OBON",
  LClabel                  = "https://github.com/GOTeDNA-OBON",
  ownerContact             = "https://github.com/GOTeDNA-OBON",
  samplingStation          = "https://github.com/GOTeDNA-OBON",
  filtrationType           = "https://github.com/GOTeDNA-OBON",
  totalDNAconc             = "https://github.com/GOTeDNA-OBON",
  unitsDNAconc             = "https://github.com/GOTeDNA-OBON",
  dateFiltration           = "https://github.com/GOTeDNA-OBON",
  timeFiltration           = "https://github.com/GOTeDNA-OBON",
  volumeFiltered           = "https://github.com/GOTeDNA-OBON",
  depthWaterTemp           = "https://github.com/GOTeDNA-OBON",
  occurrenceID             = "https://manual.obis.org/darwin_core.html"
)

emof <- dna_df %>%
  select(
    seq_id, samp_category, checkls_ver, assay_name, assay_type, targetTaxonomicAssay,
    geo_loc_name, technical_rep_id, project_contact, seq_run_id, lib_id, project_id,
    pcr_0_1, samp_store_sol, samp_store_temp, platform, instrument, tax_assign_cat,
    LClabel, occurrenceID, nucl_acid_ext_kit, filter_material
  ) %>%
  distinct() %>%
  mutate(across(-occurrenceID, as.character)) %>%
  pivot_longer(
    cols = -occurrenceID,
    names_to = "measurementType",
    values_to = "measurementValue"
  ) %>%
  mutate(
    measurementTypeID = dplyr::recode(
      measurementType,
      !!!url_map,
      .default = "https://github.com/FAIR-eDNA/FAIRe_checklist"
    ),
    measurementUnit    = NA_character_,
    measurementValueID = NA_character_,
    measurementUnitID  = NA_character_,
    measurementRemarks = NA_character_
  )

#write csv file and put into OBIS
write.csv(emof, "OBIS_MCT_SABGULFUN2024_COI_emof.csv")
