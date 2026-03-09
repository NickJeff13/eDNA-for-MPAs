# Create taxplore plots from species abundance matrices


# load libraries ----------------------------------------------------------

library(taxplore)
library(worrms)
library(dplyr)
library(patchwork)
library(obistools)


# Set krona plots options -------------------------------------------------

set_chart_display_opts(font=15, 
                       showMagnitude = T, 
                       key = F,
                       hue_range = c(40, 180))

#set taxa ranks
tax_ranks <- c("kingdom","phylum","class","order","family","genus","species")

# setup dataframe ---------------------------------------------------------
# look up phylum etc data for each species if it's not present

df <- esi22.coi.merge.filt
df <- esi23.12s.perl.merge2
df <- esi23.coi.perl.filt2

df$species <- gsub("_"," ", df$species)
#check col names
colnames(df)


species_lookup <- df %>%
  distinct(species) %>%
  mutate(
    taxonID = map_int(
      species,
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

# Join taxonID into dataframe, but DO NOT drop rows with no match
df.aphia <- df %>%
  select(-any_of(c("taxonID", "scientificNameID",
                   "kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "aphiaID"))) %>%  # we'll re-create these
  left_join(species_lookup, by = "species") %>%
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

unique_ids <- df.aphia$taxonID %>%
  unique() %>%
  discard(is.na)

taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)

#Join WoRMS taxonomy back to occurrence

df.taxonmy <- df.aphia %>%
  select(-any_of(c("kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "scientificNameID"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")

#Check that the taxa names match with WoRMS and merge any corrections
worms <- match_taxa(unique(df.taxonmy$species)) #select "y" and choose option 1 each time

#Wait.

df.merge <- left_join(df.taxonmy, worms, by=c("species"="scientificName"))
colnames(df.merge)[colnames(df.merge) == "scientificNameID.x"] = "scientificNameID"

#Check that all required fields are present in the occurrence table
#Check the uniqueness of the occurrenceID field (Want to = TRUE)

length(df.merge$occurrenceID) == length(unique(df.merge$occurrenceID))


#try plotting krona plot
plot_krona(df.merge) #Looks weird on its own

 plot_krona(df.merge[tax_ranks], interactive = T, 
           outfile = "figures/2023_Perley/ESI/2023_PerESI_COI_invertebrates.html")

# count values per species and add that to plot
df.merge.sum <- df.merge %>% 
  dplyr::select(Species,starts_with("Sample")) %>%
  group_by(Species) %>%
  summarise(across(everything(), sum)) %>%
  mutate(n_total = rowSums(across(starts_with("Sample")))) %>%
  data.frame()

#add phylum etc back in
df.merge.sums <- left_join(df.merge.sum, df.merge %>% select(c(Species, kingdom, phylum, class, order, family, genus)), by="Species", multiple = "first", relationship = "many-to-one")

plot_krona(df.merge.sums[tax_ranks], df.merge.sums$n_total)



