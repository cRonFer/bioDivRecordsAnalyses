###======================================================### 
###         Taxonomical standardization of records       ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

#### The following script contains the code to standardize the scientific
#### names of Bryopsida.
#### The checklist was obtained in script 1.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation

library(rgbif)  # tax. standardization
library(Taxonstand)  # tax. standardization
library(WorldFlora)  # tax. standardization

# Environment and data ----------------------------------------------------

rm(list = ls(all.names = TRUE))
wd <- ''  # Write your working directory
setwd(wd)
 
# Read checklist csv and check dimensions
checklist <- data.table::fread()   
dim(checklist)
checklist$scientificName <- tolower(checklist$scientificName)
# Standardization by source:
#### *TNRS* Taxonomic Name Resolution Service ####
# Due to errors with R package for TNRS tax. standardization
# We use the web API and load the output obtained:
tnrs_stand <- fread("tnrs_result.tsv", sep = "\t", header = TRUE) # All matches

# Select necessary fields, rename and create new ones
tnrs_stand <- tnrs_stand[,c('Name_submitted', 'Overall_score',
                            'Taxonomic_status',
                            'Accepted_name','Accepted_species',
                            'Accepted_name_author', 'Overall_score_order')]
tnrs_stand <- tnrs_stand[Overall_score_order == 1,]  # Best name of every source
tnrs_stand$Overall_score_order <- NULL  # Field not necessary now

tnrs_stand <- dplyr::rename(tnrs_stand, 'tnrs_score' = 'Overall_score', 
                     'tnrs_status' = 'Taxonomic_status',
                     'tnrs_accepted_name' = 'Accepted_name',
                     'tnrs_accepted_species' = 'Accepted_species',
                     'tnrs_author_name' = 'Accepted_name_author')

checklist <- merge(checklist, tnrs_stand,  # Bind to checklist
                   by.x = "scientificName", by.y = "Name_submitted", 
                   all.x = TRUE) 

#### *Tropicos* ####
# Due to errors with 'TNRS' R package for tax. standardization
# We use the web API, selecting 'Tropicos' option and load the output obtained:
trop_stand <- fread("trop_result.tsv", sep = "\t", header = TRUE) # Best matches only
# Select necessary fields, rename and create new ones
trop_stand <- trop_stand[,c('Name_submitted', 'Overall_score',
                            'Taxonomic_status',
                            'Accepted_name','Accepted_species',
                            'Accepted_name_author')]

trop_stand <- dplyr::rename(trop_stand, 'Tropicos_score' = 'Overall_score', 
                      'Tropicos_taxonomic_status' = 'Taxonomic_status',
                      'Tropicos_Accepted_name' = 'Accepted_name',
                      'Tropicos_Accepted_species' = 'Accepted_species',
                      'Tropicos_Accepted_name_author' = 'Accepted_name_author')

checklist <- merge(checklist, trop_stand,  # Bind to checklist
                   by.x = "scientificName", by.y = "Name_submitted", 
                   all.x = TRUE) 

#### *TPL - The Plant List; TaxonStand* #### 
tpl_stand <- TPL(checklist$scientificName, author = TRUE)

# Select Submitted Taxon and New. fields
colnames(tpl_stand)
tpl_stand <- tpl_stand[ , c('Taxon', 'New.Genus',  
                            'New.Hybrid.marker', 'New.Species',
                            'New.Infraspecific.rank', 'New.Infraspecific',
                            'New.Authority', 'New.Taxonomic.status')] 

tpl_stand$speciesName <- paste(tpl_stand$New.Genus,
                               tpl_stand$New.Hybrid.marker,
                               tpl_stand$New.Species)

tpl_stand$speciesName <- str_trim(tpl_stand$speciesName)
tpl_stand$speciesName <- str_squish(tpl_stand$speciesName) # Delete white spaces

# Reduce necessary fields and rename them:
tpl_stand <- tpl_stand[ , c('Taxon','New.Taxonomic.status', 'speciesName')]
tpl_stand <- dplyr::rename(tpl_stand, 'TPL_taxonStatus' = 'New.Taxonomic.status',
                               'TPL_speciesName' = 'speciesName') 
checklist <- merge(checklist, tpl_stand, 
                   by.x = 'scientificName', by.y = 'Taxon', 
                   all.x = TRUE)

#### *GBIF tax backbone* ####
list <- checklist$scientificName
gbif_stand <- data.frame()

for (i in list){
  gbif_check <- name_backbone(i)
  gbif_stand <- dplyr::bind_rows(gbif_stand, gbif_check)
}

gbif_stand <- gbif_stand[ , c('status', 'matchType', 'species')]

gbif_stand <- dplyr::rename(gbif_stand, 'gbif_status' = 'status',
                                  'gbif_matchType' = 'matchType',
                                  'gbif_species' = 'species')

checklist <- cbind(checklist, gbif_stand)

#### *WFO World Flora Online* ####
# Download and read wfo checklist
WFO.download(WFO.url = "http://104.198.143.165/files/WFO_Backbone/_WFOCompleteBackbone/WFO_Backbone.zip",
            WFO.remember = TRUE)
# WFO.remember(WFO.file = 'WFO_Backbone/classification.txt', WFO.data = "WFO.data", WFO.pos = 1)
WFO.data$scientificName <- tolower(WFO.data$scientificName)

# Join directly to WFO backbone
wfo_list <- checklist$scientificName
wfo_stand <- data.frame()  # Empty dataframe for the loop

for (i in wfo_list) {
    WFO.data1 <- WFO.data[scientificName == i, ]
    wfo_stand <- rbind(wfo_stand, WFO.data1)
}

# Check and exclude names at genus level
wfo_stand <- wfo_stand[!taxonRank == 'GENUS', ]

# Select basic fields:
wfo_stand <- wfo_stand[ ,c("taxonID", "scientificName",
                           "genus", "specificEpithet",
                           "acceptedNameUsageID", "taxonomicStatus")]
# Create species_name field:
wfo_stand$wfo.species <- paste0(wfo_stand$genus, " ", wfo_stand$specificEpithet)

# Create separate df's for accepted taxa and unchecked names:
df_accepted <- wfo_stand[taxonomicStatus == 'Accepted', ][ ,c('scientificName','taxonomicStatus', 'wfo.species')]
setnames(df_accepted, 'taxonomicStatus', 'wfo_status_acc')
df_unchecked <- wfo_stand[taxonomicStatus == 'Unchecked' , ][,c('scientificName','taxonomicStatus')]
setnames(df_unchecked, 'taxonomicStatus', 'wfo__status_unc')

# Join to original list by type of taxonom. status
wfo_list <- data.frame(wfo_list)
wfo_list <- merge(wfo_list, df_accepted, by.x = 'wfo_list', by.y = 'scientificName', all.x = TRUE) 
wfo_list <- merge(wfo_list, df_unchecked, by.x = 'wfo_list', by.y = 'scientificName', all.x = TRUE) 

# Create separate df's for synonym taxa:
df_synonym <- wfo_stand[taxonomicStatus == 'Synonym', ]

# Retrieve the id for synonyms creating a vector to find what is the name:
df_synonym_id <- df_synonym[ ,c('scientificName','acceptedNameUsageID')]
df_synonym_id_v <- df_synonym_id$acceptedNameUsageID

# Find accepted name based on previous ID:
df <- data.frame()

for (i in df_synonym_id_v) {
  WFO.data2 <- WFO.data[taxonID == i, ]
  WFO.data2 <- WFO.data2[ ,c('taxonID', 'scientificName','taxonomicStatus')]
  setnames(WFO.data2, "scientificName", "wfo_acceptedSyn")
  df <- rbind(df, WFO.data2)
}

df_synonym_id <- cbind(df_synonym_id, df)
df_synonym_id <- df_synonym_id[ ,c('scientificName','wfo_acceptedSyn','taxonomicStatus')]

wfo_list <- merge(wfo_list, df_synonym_id, by.x = 'wfo_list', by.y = 'scientificName', all.x = TRUE) 
setnames(wfo_list, 'taxonomicStatus', 'wfo_status_syn')

# Unify fields: species name for wfo and tax. status
wfo_list[is.na(wfo_list)] <- ""
wfo_list$wfo_acceptedName <- paste0(wfo_list$wfo_acceptedSyn, wfo_list$wfo.species)
wfo_list$wfo_acceptedName <- str_trim(wfo_list$wfo_acceptedName)
wfo_list$wfo_acceptedName <- tolower(wfo_list$wfo_acceptedName)

wfo_list$wfo_taxonomicStatus <- paste0(wfo_list$wfo_status_acc, wfo_list$wfo_status_unc, wfo_list$wfo_status_syn)
wfo_list$wfo_taxonomicStatus <- str_trim(wfo_list$wfo_taxonomicStatus)

# Join to original list of taxonomic standardization
checklist <- merge(checklist, wfo_list, by.x = 'scientificName', by.y = 'wfo_list', all.x = TRUE)

#### Filter only Accepted names ####
# GBIF Backbone
taxonom_GBIF <- checklist[gbif_status == 'ACCEPTED' | gbif_status == 'SYNONYM', 
                          ][,c('scientificName',"gbif_species")]

# TPL
taxonom_tpl <- checklist[TPL_taxonStatus == 'Accepted', ][,c('scientificName',"TPL_speciesName")]

# TNRS 
taxonom_tnrs <- checklist[tnrs_status == 'Accepted' | 
                tnrs_status == 'Synonym', ][,c('scientificName',"tnrs_accepted_species")]

# TROPICOS
taxonom_tropicos <- checklist[Tropicos_taxonomic_status == 'Accepted' |
                 Tropicos_taxonomic_status == 'Synonym', ][ ,c('scientificName', 'Tropicos_Accepted_species')]

# WFO
taxonom_wfo <- checklist[wfo_taxonomicStatus =='Accepted', ][ , c('scientificName', 'wfo_acceptedName')]

checklist2 <- merge(taxonom_GBIF, taxonom_tpl, all =TRUE, by = 'scientificName')
checklist2 <- merge(checklist2, taxonom_tnrs, all =TRUE, by = 'scientificName')
checklist2 <- merge(checklist2, taxonom_tropicos, all =TRUE, by = 'scientificName')
checklist2 <- merge(checklist2, taxonom_wfo, all =TRUE, by = 'scientificName')

# Create a consensus list of species names ----

# Transform species names to lower case and delete any extra white space
checklist2 <- checklist2 %>%
  mutate(across(c(gbif_species, TPL_speciesName, tnrs_accepted_species,
                  Tropicos_Accepted_species, wfo_acceptedName), tolower))
checklist2 <- checklist2 %>%
  mutate(across(c(gbif_species, TPL_speciesName, tnrs_accepted_species,
                  Tropicos_Accepted_species, wfo_acceptedName), str_squish))

checklist2 <- checklist2 %>%
  mutate(across(c(gbif_species, TPL_speciesName, tnrs_accepted_species,
                  Tropicos_Accepted_species, wfo_acceptedName), str_trim))

# Count whether same species name exists compared to the rest
checklist2$count2 <- apply(checklist2, 1, function(x) sum(x[c('TPL_speciesName', 'tnrs_accepted_species','Tropicos_Accepted_species','wfo_acceptedName')] %in% x['gbif_species']))
checklist2$count3 <- apply(checklist2, 1, function(x) sum(x[c('gbif_species', 'tnrs_accepted_species','Tropicos_Accepted_species','wfo_acceptedName')] %in% x['TPL_speciesName']))
checklist2$count4 <- apply(checklist2, 1, function(x) sum(x[c('gbif_species', 'TPL_speciesName','Tropicos_Accepted_species','wfo_acceptedName')] %in% x['tnrs_accepted_species']))
checklist2$count5 <- apply(checklist2, 1, function(x) sum(x[c('gbif_species', 'TPL_speciesName','tnrs_accepted_species','wfo_acceptedName')] %in% x['Tropicos_Accepted_species']))
checklist2$count6 <- apply(checklist2, 1, function(x) sum(x[c('gbif_species', 'TPL_speciesName','tnrs_accepted_species','Tropicos_Accepted_species')] %in% x['wfo_acceptedName']))

# Create new fields 'sppNameX' to establish accepted species name when 
# Filter by = 2, means 3 out of 5 sources give same name
data3 <- checklist2 %>% 
          filter(conteo2 == 2|conteo3 == 2|conteo4 == 2|conteo5 == 2|conteo6 == 2) %>% 
          mutate(sppName3 = case_when(count2 == 2 ~ gbif_species,
                            (count2 != 2 & count3 == 2) ~ TPL_speciesName,
                            (count2 != 2 & count3 != 2) ~ tnrs_accepted_species))

# Filter by =3, means 4 out of 5 sources give same name
data4 <- checklist2 %>% 
            filter(conteo2 == 3|conteo3 == 3|conteo == 3|conteo5 == 3|conteo6 == 3) %>%  
            mutate(sppName4 = case_when(conteo2 == 3 ~ gbif_species,
                             (conteo2 != 3 & conteo3 == 3) ~ TPL_species))

# Filter by =4, means all sources give same name
data5 <- checklist2 %>% 
            filter(conteo2 == 4|conteo3 == 4|conteo4 == 4|conteo5 == 4|conteo6 == 4) %>%  
             mutate(sppName5 = case_when(conteo2 == 4 ~ gbif_species))

# No consensus name data
data_nc <- checklist2 %>% 
        filter(conteo2 < 2 & conteo3 < 2 & conteo4 < 2 & conteo5 < 2 & conteo6 < 2)

l = list(data3, data4, data5, data_nc)
checklist2 <- rbindlist(l, use.names = TRUE, fill = TRUE)

# Create a new field 'sppName' to establish consensus species name:
checklist2$sppName <- paste0(sppName3, sppName4, sppName5)

# Join to records and save the output ----
data <- data_correct[checklist2, on = .(scientificName = scientificName)]

saveRDS(data, 'dataset')
