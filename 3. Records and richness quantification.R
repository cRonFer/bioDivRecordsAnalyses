###======================================================### 
###                 Quantitative analysis                ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

#### The following script contains the code necessary to
#### the quantification of records and observed number of species by cell 
#### in the studied cases
#### Adjust the parameters to perform your own analysis based on different filters.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation
library(CoordinateCleaner)  # data cleaning functions
# Environment and data ----------------------------------------------------

rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)

data <- readRDS('dataset')  # Dataset from script 2
# Filter dataset to only records with consensus species name:
data_cons <- data[!is.na(sppName), ]

### Count all cases in dataset by cell = 'id' ####
n <- data[ , .(records = .N), by = c('id','region')]

# Consensus species name number of records
n <- data_cons[, .(records3 = .N), by = id][n, on = .(id = id)]
# Consensus species name number of observed species
n <- data_cons[, .(richness3 = uniqueN(sppName)), by = id][n, on = .(id = id)]

### Filter Preserved specimens records ####
# Number of records
n <- data_cons[basisOfRecordRev =='PRESERVED_SPECIMEN', ][ , .(recordsPres3 = .N), by = id][n, on = .(id = id)]
# Number of observed species
n <- data_cons[basisOfRecordRev =='PRESERVED_SPECIMEN', ][ , .(richnessPres3 = uniqueN(sppName)), by = id][n, on = .(id = id)]

n[is.na(n)] <- 0 # assign 0 to cells with no data (NA)

# Proportional records
n$prop_recPres <- 1 - n$recordsPres3 / n$records3
# Proportional observed species
n$prop_richPres <- 1 - n$richnessPres3 / n$richness3

### Temporal coverage counts by cell  #### 
old <- data[, .SD[which.min(year)]$year]  # extract oldest year with data
temp <- c(old, 1900, 1970, 1990)  # list temporal ranges of study

# Loop to calculate number of records and observed species in each temporal range
for (i in temp){
  count <- data_cons[year >= i, ][, .(records = .N), by = id]
  setnames(count, 'records', paste0('rec_', i))
  n <- count[n, on = .(id = id)]
  richness <- data_cons[year >= i, ][, .(richness = uniqueN(sppName)), by = id]
  setnames(richness, 'richness', paste0('rich_', i))
  n <- richness[n, on = .(id = id)]
}
# Same for temporal range 1970-2000
n <- data_cons[year >= 1970 & year <=2000, ][, .(rec_7000 = .N), by = id][n, on = .(id = id)]
n <- data_cons[year >= 1970 & year <=2000, ][, .(rich_7000 = uniqueN(sppName)), by = id][n, on = .(id = id)]

n[is.na(n)] <- 0 # assign 0 to cells with no data (NA)

# Proportional records
n$p_rec1600 <- 1 - n$rec_1600 / n$records3
n$p_rec1900 <- 1 - n$rec_1900 / n$records3
n$p_rec1970 <- 1 - n$rec_1970 / n$records3
n$p_rec1990 <- 1 - n$rec_1990 / n$records3
n$p_rec7000 <- 1 - n$rec_7000 / n$records3
# Proportional observed species
n$p_rich1600 <- 1 - n$rich_1600 / n$richness3
n$p_rich1900 <- 1 - n$rich_1900 / n$richness3
n$p_rich1970 <- 1 - n$rich_1970 / n$richness3
n$p_rich1990 <- 1 - n$rich_1990 / n$richness3
n$p_rich7000 <- 1 - n$rich_7000 / n$richness3

### Taxonomic counts by cell ####

#### *GBIF*
unique(data$gbif_status)

# Number of observed species
n <- data[gbif_status == 'ACCEPTED' |
          gbif_status == 'SYNONYM' , ][, .(gbif_rich = uniqueN(gbif_species)), by = id][n, on = .(id = id)]

#### *TPL*
unique(data$TPL_taxonStatus)
# Number of observed species
n <- data[TPL_taxonStatus == 'accepted', ][, .(tpl_rich = uniqueN(TPL_speciesName)), by = id][n, on = .(id = id)]

#### *TNRS*
unique(data$tnrs_status)
# Number of observed species
n <- data[tnrs_status == 'Accepted' | 
          tnrs_status == 'Synonym', ][, .(tnrs_rich = uniqueN(tnrs_accepted_species)), by = id][n, on = .(id = id)]

#### *TROPICOS*
unique(data$Tropicos_taxonomic_status)
# Number of observed species
n <- data[Tropicos_taxonomic_status == 'accepted' |
          Tropicos_taxonomic_status == 'synonym', ][, .(trop_rich = uniqueN(Tropicos_Accepted_species)), by = id][n, on = .(id = id)]

#### *WFO*
unique(data$wfo_taxonomicStatus)
# Number of observed species
n <- data[wfo_taxonomicStatus == 'accepted', ][, .(WFO_rich = uniqueN(wfo_acceptedName)), by = id][n, on = .(id = id)]

# Taxonomic variations of each source against taxon. consensus observed species
n[is.na(n)] <- 0 # assign 0 to cells with no data (NA)

n$p_richGBIF <- (n$gbif_rich - n$richness3) / n$richness3
n$p_richTNRS <- (n$tnrs_rich - n$richness3) / n$richness3
n$p_richTPL <- (n$tpl_rich - n$richness3) / n$richness3
n$p_richTROP <- (n$trop_rich - n$richness3) / n$richness3
n$p_richWFO <- (n$WFO_rich - n$richness3) / n$richness3

# -1 corresponds to complete lost of richness
# 0 when richness is same using consensus checklist or a particular standardization
# 1 doubled richness compare to consensus list

### Duplicate counts by cell ####
# Species Name + Coordinates + date of collection + recorder names
data_duplic <- cc_dupl(data_cons,
                      lon = "decimalLongitude",
                      lat = "decimalLatitude",
                      species = "sppName",
                      additions = c('recBy_Rev', 'date'))
# Number of records
n <- data_duplic[ , .(dupl_date_rec = .N), by = id][n, on = .(id = id)]

### Species Name + Coordinates + date of collection
data_duplic <- cc_dupl(data_cons,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "sppName",
                       additions = 'date')
# Number of records
n <- data_duplic[ , .(dupl_dateRecBy_rec = .N), by = id][n, on = .(id = id)]

n[is.na(n)] <- 0 # assign 0 to cells with no data (NA)

# Proportional records from consensus list
n$p_date <- 1 - n$dupl_date_rec / n$records3
# Proportional records from consensus list
n$p_recByDate <- 1 - n$dupl_dateRecBy_rec / n$records3 

# Save output
fwrite(n, 'quantification.csv', sep=";")

