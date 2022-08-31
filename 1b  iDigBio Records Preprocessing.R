###======================================================### 
###        Pre-Processing records from iDigBio           ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

#### The following script contains the previous steps for clean and filter
#### biodiversity records obtained from iDigBio database..
#### Once done move to script 1 to do other filters and pre-processing.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation
library(lubridate)  # Dates
# Environment and data ----------------------------------------------------
rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)
# Read occurrences csv and check dimensions
data <- data.table::fread("idigBio.csv", header = TRUE) 

# Reduced version of idigBio dataset -------------------------------------
# 1. Filter fields 
colnames(data)
columns <- c("coreid", "dwc:family", "dwc:genus", "dwc:specificEpithet", 
             "dwc:infraspecificEpithet", "dwc:taxonRank", "dwc:taxonomicStatus", 
             "dwc:scientificName", "gbif:canonicalName",
             "dwc:country","idigbio:isoCountryCode", "dwc:county", "dwc:locality",
             "dwc:stateProvince", "dwc:verbatimEventDate",
             "idigbio:geoPoint", "dwc:coordinateUncertaintyInMeters", 
             "idigbio:eventDate", "idigbio:dataQualityScore", 
             "dwc:basisOfRecord", "dwc:recordedBy")
data <- data[ , columns, with = FALSE]
# Rename fields related to Darwin Core Terms
data <- dplyr::rename(data, "family" = "dwc:family", 
                      'genus' = "dwc:genus", 
                      'specificEpithet' = "dwc:specificEpithet", 
                      'infraspecificEpithet' = "dwc:infraspecificEpithet", 
                      'taxonRank' = "dwc:taxonRank", 
                      'taxonomicStatus' = "dwc:taxonomicStatus", 
                      'scientificName' = "dwc:scientificName", 
                      'canonicalName' = "gbif:canonicalName",
                      'country' = "dwc:country",
                      'county' = "dwc:county", 
                      'locality' = "dwc:locality",
                      'stateProvince' = "dwc:stateProvince", 
                      'verbatimEventDate' = "dwc:verbatimEventDate",
                      'coordinateUncertaintyInMeters' = "dwc:coordinateUncertaintyInMeters",
                      'basisOfRecord' = "dwc:basisOfRecord", 
                      'recordedBy' = "dwc:recordedBy")
# Previous filters --------------------------------------------------------
# 2. Check Basis of Record field: What type of occurrences do we have?
unique(data$basisOfRecord)
# And how many in each basisOfRecord class?
data[, .N, by = basisOfRecord]
# 3. Keep records with appropriate basis of record:
data <- data[basisOfRecord != "fossilspecimen", ] 

# Create new fields from idigbio ones-----------------------------------------

# 1. Extract coordinates from 'geoPoint' field
coords <- idig %>% 
          mutate(coords = str_extract_all(`idigbio:geoPoint`, "-?\\d+([.,]\\d+)?"))
coords <- coords$coords
n <- length(coords[[1]])
coords <- structure(coords, row.names = c(NA, -n), class = "data.table")
coords <- transpose(coords)
setDT(coords)
names(coords) <- c('decimalLatitude', 'decimalLongitude')

data <- cbind(data, coords)
data$decimalLatitude  <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)
rm(coords)
# 2. Filter coordinates of study area
data <- data[decimalLatitude >= 20 & decimalLatitude < 90, ]

# 3. Extract date value
date <- as.Date(as.POSIXct(data$`idigbio:eventDate`), tz = "")
date <- data.frame(date)
date$year <- year(ymd(date$date))
date$month <- month(ymd(date$date))
date$day <- day(ymd(date$date))

data <- cbind(data, date)
# Delete non reliable dates
data <- data[year >= 1700 & year <= 2022 | is.na(year), ] 
data$date <- paste(data$day, data$month, data$year, sep = '/')

# Go to 'Taxonomic Check' Section in script 1
