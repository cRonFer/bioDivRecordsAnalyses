###======================================================### 
###     Download and Pre-Processing records from BIEN    ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

#### The following script contains the previous steps for download clean 
#### and filter biodiversity records obtained from BIEN database.
#### Once done move to script 1 to do other filters and pre-processing.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation
library(lubridate)  # Dates
#devtools::install_github("bmaitner/RBIEN",force = TRUE)
library(BIEN)
# Environment and data ----------------------------------------------------
rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)

# Download Data from BIEN -----------------------------------------------------------------------
# Load a list of family names from Briopsida to download
families <- read.csv('bryo_families.csv', header = TRUE, sep = ";")
families <- families$families

occbien <- BIEN_occurrence_family(
            families,
            new.world = NULL,
            observation.type = TRUE,
            all.taxonomy = TRUE,
            native.status = TRUE,
            natives.only = FALSE, # doubtly but yes
            political.boundaries = TRUE,
            collection.info = TRUE,
            limit = 2000000
          )
setDT(occbien)
colnames(occbien)

# Reduced version of the gbif dataset -------------------------------------
# Filter fields 
columns <- c("scrubbed_family", "family_matched",             
             "verbatim_scientific_name",              
             "name_matched", "name_matched_author",         
             "higher_plant_group", "scrubbed_taxonomic_status",
             "native_status", "native_status_sources",
             "country",	"state_province", "scrubbed_species_binomial",
             "latitude",	"longitude",	"date_collected",	"datasource",
             "dataset","recorded_by", "observation_type" )

data <- occbien[ , columns, with = FALSE]
# Rename fields related to Darwin Core Terms
setnames(data, c('family_matched', 'name_matched', 
                 'recorded_by',
                 'latitude', 'longitude',
                 'observation_type', 
                 'date_collected',
                 'scrubbed_species_binomial'), 
                c('family', 'scientificName', 
                  'recordedBy', 
                  'decimalLatitude', 'decimalLongitude', 
                  'basisOfRecord', 
                  'eventDate',
                  'species') )

# Previous filters --------------------------------------------------------
# 1. Exclude GBIF records from 'datasource'
data <- data[datasource != 'GBIF', ]
# 2. Filter records from study area
data <- data[decimalLatitude >= 20 & decimalLatitude < 90, ]
# 3. observation type
unique(data$basisOfRecord)
data[, .N, by = basisOfRecord]
# 4. Temporal processing
summary(data$eventDate)
### extract year
data$year <- year(ymd(data$eventDate))
summary(data$year)
data <- data[year >= 1800, ] # Delete unreliable years
data$month <- month(data$eventDate)
data$day <- day(data$eventDate)
data$date <- paste(data$day, data$month, data$year, sep = '/')



