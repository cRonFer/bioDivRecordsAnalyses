###======================================================### 
###           Inventory Completeness analyses            ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

### In the following script you will find detailed steps to calculate 
### inventory completeness values for records by cell
### using KnowBR package and perform in each scenario

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation
library(KnowBR) # completeness analysis
library(rgdal) # gis
library(CoordinateCleaner) # detect duplicates

# Environment and data ----------------------------------------------------
wd <- "" # Set root working directory
rm(list = ls(all.names = TRUE))
setwd(wd)

data(adworld) # knowBR needs add world polygon to work 
# Load study area grid shapefile reprojected to CRS WGS84
grid <- readOGR("gridWgs84Data.shp")
# Load records
data <- readRDS('finalDT')  # Dataset from script 2
# Filter dataset to only records with consensus species name:
data <- data[!is.na(sppName), ]
# Add field of abundance for knowBR
data$abundance <- 1

# Function to calculate inventory completeness using KnowBR ####
knowbr <- funtion(d = d){
  KnowBPolygon(data = d, shape = grid, admAreas = FALSE,  # Use predefined grid as personalized polygons  
               shapenames = "id",  # unique "id" cell from grid
               inLon = -180, maxLon = 180, minLat = 20, maxLat = 90, # study area extent
               jpg = TRUE, dec = ".")
}

# Filtering Preserved Specimen records####
dir.create('presSp')
setwd("presSp") 
d <- data[basisOfRecord =='PRESERVED_SPECIMEN'| basisOfRecord =='preservedspecimen'|
          basisOfRecord =='Preserved Specimen'| basisOfRecord =='specimen'|
          basisOfRecord =='PreservedSpecimen'| basisOfRecord =='Preservedspecimen'|
          basisOfRecord == 'in specimen envelope',][, c('sppName','decimalLongitude','decimalLatitude','abundance')]
knowbr()

# Filtering by Temporal Coverage different ranges ####
temp <- c(1600, 1990, 1970, 1900)

for (i in temp){
  setwd(d) # set root working directory 
  dir.create(as.character(i))
  setwd(as.character(i))
  d <- data[year >= i, ][, c('sppName', 'decimalLongitude', 'decimalLatitude','abundance')]
  knowbr()
}
# 1970-2000 range 
setwd(wd) 
dir.create('70_00')
setwd('70_00')
d <- data[year >= 1970 & year <= 2000, ][ , c('sppName','decimalLongitude','decimalLatitude','abundance')]
knowbr()

# Discarding Duplicate records ####
data_correct <- data[, c('id','sppName', 'decimalLongitude', 'decimalLatitude',
                         'date','recBy_Rev','abundance')]

# A. Species Name + Coordinates + date of collection + recorder's name
data_dupl <- cc_dupl(data_correct,
                     lon = "decimalLongitude",
                     lat = "decimalLatitude",
                     species = "sppName",
                     additions = c('recBy_Rev', 'date'))
d <- data_dupl [, c('sppName', 'decimalLongitude', 'decimalLatitude','abundance')]

setwd(wd) 
dir.create('DuplRecByDate')
setwd('DuplRecByDate')
knowbr()

# B. Species Name + Coordinates + date of collection 
dataDupl <- cc_dupl(data_correct,
                         lon = "decimalLongitude",
                         lat = "decimalLatitude",
                         species = "sppName",
                         additions = 'date')
d <- dataDupl[, c('sppName', 'decimalLongitude', 'decimalLatitude','abundance')]

setwd(wd) 
dir.create('DuplDate')
setwd('DuplDate')
knowbr()

# Combination of filters - Section 4 ####
# A. Records collected after 1970 and delete duplicates by date of collection
d <- data[year >= 1970, ][, c('id','sppName', 
                            'decimalLongitude', 'decimalLatitude',
                             'date','abundance')]
d <- cc_dupl(d, lon = "decimalLongitude", lat = "decimalLatitude",
                     species = "sppName", additions = 'date')
d <- d[, c('sppName', 'decimalLongitude', 'decimalLatitude','abundance')]

setwd(wd) 
dir.create('dupl1970')
setwd('dupl1970')
knowbr()

# B. Filtering Preserved specimens without duplicates by date of collection
d <- data[basisOfRecord =='PRESERVED_SPECIMEN' | basisOfRecord =='preservedspecimen'|
          basisOfRecord =='Preserved Specimen'| basisOfRecord == 'in specimen envelope'|
          basisOfRecord =='PreservedSpecimen'| basisOfRecord =='Preservedspecimen'|
          basisOfRecord =='specimen', ][, c('id','sppName',  'date','abundance',
                                        'decimalLongitude', 'decimalLatitude')]

d <- cc_dupl(d, lon = "decimalLongitude", lat = "decimalLatitude",
                     species = "sppName", additions = 'date')
d <- d[, c('sppName', 'decimalLongitude', 'decimalLatitude','abundance')]

setwd(wd)
dir.create('duplPres')
setwd('duplPres')
knowbr()

# C. Using TNRS taxonomic standardization and discarding duplicates by date of collection
data <- bryo[tnrs_accepted_species != "", ]
data$abundance <- 1
d <- data[year >=1970, ][, c('id','tnrs_accepted_species', 
                            'decimalLongitude', 'decimalLatitude',
                            'date','abundance')]
d <- cc_dupl (d, lon = "decimalLongitude", lat = "decimalLatitude",
                     species = "tnrs_accepted_species", additions = 'date')
d <- d[, c('tnrs_accepted_species', 'decimalLongitude', 'decimalLatitude','abundance')]

setwd(wd)
dir.create('tnrsDupl')
setwd('tnrsDupl')
knowbr()
