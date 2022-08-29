###======================================================### 
###          Pre-Processing occurrences records          ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

#### The following script contains basic steps for clean and filter
#### biodiversity records after the downloading process.
#### Adjust the parameters to perform your analysis.
#### Field names have to be adapted to the dataset type,
#### here it's an example of GBIF data using Darwin Core terms.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse)  # data manipulation
library(GADMTools)  # administrative units and GIS
library(CoordinateCleaner)  # data cleaning functions
library(lubridate) # dates
library(countrycode)  # country names standardization
# Environment and data ----------------------------------------------------

rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)

# Read occurrences csv and check dimensions
data <- data.table::fread(".csv", encoding ='UTF-8', sep="\t") 

dim(data)
head(data)

# Reduced version of the gbif dataset -------------------------------------
# Filter fields 
colnames(data)
fields <- c("gbifID", "family", "genus", "species", "infraspecificEpithet", 
             "taxonRank", "scientificName",
             "countryCode", "locality", "stateProvince", 
             "decimalLatitude",	"decimalLongitude", "coordinatePrecision", 
             "day", "month", "year", 
             "occurrenceStatus", "basisOfRecord", "recordedBy", "issue")
data <- data[ , fields, with = FALSE]
colnames(data)
rm(fields)
data$date <- paste(data$day, data$month, data$year, sep = '/')
# Previous filters --------------------------------------------------------

# 1. Check and filter only presences (discard absences)
data <- data[occurrenceStatus == 'PRESENT', ]
# 2. Check Basis of Record field
# 3. What type of occurrences do we have?
unique(data$basisOfRecord)
# And how many in each basisOfRecord class?
data[, .N, by = basisOfRecord]
# 4. Remove records without appropriate basis of record:
data <- data[!basisOfRecord == 'FOSSIL_SPECIMEN', ]  # ! means exclude

####------------------ Taxonomic check ------------------#####

# 1. Is the record identified at a proper taxonomic rank?
# Let's see what do we have in our dataset:
unique(data$taxonRank)
# And how many of each taxonRank:
data[, .N, by = taxonRank]
# 2. Filter those records with appropriate taxonomic rank
data <- data[taxonRank =='SUBSPECIES' |
             taxonRank == 'VARIETY'|
             taxonRank == 'FORM'|
             taxonRank =='SPECIES', ]

# 3. Create a checklist with scientific names
checklist <- data[ ,c('family', 'genus', 'species', 'infraspecificEpithet', 
                      'taxonRank', 'scientificName') ]
checklist <- data.frame(unique(checklist)) # remove duplicated entries

write.table(checklist, 'checklist.csv', sep="\t", row.names=FALSE) 

### Go to Script 2 for taxonomic standardization of this checklist

####------------------ Geographical check -------------------------#####
# 1. Check coordinates values: ----
# Discard records with latitude/longitude values equals to zero, exact same value or NULL
data <- data[decimalLatitude!=0 & decimalLongitude!=0, ][decimalLatitude!=decimalLongitude, ]

# 2. Check coordinates precision: ----
# Here decimal digits of coordinates as a measure of precision 

## Function to count number of decimals:
decimal.num.count <- function(x) {
    stopifnot(class(x) == "character")
    x <- ifelse(grepl("\\.", x), sub(".*\\.(.*)", "\\1", x), "")
    nchar(x)
}
# Apply function to latitude and longitude fields
latDigits <- decimal.num.count(as.character(data$decimalLatitude))
data <- cbind(data, latDigits)
lonDigits <- decimal.num.count(as.character(data$decimalLongitude))
data <- cbind(data, lonDigits)

# Filter coordinates with 1 or more digit places
data <- data[latDigits > 0 & lonDigits > 0, ] 

# 3. Check records that don't meet location criteria ----
# Label coordinates placed in centroids of the country
cap <- cc_cap (data,
               lon = "decimalLongitude",
               lat = "decimalLatitude",
               value = "flagged")
data <- cbind(data, cap)

# Label coordinates placed in gbif headquarters
gbif <- cc_gbif(data,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                value = "flagged")
data <- cbind(data, gbif)

# Label coordinates from museums, gardens, institutions, zoo's... 
inst <- cc_inst(data,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                value = "flagged")
data <- cbind(data, inst)

# Filter and exclude centroids:
data <- data[cap == 'TRUE', ]
data <- data[inst == 'TRUE', ]
data <- data[gbif == 'TRUE', ]

data <- data[, 1:23]

# 4. Check if coordinates are placed in the assigned country ---- 
# Point in polygon by country analysis

# Import gpkg of the world with Adm. Units borders at country level (GADM):
countries <- read_sf('gis_layers/gadm36_levels.gpkg')
# Set the correct projection
st_crs(countries) <- "+proj=longlat +datum=WGS84 +no_defs"
countries$GID_0 <- NULL  # keep only 'Name_0' (country names of GADM)
# Transform occurrences into spatial points and project:
data$x <- data$decimalLongitude
data$y <- data$decimalLatitude

datapoints  <- st_as_sf(x = data,                         
                   coords = c("x", "y"),
                   crs = "+proj=longlat +datum=WGS84 +no_defs")

# Point in polygon: join each point to a polygon based on position:
data <- st_join(datapoints, countries)

setnames(data,"NAME_0", "GADM_Location")  # rename new field

# Translate 'countryCode' information (ISO 3166-1-alpha-2 = "iso2c") 
# into country names of GADM 'countryCode'
data$countryName <- countrycode(data$countryCode, 
                                origin = "iso2c", 
                                destination = "country.name")

# Label those that didn't fall in any country polygon as 'SEA'
data$GADM_Location[is.na(data$GADM_Location)] <- "SEA"
# Check the names obtained
unique(data$countryName)
unique(data$GADM_Location)
# Make changes to match countryName and GADM-Location
data$countryName[data$countryName == "Czechia"] <- "Czech Republic" # Change 'Czechia'='Czech Republic'
data$countryName[data$countryName == "Svalbard & Jan Mayen"] <- 'Svalbard and Jan Mayen'
data$countryName[data$countryName == "Eswatini"] <- "Swaziland"
data$countryName[data$countryName == "St. Pierre & Miquelon"] <- "Saint Pierre and Miquelon"
data$countryName[data$countryName == "Turks & Caicos Islands"] <- "Turks and Caicos Islands"
data$countryName[data$countryName == "Åland Islands"] <- "Finland"
data$countryName[data$countryName == "Palestinian Territories"] <- "Israel"
data$countryName[data$countryName == "Myanmar (Burma)"] <- "Myanmar"
data$countryName[data$countryName == "North Macedonia"] <- "Macedonia"
data$countryName[data$countryName == "Hong Kong SAR China"] <- "Hong Kong"
data$countryName[data$countryName == "Bosnia & Herzegovina"] <- "Bosnia and Herzegovina"

data$GADM_Location[data$GADM_Location == "Akrotiri and Dhekelia"] <- "Cyprus"
data$GADM_Location[data$GADM_Location == "Åland"] <- "Finland"
data$GADM_Location[data$GADM_Location == "Palestina"] <- "Israel"
data$GADM_Location[data$GADM_Location == "Northern Cyprus"] <- "Cyprus"

# Match country names and label as 'FALSE' errors of location
data <- data %>% mutate(countryCheck = case_when(
                        GADM_Location != data$countryName ~ FALSE,
                        TRUE ~ TRUE))
# Subset and extract records located in country assigned by collector ('correct')
setDT(data)
data_correct <- data[countryCheck == 'TRUE', ]

# 5. Check coordinates placed in incorrect habitat (sea-land) ----

# Filter occurrences that didn't fall in the correct country polygon
data_1 <- data[countryCheck == 'FALSE', ]

# Subset and check points that fall out of their habitat GADM_Location =='SEA'
data_sea_points <- data_1[GADM_Location == 'SEA', ]

# Transform them into spatialPoints
data_sea_points$x <- data_sea_points$decimalLongitude
data_sea_points$y <- data_sea_points$decimalLatitude

datapoints  <- st_as_sf(x = data_sea_points,                         
                        coords = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs")

# Label points that fall in a fixed distance from the coastline 
# Load buffer shapefile (here 0.1 degrees - ca. 10 km)
landBuff <- read_sf('buffer_0.1.shp')
st_crs(landBuff) <- "+proj=longlat +datum=WGS84 +no_defs"
sf::sf_use_s2(FALSE)  # Spherical geometry (s2) switched off

# Function to overlap points with buffer
sea_0.1 <- st_join(datapoints, landBuff)
sea_0.1 <- unique(sea_0.1)
# If the occurrence fall in our buffer  
# (= Coastline) occurrences out of the buffer have NA
setDT(sea_0.1) # transform df to data table
sea_0.1$scalerank <- NULL
setnames(sea_0.1, "featurecla", "shoreLine0.1")

# Subset those points that fall in our buffer
dataSeaPointsGood <- sea_0.1[shoreLine0.1 == 'Coastline', ]
# 'Join attributes by nearest'
# Assign the nearest country to each of these 'coastal' occurrences 
dataSeaPointsGood$x <- dataSeaPointsGood$decimalLongitude
dataSeaPointsGood$y <- dataSeaPointsGood$decimalLatitude

datapoints  <- st_as_sf(x = dataSeaPointsGood,                         
                        coords = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs")
datapoints <- st_join(datapoints, countries, join = st_nearest_feature)

setnames(datapoints, "NAME_0", "Coastal_Country")

# Check where are located the coastal records
unique(datapoints$Coastal_Country)

# Transform some country names:
datapoints$Coastal_Country[datapoints$Coastal_Country == "Åland"] <- "Finland"
datapoints$Coastal_Country[datapoints$Coastal_Country == "St. Pierre and Miquelon"] <- "St. Pierre & Miquelon"
datapoints$Coastal_Country[datapoints$Coastal_Country == "Turks and Caicos Islands"] <- "Turks & Caicos Islands"
datapoints$Coastal_Country[datapoints$Coastal_Country == "Isle of Man"] <- "United Kingdom"
datapoints$Coastal_Country[datapoints$Coastal_Country == "Gibraltar"] <- "Spain"

# Check if they are placed in the assigned country and subset those who are correctly located
datapoints$country_shore_check <- datapoints$Coastal_Country == datapoints$countryName
setDT(datapoints)
# Filter correct records:
data_sea_correct <- datapoints[country_shore_check == 'TRUE', ]
data_sea_correct$country_shore_check <- NULL
data_sea_correct$geometry <- NULL

# Merge with previous correct dataset
data_correct <- rbind(data_correct, data_sea_correct, fill = TRUE)
data_correct <- data_correct[, 1:23]  # Discard fields not needed
# Merge datasets with incorrect records:
# Subset points that fall in wrong country and check manually:
incorrect_points <- data_1[!GADM_Location == 'SEA', ]  
# If they are correct add to data_correct dataset

# 6. Intersect correct records with predefined grid ---------------------------

# Load grid shapefile (Equal Area, width = 100km)
grid <- read_sf('shapefiles/GRID_EA100.shp') 
st_crs(grid)

# Transform occurrences into spatial points
data_correct$x <- data_correct$decimalLongitude
data_correct$y <- data_correct$decimalLatitude

data_correct <- st_as_sf(x = data_correct,                         
                        coords = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs")

# Reproject records to grid projection
data_correct <- st_transform(data_correct, crs = st_crs(grid))

# Point in polygon: join each record to a cell of the grid:
data_correct <- st_join(data_correct, grid)
data_correct$geometry <- NULL

# Save output
fwrite(data_correct, 'Data.csv', sep='\t') 
