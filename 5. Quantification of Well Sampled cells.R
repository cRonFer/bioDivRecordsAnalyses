###======================================================### 
###         Quantification of well-surveyed cells        ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

### The following script describe the steps to quantify the number of
### well-sampled cells after the inventory completeness analyses of each
### filter considered

# Load packages ---------------------------------------------------------
library(tidyverse)  # data manipulation

# Environment and data ----------------------------------------------------

rm(list = ls(all.names = TRUE))
wd <- ''  # Write the working directory path
setwd(wd)

# Load each output of estimators from knowBR analyses:
d1600 <- read.csv('1600/Estimators.CSV', header = TRUE, sep = ",")
d1900 <- read.csv('1900/Estimators.CSV', header = TRUE, sep = ",")
d1970 <- read.csv('1970/Estimators.CSV', header = TRUE, sep = ",")
d1990 <- read.csv('1990/Estimators.CSV', header = TRUE, sep = ",")
d7000 <- read.csv('70_00/Estimators.CSV', header = TRUE, sep = ",")
date <- read.csv('DuplDate/Estimators.CSV', header = TRUE, sep = ",")
recByDate <- read.csv('DuplRecByDate/Estimators.CSV', header = TRUE, sep = ",")
pres <- read.csv('presSp/Estimators.CSV', header = TRUE, sep = ",")

names_of_dataframes <- ls.str(mode = "list")

# Filter estimators based on previously established thresholds ####
# (choose the dataset with the highest value of well-surveyed cells)
dt <- d1600 %>% 
        filter(Completeness >= 70 & Ratio >= 5 & Records >= 100 & Slope <= 0.1) %>%  
        select(Area) 

dt$WellS <- 1 # Add count field to well surveyed cells
setnames(dt, 'WellS', names_of_dataframes[1])

# List datasets
l <- list(d1900, d1970, d1990, d7000, pres, date, recByDate)

# Loop for filter cells of each dataset
for (i in 1:7){
  d <- l[[i]]
  d <- d %>% filter(Completeness >= 70 & Ratio >= 5 & Records >= 100 & Slope <= 0.1)
  d <- d %>%  select(Area)
  d$WellS <- 1
  setnames(d, 'WellS', names_of_dataframes[i + 1])
  dt <- merge(dt, d, all.x = TRUE, by = 'Area')
}
dt[is.na(dt)] <- 0 # Add 0 to NA cells

# Sum-up all filters and count well-surveyed freq
dt <- dt %>% mutate(total = d1600 + d1900 + d1970 + d1990 + 
                            d7000 + date + pres + recByDate)
# Save output
write.table(dt, '.csv', sep=";", row.names = FALSE)
