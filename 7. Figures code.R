# The following script describe how to create figures and of the paper

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(tidyverse) # data manipulation and plot
library(patchwork)  # plot
library(sf)  # GIS 
library(rnaturalearth)  # GIS 
library(rnaturalearthdata)  # GIS 

# Environment and data ----------------------------------------------------
rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)

# Load 100x100km grid 
grid <- read_sf('shapefiles/GRID_EA100.shp')
# Load base map of the world
world <- ne_countries(scale = "medium", returnclass = "sf")  

# Load dataset obtained in script 3 of records and richness by cell
data <- read.csv("quantification.csv", sep=";") 

# Create shapefile to plot data in grid maps
data_shp <- merge(grid, data, by = "id")

### Map template
template <- function(data, var, title){
  ggplot(data) +
    geom_sf(data = grid, fill = 'lightblue', color ="transparent") +
    geom_sf(data = world, fill ='lightgrey', color = 'darkgrey') +
    geom_sf(aes(fill = var), color ="transparent") +
    coord_sf(crs = st_crs(8857))+
    scale_fill_viridis_c(limits = c(0, 1), option = "plasma", name = "", direction = -1) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(size = 16, face = "bold", hjust = 1.1)) +
    ylim(2800000, 8500000) +
    xlim(-14500000, 14500000) +
    labs(title = title)
} 
# Figure 1. Records distribution ####
fig1 <- template(data_shp, data_shp$richness3, "") +
  scale_fill_viridis_c(option = "plasma", name = "", direction = -1) +
  theme(legend.position = 'right')

ggsave('figure1.png', fig1, width = 30, height = 8, units = "cm")

# Figure 3 Preserved specimens filter ####
rPres <- template(data_shp, data_shp$prop_recPres, '')  # Number of records
sPres <- template(data_shp, data_shp$prop_richPres, '')  # Number of observed species

presSp <- rPres + sPres +
  plot_layout(guides = "collect", nrow = 2) & theme(legend.position = "right")
ggsave('figure3.png', presSp, width = 30, height = 15, units = "cm")

# Figure 4. Temporal distribution of records and observed species filter ####
r1600 <- template(data_shp, data_shp$p_rec1600, expression('' >= '1600'))
s1600 <- template(data_shp, data_shp$p_rich1600, '1600') + theme(plot.title = element_blank())
r1900 <- template(data_shp, data_shp$p_rec1900, expression('' >= '1900'))
s1900 <- template(data_shp, data_shp$p_rich1900, '1900') + theme(plot.title = element_blank())
r1970 <- template(data_shp, data_shp$p_rec1970, expression('' >='1970'))
s1970 <- template(data_shp, data_shp$p_rich1970, '1970') + theme(plot.title = element_blank())
r7000 <- template(data_shp, data_shp$p_rec7000, expression('  1970 - 2000'))
s7000 <- template(data_shp, data_shp$p_rich7000, '1970-2000') + theme(plot.title = element_blank())
r1990 <- template(data_shp, data_shp$p_rec1990, expression('' >= '1990'))
s1990 <- template(data_shp, data_shp$p_rich1990, '1990') + theme(plot.title = element_blank())


fig4 <- r1600+s1600+r1900+s1900+r1970+s1970+r7000+s7000+r1990+s1990+
      plot_layout(guides = "collect", ncol = 2, nrow = 5) & theme(legend.position = "bottom")
ggsave('figure4.png', fig4, width = 30, height = 30, units = "cm")

# Supp. Mat. 3 Duplicates records filter ####
dup1 <- template(data_shp, data_shp$p_Date, 'A') + 
        theme(plot.title = element_text(hjust = 0.1))

dup2 <- template(data_shp, data_shp$p_recByDate, 'B') + 
        theme(plot.title = element_text(hjust = 0.1)) 

dupPlot <- dup1 + dup2 +
  plot_layout(guides = "collect", nrow = 2) & theme(legend.position = "right")
ggsave('SuppMat3.png', dupPlot, width = 30, height = 30, units = "cm")

# Figure 5. Geographic location of well-sampled cells####
# Load dataframe obtained in script 5
data <- read.csv("countWellSurvey.csv", sep = ";")
data_shp <- merge(grid, data, by.x = "id", by.y = "Area")

fig5 <- template(data_shp, data_shp$total, "") +
        scale_fill_viridis_c(option = "plasma", name = "", direction = -1) +
        theme(legend.position = 'right')
ggsave('figure5.png', fig5, width = 30, height = 15, units = "cm")
# Figure 2. Histograms of cells' proportion of each class by region ####
# Filter fields of proportions by tax. source
data <- data %>% 
        select(c('id', 'region', 'p_richGBIF', 'p_richTPL', 'p_richTNRS', 'p_richTROP', 'p_richWFO')) %>%
        mutate(region = case_when(region == 'am' ~ 1,  # reclassify regions for loop
                                  region == 'eu' ~ 2,
                                  region == 'as' ~ 3))
data <- rename(data, 'GBIF'='p_richGBIF', 'TPL'='p_richTPL',
               'TNRS'='p_richTNRS', 'TROP'='p_richTROP', 'WFO'='p_richWFO')
data[is.na(data)] <- 0 # assign 0 to cells with no data (NA)
list <- colnames(data[ ,3:7])  # extract name from tax. sources

# Create empty 3x3 matrix to calculate values:
n <- data.frame(matrix(, ncol = 3, nrow = 3)) # 3 regions and 3 levels of change
colnames(n) <- c('am', 'eu', 'as') # column by region
class <- c('Decrease', 'Equal', 'Increase') # classes of change rate
rownames(n) <- class # row by class

# Loop to reclassify proportional changes in species names by cell to 3 levels:
# Increases -> '> 0'
# Decreases -> ' < 0'
# No change -> '= 0'
# And quantify their proportion by region in each source 
for (j in 1:5){
    a <- list[[j]] # select name of tax. source by loop
    
    for (i in 1:length(n)){
          bot <- data[region == i, ] # subset by region
          n[1,i] <- nrow(bot[get(a) < 0, ]) / nrow(bot) #'Decrease'
          n[2,i] <- nrow(bot[get(a) == 0, ]) / nrow(bot) # 'equal'
          n[3,i] <- nrow(bot[get(a) > 0, ]) / nrow(bot) #'Increase'
    }
    n$class <- rownames(n)
    assign(a, n) # rename matrix to tax. source
}

# List output: 5 dataframes
dfs <- list(GBIF, TPL, TNRS, TROP, WFO)

# Make a template of histograms' plot function:
plot <- function(source, region){
         ggplot(source) +
          geom_bar(aes(x = class, fill = class, weight = source[[region]]), color='black') +
          scale_fill_manual(values = list(Equal = "#f7f7f7", 
                                          Increase = "#b8e186", 
                                          Decrease = "#f1b6da")) +
          labs(x = '', y = '') +
          scale_x_discrete(limits = class) +
          theme(panel.grid.minor = element_blank(),
                legend.position = "none",
                axis.text.x = element_blank(), 
                axis.text.y = element_text(size = 10, face = "bold"),
                panel.background = element_rect(fill = "transparent"),  # bg of the panel
                plot.background = element_rect(fill = "transparent", color = NA))  # bg of the plot
}

# Loop each source and region to create each histogram and save them:
for (j in 1:5){
     source <- dfs[[j]]
     source_name <- list[j]
     myplot <- plot(source, 1)
     assign(paste0(source_name,'_AM'), myplot)
     myplot <- plot(source, 2)
     assign(paste0(source_name,'_EU'), myplot)
     myplot <- plot(source, 3)
     assign(paste0(source_name,'_AS'), myplot)
 }

newWFOhistAM <- WFO_AM +
  labs(caption ="AM") + 
  theme(axis.text.x = element_text(angle = 40),
        plot.caption = element_text(hjust = 0.5, size = rel(1.2)))
newWFOhistEU <- WFO_EU + 
  labs(caption ="EU / AF") + 
  theme(axis.text.x = element_text(angle = 40),
        plot.caption = element_text(hjust = 0.5, size = rel(1.2)))
newWFOhistAS <- WFO_AS +
  labs(caption ="AS") + 
  theme(axis.text.x = element_text(angle = 40),
        plot.caption = element_text(hjust = 0.5, size = rel(1.2)))

hists <- (GBIF_AM | GBIF_EU | GBIF_AS ) /
        (TPL_AM | TPL_EU | TPL_AS) /
        (TNRS_AM | TNRS_EU | TNRS_AS) /
        (TROP_AM | TROP_EU | TROP_AS )/
        (newWFOhistAM | newWFOhistEU | newWFOhistAS)
ggsave('figure2hists.png', hists, width = 12, height = 30, units = "cm")

# Maps of richness change by tax. source
data_shp <- merge(grid, data, by = "id")  # Join dataset info with grid shapefile by id of cell
# Template to create maps
template <- function(var, title){
  g <- data_shp %>%  mutate(val = case_when(data_shp[[var]] < 0 ~ 'dec',
                                            data_shp[[var]] == 0 ~ 'eq',
                                            data_shp[[var]] > 0 ~ 'inc'))    
  ggplot() +
    geom_sf(data = grid, fill = 'lightblue', color = "transparent") +
    geom_sf(data = world, fill ='lightgrey', color = 'darkgrey') +
    geom_sf(data = g, aes(fill = val), color = "transparent") +
    scale_fill_manual(values =  list(eq = "#f7f7f7", 
                                     inc = "#b8e186", 
                                     dec = "#f1b6da")) +
    coord_sf(crs = st_crs(8857))+
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          legend.position = "none", 
          plot.title = element_blank()) +
    ylim(2550000, 8500000) +
    xlim(-15000000, 15000000) +
    labs(title = title)
}

# Apply template to each tax. source
gbifMap <- template(3, 'GBIF')
tplMap <- template(4, 'TPL')
tnrsMap <- template(5, 'TNRS')
tropMap <- template(6, 'Tropicos')
wfoMap <- template(7, 'WFO')

# Create the composition of maps and save
maps <- gbifMap / tplMap / tnrsMap / tropMap / wfoMap
ggsave('figure2maps.png', maps, width = 30, height = 30, units = "cm")

# Supp. Mat. 2: Density plots of richness changes rates by cell and source ####
# Plot by region
template <- function(var, color){
            ggplot(data) +  
            geom_density(aes(x = var), color = color, lwd = 1, linetype = 1) +
            xlim(c(-0.5, 0.5)) +  # x limits zoomed to 0.5 range
            facet_grid(region ~ ., scales = "free") + 
            theme_minimal() + 
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text = element_text(size = 12, face = "bold"),
                  strip.text.y = element_blank())
}
a <- template(data$p_richGBIF, '#d7191c') # GBIF
b <- template(data$p_richTPL, '#fdae61')  # TPL
c <- template(data$p_richTNRS, '#fee090')  # TNRS
d <- template(data$p_richTROP, '#abd9e9')  # TROPICOS
e <- template(data$p_richWFO, '#2c7bb6')  # WFO

# Plot global data by source
template_glob <- function(var, color, title){
              ggplot(data) + 
              geom_density(aes(x = var), color = color, lwd = 1, linetype = 1) +
              xlim(c(-0.5, 0.5)) + 
              theme_minimal() + 
              theme(strip.text.y = element_blank(),
                    axis.text = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 12, face = "bold"),
                    axis.title.y = element_blank())+
              xlab(title) 
  }
aG <- template_glob(data$p_richGBIF, '#d7191c', "GBIF") # GBIF
bG <- template_glob(data$p_richTPL, '#fdae61', "TPL")  # TPL
cG <- template_glob(data$p_richTNRS, '#fee090', "TNRS")  # TNRS
dG <- template_glob(data$p_richTROP, '#abd9e9', "TROPICOS")  # TROPICOS
eG <- template_glob(data$p_richWFO, '#2c7bb6', "WFO")  # WFO

# Create the composition of plots
sm2 <- (a|b|c|d|e)/(aG|bG|cG|dG|eG) + plot_layout(heights = unit(c(15, 5), c('cm', 'cm')))
ggsave('suppMat2.png', sm2)

# Supp. Mat. 4: Distribution of completeness values based on Ratio, Records and Slope ####
# Load knowBR output estimators from script 4
data <- fread('knowBR/1600/Estimators.csv', header = TRUE, sep = ",")
# Create template to plot values
plot <- function(var, x, xtitle){
          ggplot(dt) +
          geom_point(aes(var, Completeness), pch = 19, size = 1.5) +
          geom_vline(xintercept = x, col = 'red', lwd = 1, lty = 2) +
          geom_hline(yintercept = 70, col = 'grey', lwd = 1, lty = 2) +
          theme_minimal() + 
          ylab('') +
          xlab(xtitle) +
          theme(strip.text.y = element_blank(),
                axis.text = element_text(size = 16, face = "bold"),
                axis.title = element_text(size = 16, face = "bold"))
 }
# Plot records value vs completeness (and threshold set in 100):
p1 <- plot(dt$Records, 100, 'Records') + ylab('Completeness')
# Plot ratio value vs completeness (and threshold set min 5):
p2 <- plot(dt$Ratio, 5, 'Ratio')
# Plot slope value vs completeness (and threshold set in 0.1):  
p3 <- plot(dt$Slope, 0.1, 'Slope')

sm4 <- p1|p2|p3
ggsave('SuppMat4.png', sm4, width = 30, height = 10, units = "cm")

# Supp. Mat. 5 and 6 Maps of 'well-surveyed cells for section 3 and 4' ####
template <- function(data, title){
  data <- data[Records >= 100, ][Completeness >= 70, ][Slope <= 0.1, ][Ratio >= 5, ][!is.na(Observed.richness), ][,'Area']
  data_shp <- merge(grid, data, by.x = "id", by.y = 'Area')
  n <- nrow(data_shp)        
  ggplot() +
  geom_sf(data = grid, fill = 'lightblue', color = "transparent") +
  geom_sf(data = world, fill = 'lightgrey', color = 'darkgrey') +
  geom_sf(data = data_shp, fill = "#B03060", color = 'transparent') +
  coord_sf(crs = st_crs(8857))+
  annotate("text", x = 13700000, y = 2850000, label = "n =") +
  annotate("text", x = 14500000, y = 2850000, label = n) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))+
  ylim(2800000, 8500000) +
  xlim(-14500000, 14500000) +
  labs(title = title) 
}
# Load estimators datasets obtained in Script 4 for the 8 filters:
kbr1600 <- fread('1600/Estimators.CSV', header = TRUE, sep = ",")
g1600 <- template(kbr1600, expression(bold('' >= '1600')))

kbr1900 <- fread('1900/Estimators.CSV', header = TRUE, sep = ",")
g1900 <- template(kbr1900, expression(bold('' >= '1900')))

kbr1970<- fread('1970/Estimators.CSV', header = TRUE, sep = ",")
g1970 <-  template(kbr1970, expression(bold('' >= '1970')))

kbr1990<- fread('1990/Estimators.CSV', header = TRUE, sep = ",")
g1990 <- template(kbr1990, expression(bold('' >= '1990')))

kbr70_00<- fread('70_00/Estimators.CSV', header = TRUE, sep = ",")
g70_00 <- template(kbr70_00, '1970 - 2000')

kbrpresSp <- fread('presSp/results.CSV', header = TRUE, sep = ",")
gpresSp <- template(kbrpresSp, 'Preserved Specimen')

kbrDuplDate <- fread('DuplDate/Estimators.CSV', header = TRUE, sep = ",")
gDuplDate <-template(kbrDuplDate, 'Duplicates with identical date of collection')

kbrDuplDateRec <- fread('DuplRecByDate/results.CSV', header = TRUE, sep = ",")
gDuplDateRec <- template(kbrDuplDateRec, "Duplicates with identical date of collection and recorders' names", n)

sm5 <- gDuplDateRec/gDuplDate/gpresSp/g1600/g1900/g1970/g1990/g70_00
ggsave('SuppMat5.png', sm5, width = 30, height = 60, units = "cm")

# Load estimators dataset obtained in Script 4 for the 3 filter combinations:
kbr1 <- fread('tnrs/Estimators.CSV', header = TRUE, sep = ",")
g1 <- template(kbr1, 'Dupl. by date + TNRS tax. standardization')

kbr2 <- fread('duplPres/Estimators.CSV', header = TRUE, sep = ",")
g2 <- template(kbr2, 'Dupl. by date + Preserved specimens')

kbr3 <- fread('dupl1970/Estimators.CSV', header = TRUE, sep = ",")
g3 <- template(kbr3, 'Dupl. by date + collected after 1970')

# Create map composition and save
sm6 <- g3/g2/g1
ggsave('SuppMat6.png', sm6, width = 30, height = 45, units = "cm")
