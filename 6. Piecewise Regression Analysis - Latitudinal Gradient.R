###======================================================### 
###             Piecewise regression analysis            ###
###              By Cristina Ronquillo 2022              ###
###======================================================### 

# The following script describe how to calculate piecewise regressions 
# and breakpoints of the latitudinal gradient of richness
# for well-sampled cells.
# It also includes the code for create figure 6.

# Load packages ---------------------------------------------------------
library(data.table)  # big dataset manipulation
library(segmented) # piecewise reg. analysis
library(tidyverse)  # plot
# Environment and data ----------------------------------------------------
rm(list = ls(all.names = TRUE))
wd <- ''  # Write working directory path
setwd(wd)

# Load latitudinal bands + region info by id cell
bands <- fread('LatBands_idCell.csv', header = TRUE, sep = ";")

# Load knowbr output by filter (obtained in script 4)
# (Change and run all script in each scenario)
i <- '' # Personalized objects generated
kbr <- fread('Estimators.csv', header = TRUE, sep = ";")

# Filter well-sampled cells with thresholds
dt <- kbr[Records >= 100, ][Completeness >= 70, ][Slope <= 0.1, ][Ratio >= 5, ][!is.na(Observed.richness), ]

# Join latitudinal band number and region to inventory completeness info
dt <- bands[dt, on = .(id = Area)]
# Number of well-sampled cells by region:
dt2 <- dt[, .(freq = .N), by = 'region']

# Create a function to make plots (data= dataset, p = max. number of cells per lat. band 
# q =  max. value of observed richness; s = p/q; text = x axis value to write leyend 
plot <- function(data, p, q, s, text){  
    ggplot(data) +
      aes(x = latitudeBands) +
      geom_linerange(aes(ymin = min, ymax = max)) +
      geom_point(aes(y = obs), size  = 1.5, colour = "#112446") +
      geom_point(aes(y = pred), size  = 1.5, colour = '#c51b8a') +
      geom_line(aes(y = bl_pred), colour = '#c51b8a') +
      geom_line(aes(y = bl_obv), colour = "#112446") + 
      geom_point(aes(y = freq*s), shape = 4, size  = 1.5, colour = "#ffa600") +
      geom_line(aes(y = bl_Freq*s), colour = "#ffa600") + 
      scale_y_continuous(name = "Richness", breaks = seq(0, q, by = 200), limits = c(0, q),
                         sec.axis = sec_axis(~./s, name = "Freq", breaks = seq(0, p, by = 10))) +
      geom_text(x = 9, y = text, label = paste0("bp freq = ", bp_Freq), size = 6) +
      geom_text(x = 6, y = text, label = paste0("bp obv = ", bp_obv), size = 6) +
      geom_text(x = 3, y = text, label = paste0("bp pred = ", bp_pred), size = 6) + 
      xlim(1,50)+
      theme_minimal() +    
      theme(plot.title = element_text(size = 22, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 12, face = "bold")) +
      coord_flip()
}

# 1. Calculate global median richness by latitudinal band ----------------------
dtGlob <- dt

t <- dtGlob[, list(pred = median(Richness)), by = c('latitudeBands')]  # predicted richness
t$min <- dtGlob[, list(min = min(Observed.richness)), by = c('latitudeBands')][,2]  # minimum observed richness
t$max <- dtGlob[, list(max = max(Observed.richness)), by = c('latitudeBands')][,2]  # maximum observed richness
t$obs <- dtGlob[, list(obs = median(Observed.richness)), by = c('latitudeBands')][,2]  # observed richness
t$freq <- dtGlob[, .(freq = .N), by = c('latitudeBands')][ ,2]  # number of well-sampled cells by lat. band

# Visualize to check data
plot(dtGlob$latitudeBands, dtGlob$Richness)  # SEE WHERE CURVE BREAKS APROX

# Calculate break points based on lm for:

### a) PREDICTED RICHNESS
fit <- lm(Richness ~ latitudeBands, data = dtGlob)
sg.fit <- segmented(fit, seg.Z =  ~ latitudeBands, psi = 30) # psi = value where curve breaks aprox
summary(sg.fit, short = TRUE)
bl_pred <- broken.line(sg.fit)  # regression line 
bl_pred <- bl_pred$fit  # prediction line values 
dtGlob <- cbind(dtGlob, bl_pred)  # add line values to data
bp_pred <- round(sg.fit$psi[1,2], 2)  # breakpoint value

# Save summary
sink(paste0("lmGlobalPred", i, ".txt"))
print(summary(sg.fit))
sink() 
# Visualization of regression
plot(dtGlob$latitudeBands, dtGlob$Richness)
plot(sg.fit, add = T, col = 'red')

### b) OBSERVED RICHNESS
fitOb <- lm(Observed.richness ~ latitudeBands, data = dtGlob)
sg.fitObv <- segmented(fitOb, seg.Z =  ~ latitudeBands, psi = 30)
summary(sg.fitObv, short = TRUE)
bl_obv <- broken.line(sg.fitObv)  # regression line 
bl_obv <- bl_obv$fit  # regression line values
dtGlob <- cbind(dtGlob, bl_obv)  # add line values to data
bp_obv <- round(sg.fitObv$psi[1,2], 2)  # breakpoint value

# Save summary
sink(paste0("lmGlobalObv", i, ".txt"))
print(summary(sg.fitObv))
sink() 
# Visualization of regression
plot(dtGlob$latitudeBands, dtGlob$Observed.richness)
plot(sg.fitObv, add = T, col = 'red')

### c) FREQUENCE of well-surveyed cells
fitFreq <- lm(freq ~ latitudeBands, data = t)
sg.fitFreq <- segmented(fitFreq, seg.Z =  ~ latitudeBands, psi = 30)
summary(sg.fitFreq, short=TRUE)
bl_Freq <- broken.line(sg.fitFreq)  # regression line 
bl_Freq <- bl_Freq$fit  # regression line values
bp_Freq <- round(sg.fitFreq$psi[1,2], 2)  # breakpoint value
t <- cbind(t, bl_Freq)  # add line values to data
# Save summary
sink(paste0("lmGlobalFreq", i, ".txt"))
print(summary(sg.fitFreq))
sink() 

dtGlob <- merge(dtGlob, t, by = 'latitudeBands', all.x = TRUE)

# Plot
p <- max(dtGlob$freq)  # max. number of cells per lat. band
q <- round(max(dtGlob$max))  # max. value of observed richness
s <- round(q/p)

(glob <- plot(dtGlob, p, q, s, 400) +
        labs(title = 'Global', x = 'Latitudinal band'))
ggsave('.png', glob, width = 10, height = 15, units = "cm")
# 2. Calculate median richness by latitudinal band by region -------------------
r <- dt[, list(pred = median(Richness)), by = c('latitudeBands', 'region')]
r$min <- dt[, list(min = min(Observed.richness)), by = c('latitudeBands', 'region')][ ,3]
r$max <- dt[, list(max = max(Observed.richness)), by = c('latitudeBands', 'region')][ ,3]
r$obs <- dt[, list(obs = median(Observed.richness)), by = c('latitudeBands', 'region')][ ,3]
r$freq <- dt[, .(freq = .N), by = c('latitudeBands', 'region')][,3]
r$region <- as.factor(r$region)

# Visualize to check break points and values
plot(r$pred, r$latitudeBands, col = r$region)

# 2a. AMERICA ------------------------------------------------------------------
rAm <- r[region == 'am', ]
dtAm <- dt[region == 'am', ]
# plot(dtAm$Richness, dtAm$latitudeBands)

# Calculate break points based on lm for:
###-- PREDICTED RICHNESS --###
fit <- lm(Richness ~ latitudeBands , data = dtAm)
sg.fit <- segmented(fit, seg.Z =  ~ latitudeBands, psi = 30) 
summary(sg.fit, short = TRUE)

bl_pred <- broken.line(sg.fit) # line 
bl_pred <- bl_pred$fit # 
bp_pred <- round(sg.fit$psi[1,2], 2) # breakpoint value 
dtAm <- cbind(dtAm, bl_pred) # extract line to plot in ggplot

# save summary
sink(paste0("lmAM_Pred",i,".txt"))
print(summary(sg.fit))
sink() 

# Visualization
plot(dtAm$latitudeBands, dtAm$Richness)
plot(sg.fit, add = T, col = 'red')

###-- OBSERVED RICHNESS --###
fitOb <- lm(Observed.richness ~ latitudeBands, data = dtAm)
sg.fitObv <- segmented(fitOb, seg.Z =  ~ latitudeBands, psi = 30)
summary(sg.fitObv, short = TRUE)
bl_obv <- broken.line(sg.fitObv) # line observed
bl_obv <- bl_obv$fit #line values observed
bp_obv <- round(sg.fitObv$psi[1,2], 2)# breakpoint value for observed
dtAm <- cbind(dtAm, bl_obv) #extract line to plot in ggplot

# save summary
sink(paste0("lmAM_Obv",i,".txt"))
print(summary(sg.fitObv))
sink() 
# visualization
plot(dtAm$latitudeBands, dtAm$Observed.richness) 
plot.segmented(sg.fitObv, add = TRUE, col = 'red')

###-- FREQUENCE of well-surveyed cells --###
fitFreq <- lm(freq ~ latitudeBands, data = rAm)
sg.fitFreq <- segmented(fitFreq, seg.Z =  ~ latitudeBands, psi = 30)
bl_Freq <- broken.line(sg.fitFreq)
bl_Freq <- bl_Freq$fit
bp_Freq <- round(sg.fitFreq$psi[1,2], 2) 
rAm <- cbind(rAm, bl_Freq)
# Save summary
sink(paste0("lmAmFreq_", i, ".txt"))
print(summary(sg.fitFreq))
sink()

dtAm <- merge(dtAm, rAm, by = 'latitudeBands', all.x = TRUE)

p1 <- max(rAm$freq)  # max. number of cells per lat. band
q1 <- round(max(rAm$max))  # max. value of observed richness
s1 <- round(q1/p1)

# Plot
(am <- plot(dtAm, p1, q1, s1, 350) +
       labs(title = 'America', x = ''))
ggsave('.png', am, width = 10, height = 15, units = "cm")

# 2b. EUROPE ------------------------------------------------------------------
rEu <- r[region == 'eu',]
dtEu <- dt[region == 'eu',]
plot(dtEu$Richness, dtEu$latitudeBands)

# Calculate break points based on lm for:
###-- PREDICTED RICHNESS --###
fit <- lm(Richness ~ latitudeBands, data = dtEu)
sg.fit <- segmented(fit, seg.Z =  ~ latitudeBands, psi = 30)
summary(sg.fit, short = TRUE)
bl_pred <- broken.line(sg.fit)
bl_pred <- bl_pred$fit 
bp_pred <- round(sg.fit$psi[1,2], 2)
dtEu <- cbind(dtEu, bl_pred) 

# save summary
sink(paste0("lmEU_Pred", i, ".txt"))
print(summary(sg.fit))
sink() 

# Visualization
plot(dtEu$latitudeBands, dtEu$Richness)
plot(sg.fit, add = T, col = 'red')

###-- OBSERVED RICHNESS --###
fitOb <- lm(Observed.richness ~ latitudeBands, data = dtEu)
sg.fitObv <- segmented(fitOb, seg.Z =  ~ latitudeBands, psi = 30)
summary(sg.fitObv, short = TRUE)

bl_obv <- broken.line(sg.fitObv)
bl_obv <- bl_obv$fit
dtEu <- cbind(dtEu, bl_obv)
bp_obv <- round(sg.fitObv$psi[1,2], 2)

# save summary
sink(paste0("lmEU_Obv", i, ".txt"))
print(summary(sg.fitObv))
sink() 
# Visualization
plot(dtEu$latitudeBands, dtEu$Observed.richness)
plot.segmented(sg.fitObv, add = TRUE, col = 'red')

###-- FREQUENCE of well-surveyed cells --###
fitFreq <- lm(freq ~ latitudeBands, data = rEu)
sg.fitFreq <- segmented(fitFreq, seg.Z =  ~ latitudeBands, psi = 30)
bl_Freq <- broken.line(sg.fitFreq)
bl_Freq <- bl_Freq$fit
rEu <- cbind(rEu, bl_Freq)
bp_Freq <- round(sg.fitFreq$psi[1,2], 2) 
summary(sg.fitFreq, short = TRUE)
# Save summary
sink(paste0("lmEuFreq_", i, ".txt"))
print(summary(sg.fitFreq))
sink()

dtEu <- merge(dtEu, rEu, by = 'latitudeBands', all.x = TRUE)

p2 <- max(rEu$freq)
q2 <- round(max(rEu$max))
s2 <- round(q2/p2)

# Plot
(eu <- plot(dtEu, p2, q2, s2, 350) +
       labs(title = 'Europe', x = ''))
ggsave('.png', eu, width = 10, height = 15, units = "cm")
# 2c. ASIA --------------------------------------------------------------------
rAs <- r[region == 'as',]
dtAs <- dt[region == 'as',]
plot(dtAs$Richness, dtAs$latitudeBands)

# Calculate break points based on lm for:
###-- PREDICTED RICHNESS --###
fit <- lm(Richness ~ latitudeBands , data = dtAs)
sg.fit <- segmented(fit, seg.Z =  ~ latitudeBands, psi = 20)
bl_pred <- broken.line(sg.fit)
bl_pred <- bl_pred$fit 
dtAs <- cbind(dtAs, bl_pred)
bp_pred <- round(sg.fit$psi[1,2], 2) 

summary(sg.fit, short=TRUE)
# save summary
sink(paste0("lmAS_Pred", i, ".txt"))
print(summary(sg.fit))
sink() 

# Visualization
plot(rAs$latitudeBands, rAs$pred)
plot(sg.fit, add = T, col = 'red')

###-- OBSERVED RICHNESS --###
fitOb <- lm(Observed.richness ~ latitudeBands, data = dtAs)
sg.fitObv <- segmented(fitOb, seg.Z =  ~ latitudeBands, psi = 20)
bl_obv <- broken.line(sg.fitObv)
bl_obv <- bl_obv$fit 
dtAs <- cbind(dtAs, bl_obv)
bp_obv <- round(sg.fitObv$psi[1,2], 2)

summary(sg.fitObv, short = TRUE)
# save summary
sink(paste0("lmAS_Obv", i, ".txt"))
print(summary(sg.fitObv))
sink() 
# Visualization
plot(rAs$latitudeBands, rAs$obs)
plot.segmented(sg.fitObv,add = TRUE, col = 'red')

###-- freq of well surveyed cells --###
fitFreq <- lm(freq ~ latitudeBands, data=rAs)
sg.fitFreq <- segmented(fitFreq, seg.Z =  ~ latitudeBands, psi=30)
bl_Freq <- broken.line(sg.fitFreq)
bl_Freq <- bl_Freq$fit
bp_Freq <- round(sg.fitFreq$psi[1,2], 2) 
#rAs <- cbind(rAs, bl_Freq)
# Save summary
sink(paste0("lmAsFreq_", i, ".txt"))
print(summary(sg.fitFreq))
sink()

dtAs <- merge(dtAs, rAs, by = 'latitudeBands', all.x = TRUE)

p3 <- max(rAs$freq)
q3 <- round(max(rAs$max))
s3 <- round(q3/p3)

(as <- plot(dtAs, p3, q3, s3, 250) +
      labs(title = 'Asia', x = ''))
ggsave('.png', as, width = 15, height = 15, units = "cm")

