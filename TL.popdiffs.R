

#### This code compares TL across the different YEWA breeding populations ####
#### TL vs GV ####


library(AICcmodavg)
library(lme4)
library(dplyr)
library(ggplot2)
library(tidybayes)
library(pbkrtest)


setwd("~/Desktop/CH1_TLvsGV/Telomere-Length-and-Genomic-Offset")

yewa=read.csv("YEWA_TL_data_2023_noHY.csv")

# Make separate "Year," "Month," and "Day" into one "Date" variable
yewa$Date <- paste(yewa$Year, yewa$Month, yewa$Day, sep=",") %>% ymd() %>% as.Date()

# Make "Date" variable into Julian Date ("Jdate") variable
yewa$Jday <- yday(yewa$Date)

# Standardize continuous variables
yewa$Plate <- as.factor(yewa$Plate)
yewa$Location <- as.factor(yewa$Location)
yewa$zElevation <- scale(yewa$Elevation)
yewa$zLat <- scale(yewa$Lat)
yewa$zJday <- scale(yewa$Jday)
yewa$zTarsus <- scale(yewa$Tarsus)

## telomere length is right-skewed, so log-transforming
TS_log = log(yewa$TS)
plotNormalHistogram(TS_log) # normal distribution
qqnorm(TS_log)
qqline(TS_log,
       col="red")

yewa$zTS <- scale(log(yewa$TS)) # will use log(TL) as response variable



##########################################################
### Test for population differences in telomere length ###
TL.pops <- lmer(zTS ~ Location + (1|Plate),data = yewa, REML=F)
TL.nopops <- lmer(zTS ~ (1|Plate),data = yewa, REML=F)

summary(TL.pops)
sjPlot::tab_model(TL.pops)

ggplot(yewa, aes(x=Location, y=zTS)) +
  geom_boxplot() +
  labs(y="Standardized Relative Log Telomere Length", x="Population") +
  theme(
    axis.title.x = element_text(size = 16),   # X-axis title text size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 10),    # X-axis tick labels text size
    axis.text.y = element_text(size = 14))





### Genomic Vulnerabilty
#### Candidate model set using all combinations of variables up to 4 variables ####
# explanatory variables: fut.offset, Age, Sex, Lat, Location, Elevation and Jday
# Response variable is zTS

fut.offset.mod1 <- lmer(zTS ~ GO + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod2 <- lmer(zTS ~ Sex + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod3 <- lmer(zTS ~ Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod4 <- lmer(zTS ~ zLat + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod5 <- lmer(zTS ~ zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod6 <- lmer(zTS ~ zElevation + (1|Plate) + (1|Location), data=yewa, REML = F)
fut.offset.mod7 <- lmer(zTS ~ zTarsus + (1|Plate) + (1|Location), data=yewa, REML = F)
fut.offset.mod8 <- lmer(zTS ~ GO * Sex + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod9 <- lmer(zTS ~ GO * Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod10 <- lmer(zTS ~ GO * zLat + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod11 <- lmer(zTS ~ GO * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod12 <- lmer(zTS ~ Sex * Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod13 <- lmer(zTS ~ Sex * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod14 <- lmer(zTS ~ Age * zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod15 <- lmer(zTS ~ Age * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod16 <- lmer(zTS ~ zLat * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod17 <- lmer(zTS ~ GO + Sex + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod18 <- lmer(zTS ~ GO + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod19 <- lmer(zTS ~ GO + zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod20 <- lmer(zTS ~ GO + zLat + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod21 <- lmer(zTS ~ GO + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod22 <- lmer(zTS ~ GO + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod23 <- lmer(zTS ~ Sex + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod24 <- lmer(zTS ~ Sex + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod25 <- lmer(zTS ~ Sex + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod26 <- lmer(zTS ~ Age + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod27 <- lmer(zTS ~ Age + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod28 <- lmer(zTS ~ Age + zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod29 <- lmer(zTS ~ zJday + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod30 <- lmer(zTS ~ zJday + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod31 <- lmer(zTS ~ zLat + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod32 <- lmer(zTS ~ zLat + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod33 <- lmer(zTS ~ GO * Sex * zJday +  (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod34 <- lmer(zTS ~ GO * Sex * zElevation +(1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod35 <- lmer(zTS ~ GO * Age * zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod36 <- lmer(zTS ~ GO * zLat * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod37 <- lmer(zTS ~ GO * zJday * zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod38 <- lmer(zTS ~ Age * zJday * zElevation +  (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod39 <- lmer(zTS ~ Sex * zElevation * zTarsus +  (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod40 <- lmer(zTS ~ Sex * zElevation * zJday +  (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod41 <- lmer(zTS ~ GO + Sex + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod42 <- lmer(zTS ~ GO + Sex + zLat + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod43 <- lmer(zTS ~ GO + Sex + zJday +  (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod44 <- lmer(zTS ~ GO + Sex + zElevation +(1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod45 <- lmer(zTS ~ GO + Sex + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod46 <- lmer(zTS ~ GO + Age + zLat + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod47 <- lmer(zTS ~ GO + Age + zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod48 <- lmer(zTS ~ GO + Age + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod49 <- lmer(zTS ~ GO + Age + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod50 <- lmer(zTS ~ GO + zLat + zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod51 <- lmer(zTS ~ GO + zLat + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod52 <- lmer(zTS ~ GO + zLat + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod53 <- lmer(zTS ~ GO + zJday + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod54 <- lmer(zTS ~ GO + zJday + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod55 <- lmer(zTS ~ GO + zElevation + zTarsus+ (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod56 <- lmer(zTS ~ GO * zJday + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod57 <- lmer(zTS ~ GO * zJday + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod58 <- lmer(zTS ~ GO * zJday + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod59 <- lmer(zTS ~ GO * zElevation + Sex + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod60 <- lmer(zTS ~ GO * zElevation + zJday + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod61 <- lmer(zTS ~ GO * zElevation + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod62 <- lmer(zTS ~ GO * zElevation + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod63 <- lmer(zTS ~ GO * Sex + zElevation + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod64 <- lmer(zTS ~ GO * Sex + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod65 <- lmer(zTS ~ GO * Sex + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod66 <- lmer(zTS ~ zJday* zElevation + GO + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod67 <- lmer(zTS ~ zJday* zElevation + Sex  + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod68 <- lmer(zTS ~ zJday* zElevation + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod69 <- lmer(zTS ~ zJday* zElevation + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod70 <- lmer(zTS ~ zJday* Sex + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod71 <- lmer(zTS ~ zJday* Sex + GO + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod72 <- lmer(zTS ~ zJday* Sex + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod73 <- lmer(zTS ~ zElevation * Sex + GO + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod74 <- lmer(zTS ~ zElevation * Sex + zTarsus + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod75 <- lmer(zTS ~ zElevation * Sex + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod76 <- lmer(zTS ~ zElevation * zTarsus + GO + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod77 <- lmer(zTS ~ zElevation * zTarsus + Age + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod78 <- lmer(zTS ~ zElevation * zTarsus + Sex + (1|Plate) + (1|Location), data = yewa, REML = F)
fut.offset.mod79 <- lmer(zTS ~ (1|Plate) + (1|Location), data = yewa, REML = F)


fut.offsetmodel_list <- list(fut.offset.mod1, fut.offset.mod2, fut.offset.mod3, fut.offset.mod4, fut.offset.mod5, fut.offset.mod6, fut.offset.mod7, fut.offset.mod8, fut.offset.mod9, fut.offset.mod10,
                             fut.offset.mod11, fut.offset.mod12, fut.offset.mod13, fut.offset.mod14, fut.offset.mod15, fut.offset.mod16, fut.offset.mod17, fut.offset.mod18, fut.offset.mod19, fut.offset.mod20,
                             fut.offset.mod21, fut.offset.mod22, fut.offset.mod23, fut.offset.mod24, fut.offset.mod25, fut.offset.mod26, fut.offset.mod27, fut.offset.mod28, fut.offset.mod29,fut.offset.mod30,
                             fut.offset.mod31, fut.offset.mod32, fut.offset.mod33, fut.offset.mod34, fut.offset.mod35, fut.offset.mod36, fut.offset.mod37, fut.offset.mod38, fut.offset.mod39,
                             fut.offset.mod40, fut.offset.mod41, fut.offset.mod42, fut.offset.mod43, fut.offset.mod44, fut.offset.mod45, fut.offset.mod46, fut.offset.mod47, fut.offset.mod48, fut.offset.mod49,
                             fut.offset.mod50, fut.offset.mod51, fut.offset.mod52, fut.offset.mod53, fut.offset.mod54, fut.offset.mod55, fut.offset.mod56, fut.offset.mod57, fut.offset.mod58, fut.offset.mod59,
                             fut.offset.mod60, fut.offset.mod61, fut.offset.mod62, fut.offset.mod63, fut.offset.mod64, fut.offset.mod65, fut.offset.mod66, fut.offset.mod67, fut.offset.mod68,
                             fut.offset.mod69, fut.offset.mod70, fut.offset.mod71, fut.offset.mod72, fut.offset.mod73, fut.offset.mod74, fut.offset.mod75, fut.offset.mod76, fut.offset.mod77, fut.offset.mod78,
                             fut.offset.mod79)

## Create AIC model selection table
fut.offset_AIC <- aictab(fut.offsetmodel_list, sort = TRUE, second.ord = TRUE)
fut.offset_AIC # Model 44 is AICw = 0.86

summary(fut.offset.mod61)

sjPlot::tab_model(fut.offset.mod61) # Conditional R2 = 0.366



###### Differnces in TL across populations considering various covariates

GO.mod1 <- lmer(zTS ~ Location +  (1|Plate), data=yewa)
GO.mod2 <- lmer(zTS ~ Age + (1|Plate) + (1|Location), data=yewa)
GO.mod3 <- lmer(zTS ~ Sex + (1|Plate) + (1|Location), data=yewa)
GO.mod4 <- lmer(zTS ~ GO + Age + (1|Plate) + (1|Location), data=yewa)
GO.mod5 <- lmer(zTS ~ GO + Sex + (1|Plate) + (1|Location), data=yewa)
GO.mod6 <- lmer(zTS ~ Age + Sex + (1|Plate) + (1|Location), data=yewa)
GO.mod7 <- lmer(zTS ~ GO + Age + Sex + (1|Plate) + (1|Location), data=yewa)
GO.mod8 <- lmer(zTS ~ GO * Age + (1|Plate) + (1|Location), data=yewa)
GO.mod9 <- lmer(zTS ~ GO * Sex + (1|Plate) + (1|Location), data=yewa)
GO.mod10 <- lmer(zTS ~ Age * Sex * (1|Plate) + (1|Location), data=yewa)
GO.mod11 <- lmer(zTS ~ (1|Plate) + (1|Location), data=yewa)

model_list <- list(GO.mod1,GO.mod2,GO.mod3,GO.mod4,GO.mod5,GO.mod6,GO.mod7,GO.mod8,GO.mod9,GO.mod10,GO.mod11)


## Create AIC model selection table
AIC <- aictab(model_list, sort = TRUE, second.ord = TRUE)
AIC # model 44 is AICw = 0.86


sjPlot::tab_model(GO.mod1)
summary(GO.mod1)


ggplot(yewa, aes(x=Elevation, y=zTS)) +
  geom_point()+
  geom_smooth(method = "lm", fullrange=F, alpha = 0.2)

ggplot(yewa, aes(x=GO, y=zTS)) +
  geom_point()+
  geom_smooth(method = "lm", fullrange=F, alpha = 0.2)

ggplot(yewa, aes(x=Tarsus, y=zTS)) +
  geom_point()+
  geom_smooth(method = "lm", fullrange=F, alpha = 0.2)



##### variation in qPCR values across populations

qPCR.mod1 <- lm(GAP.Cq ~ Plate, data=yewa)
sjPlot::tab_model(qPCR.mod1)
summary(qPCR.mod1)

qPCR.mod2 <- lm(GAP.Cq ~ Location, data=yewa)
sjPlot::tab_model(qPCR.mod2)
summary(qPCR.mod2)


qPCR.mod3 <- lm(TL.Eff ~ Location, data=yewa)
sjPlot::tab_model(qPCR.mod3)
summary(qPCR.mod3)


qPCR.mpd4 <- lm(GAP.eff ~ Location, data=yewa)
sjPlot::tab_model(qPCR.mpd4)
summary(qPCR.mpd4)



yewa.qPCR <- yewa %>%
  group_by(Plate) %>%
  summarize(
    Avg_TL.Cq = mean(TL.Cq, na.rm = TRUE),
    sd_TL.Cq = sd(TL.Cq, na.rm = TRUE),
    Avg_GAP.Cq = mean(GAP.Cq, na.rm = TRUE),
    sd_GAP.Cq = sd(GAP.Cq, na.rm = TRUE)
  )


mean(yewa.qPCR$Avg_TL.Cq)
sd(yewa.qPCR$Avg_TL.Cq)

mean(yewa.qPCR$Avg_GAP.Cq)
sd(yewa.qPCR$Avg_GAP.Cq)
