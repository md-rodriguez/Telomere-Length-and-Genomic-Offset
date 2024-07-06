##### Analyzing relationship between climate data (past, present, and future) and genomic offset values (past and future)


library(dplyr)
library(ggplot2)
library(tidyverse)
library(readr)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(gridExtra)
library(gtable)
library(grid)
library(gridtext)


setwd("~/Desktop/CH1_TLvsGV")

#Read in climate datasets for current, future, and past 

yewa <- read_csv("Data/Analysis_redo_2023/YEWA_TL_GV_ElevAb_2021.csv")


##### Telomere length and precip ######

yewa$LogTS <- log(yewa$TS) # will use log(TL) as response variable
yewa$zTS <- scale(yewa$TS)

yewa$pr_hist_ave2 <- yewa$pr_hist_ave^2
yewa$pr_hist_slope2 <- yewa$pr_hist_slope^2

#average TL across locations

yewa_sub <- yewa[, -1:-4]

TL_bySite <- yewa_sub %>%
  group_by(loc_num) %>%
  summarise_each(funs(mean))
  



TL_bySite <- yewa %>% group_by(loc_num) %>% summarise_each(funs(mean))

TL_bySite <- TL_bySite[,-c(2:3, 7:11)]




###### Past vs future climate comparisons ########

# Correlation between past and future Bio13 climate anomalies
Bio13_time_plot <- ggplot(TL_bySite, aes(x=bio13_Hdiffs, y=bio13_Fdiff)) +
  geom_point(size=3, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  labs(x="Historical Climate Change", y="Future Climate Change") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17)

Bio13_time_mod <- lm(bio13_Hdiffs ~ bio13_Fdiff, data = TL_bySite)
sjPlot::tab_model(Bio13_time_mod) # pvalue = <0.001  R2=0.76


# Correlation between past and future Bio15 climate anomalies
Bio15_time_plot <- ggplot(TL_bySite, aes(x=bio15_Hdiffs, y=bio15_Fdiff)) +
  geom_point(size=1.5, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  labs(x="Historical Climate Change (Bio15)", y="Future Climate Change (Bio15)") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17)

Bio15_time_mod <- lm(bio15_Hdiffs ~ bio15_Fdiff, data = TL_bySite)
sjPlot::tab_model(Bio15_time_mod) # pvalue = <0.001  R2=0.76


# Correlation between past and future Bio18 climate anomalies
Bio18_time_plot <- ggplot(TL_bySite, aes(x=bio18_Hdiffs, y=bio18_Fdiff)) +
  geom_point(size=1.5, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  labs(x="Historical Climate Anomaly (Bio18)", y="Future Climate Anomaly (Bio18)") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17)

Bio18_time_mod <- lm(bio18_Hdiffs ~ bio18_Fdiff, data = TL_bySite)
sjPlot::tab_model(Bio18_time_mod) # pvalue = <0.001  R2=0.76




######### Telomere/abundance trends vs bioclim trends ########

#TL vs bio13 ---- MS version
TL_prTrend1 <- ggplot(yewa, aes(x=bio13_Hdiffs, y=zTS)) +
  geom_point(size=1.5, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  labs(x="Change in Precipitation of the Wettest Month (Bio13)", y="Standardized Log Telomere Length") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17) +
  ylim(-2.5, 5)

#TL vs bio15 ---- MS version
TL_prTrend2 <- ggplot(yewa, aes(x=bio15_Hdiffs, y=zTS)) +
  geom_point(size=1.5, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  labs(x="Change in Precipitation Seasonality (Bio15)", y="Standardized Log Telomere Length") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17) +
  ylim(-2.5, 5)

TL_prTrend_mod <- lmer(zTS ~ bio13_Hdiffs + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
sjPlot::tab_model(TL_prTrend_mod) # bio13: pvalue = <0.023 bio15: pvalue = 0.009 R2 = 0.18
summary(TL_prTrend_mod)

# abundance trends vs past precip ave ---- MS version

ab_prtrend <- ggplot(TL_bySite, aes(y=abund, x=bio15_Hdiffs)) +
  geom_point(size=1.5, alpha=0.8, color="turquoise3") +
  geom_smooth(method = "lm",se=T, fullrange=T, color="turquoise3", alpha=0.3) + 
  theme_ipsum(base_size = 15)  +
  labs(y="Abundance Trends", x="Change in Preciptation Seasonality (Bio15)") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 17) +
  ylim(-3, 5)

ab_prAve_mod <- lm(abund ~ bio15_Hdiffs * Lat, data = TL_bySite)
summary(ab_prAve_mod)
sjPlot::tab_model(ab_prAve_mod) # bio15:pvalue = <0.022 lat:pvalue=0.015 int: pvalue=0.044 R2=38


### TL and precipitations
#### Candidate model set using all combinations of variables ####
# explanatory variables: precip, Lat, Loc_Num, Elevation, tmax, tmin
# Response variable is LogTS

TL.mod1 <- lmer(zTS ~ (1|loc_num), data = yewa, REML = F)
TL.mod2 <- lmer(zTS ~ bio18_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod3 <- lmer(zTS ~ Lat + (1|loc_num), data = yewa, REML = F)
TL.mod4 <- lmer(zTS ~ Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod5 <- lmer(zTS ~ bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod6 <- lmer(zTS ~ bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod7 <- lmer(zTS ~ bio18_Hdiffs + Lat + (1|loc_num), data = yewa, REML = F)
TL.mod8 <- lmer(zTS ~ bio18_Hdiffs + Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod9 <- lmer(zTS ~ bio18_Hdiffs + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod10 <- lmer(zTS ~ bio18_Hdiffs + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod11 <- lmer(zTS ~ Lat + Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod12 <- lmer(zTS ~ Lat + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod13 <- lmer(zTS ~ Lat + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod14 <- lmer(zTS ~ Elevation_1 + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod15 <- lmer(zTS ~ Elevation_1 + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod16 <- lmer(zTS ~ bio15_Hdiffs + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod17 <- lmer(zTS ~ bio18_Hdiffs + Lat + Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod18 <- lmer(zTS ~ bio18_Hdiffs + Lat + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod19 <- lmer(zTS ~ bio18_Hdiffs + Lat + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod20 <- lmer(zTS ~ Lat + Elevation_1 + bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod21 <- lmer(zTS ~ Lat + Elevation_1 + bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)

TL.mod22 <- lmer(zTS ~ bio18_Hdiffs * Lat + (1|loc_num), data = yewa, REML = F)
TL.mod23 <- lmer(zTS ~ bio18_Hdiffs * Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod24 <- lmer(zTS ~ bio18_Hdiffs * bio15_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod25 <- lmer(zTS ~ bio18_Hdiffs * bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod26 <- lmer(zTS ~ bio15_Hdiffs * Lat + (1|loc_num), data = yewa, REML = F)
TL.mod27 <- lmer(zTS ~ bio15_Hdiffs * Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod28 <- lmer(zTS ~ bio15_Hdiffs * bio13_Hdiffs + (1|loc_num), data = yewa, REML = F)
TL.mod29 <- lmer(zTS ~ bio13_Hdiffs * Lat + (1|loc_num), data = yewa, REML = F)
TL.mod30 <- lmer(zTS ~ bio13_Hdiffs * Elevation_1 + (1|loc_num), data = yewa, REML = F)
TL.mod31 <- lmer(zTS ~ bio15_Hdiffs + bio13_Hdiffs + bio18_Hdiffs + (1|loc_num), data = yewa, REML = F)




TLmodel_list <- list(TL.mod1, TL.mod2, TL.mod3, TL.mod4, TL.mod5, TL.mod6, TL.mod7, TL.mod8, TL.mod9, TL.mod10,
                             TL.mod11, TL.mod12, TL.mod13, TL.mod14, TL.mod15, TL.mod16, TL.mod17, TL.mod18, TL.mod19, 
                             TL.mod20, TL.mod21, TL.mod22, TL.mod23, TL.mod24, TL.mod25, TL.mod26, TL.mod27, TL.mod28, 
                     TL.mod29, TL.mod30, TL.mod31)

## Create AIC model selection table
TL_AIC <- aictab(TLmodel_list, sort = TRUE, second.ord = TRUE)
TL_AIC # Model 2 is AICw = 0.36

#### Results from top model: TL.mod2 ####

sjPlot::tab_model(TL.mod16,
                  show.re.var= TRUE, 
                  dv.labels= "Precip vs Telomere Length") # Conditional R2 = 0.466


### Abundance vs past precip
#### Candidate model set using all combinations of variables ####
# explanatory variables: precip, Lat, Loc_Num, Elevation, tmax, tmin
# Response variable is abundance

ab.mod1 <- lm(abund ~ 1, data = TL_bySite)
ab.mod2 <- lm(abund ~ bio18_Hdiffs, data = TL_bySite)
ab.mod3 <- lm(abund ~ Lat, data = TL_bySite)
ab.mod4 <- lm(abund ~ Elevation_1, data = TL_bySite)
ab.mod5 <- lm(abund ~ bio15_Hdiffs, data = TL_bySite)
ab.mod6 <- lm(abund ~ bio13_Hdiffs, data = TL_bySite)

ab.mod7 <- lm(abund ~ bio18_Hdiffs + Lat, data = TL_bySite)
ab.mod8 <- lm(abund ~ bio18_Hdiffs + Elevation_1, data = TL_bySite)
ab.mod9 <- lm(abund ~ bio18_Hdiffs + bio15_Hdiffs, data = TL_bySite)
ab.mod10 <- lm(abund ~ bio18_Hdiffs + bio13_Hdiffs, data = TL_bySite)

ab.mod11 <- lm(abund ~ Lat + Elevation_1, data = TL_bySite)
ab.mod12 <- lm(abund ~ Lat + bio15_Hdiffs, data = TL_bySite)
ab.mod13 <- lm(abund ~ Lat + bio13_Hdiffs, data = TL_bySite)

ab.mod14 <- lm(abund ~ Elevation_1 + bio15_Hdiffs, data = TL_bySite)
ab.mod15 <- lm(abund ~ Elevation_1 + bio13_Hdiffs, data = TL_bySite)

ab.mod16 <- lm(abund ~ bio15_Hdiffs + bio13_Hdiffs, data = TL_bySite)

ab.mod17 <- lm(abund ~ bio18_Hdiffs + Lat + Elevation_1, data = TL_bySite)
ab.mod18 <- lm(abund ~ bio18_Hdiffs + Lat + bio15_Hdiffs, data = TL_bySite)
ab.mod19 <- lm(abund ~ bio18_Hdiffs + Lat + bio13_Hdiffs, data = TL_bySite)

ab.mod20 <- lm(abund ~ Lat + Elevation_1 + bio15_Hdiffs, data = TL_bySite)
ab.mod21 <- lm(abund ~ Lat + Elevation_1 + bio13_Hdiffs, data = TL_bySite)

ab.mod22 <- lm(abund ~ bio18_Hdiffs * Lat, data = TL_bySite)
ab.mod23 <- lm(abund ~ bio18_Hdiffs * Elevation_1, data = TL_bySite)
ab.mod24 <- lm(abund ~ bio15_Hdiffs * Lat, data = TL_bySite)
ab.mod25 <- lm(abund ~ bio15_Hdiffs * Elevation_1, data = TL_bySite)
ab.mod26 <- lm(abund ~ bio15_Hdiffs * bio13_Hdiffs, data = TL_bySite)
ab.mod27 <- lm(abund ~ bio13_Hdiffs * Lat, data = TL_bySite)
ab.mod28 <- lm(abund ~ bio13_Hdiffs * Elevation_1, data = TL_bySite)
ab.mod29 <- lm(abund ~ bio18_Hdiffs * bio15_Hdiffs * bio13_Hdiffs, data = TL_bySite)

ab.mod30 <- lm(abund ~ bio15_Hdiffs + bio18_Hdiffs + Elevation_1, data = TL_bySite)
ab.mod31 <- lm(abund ~ bio13_Hdiffs + bio18_Hdiffs + Elevation_1, data = TL_bySite)
ab.mod32 <- lm(abund ~ bio13_Hdiffs + bio15_Hdiffs + Elevation_1, data = TL_bySite)

ab.mod33 <- lm(abund ~ bio15_Hdiffs * bio18_Hdiffs + Elevation_1, data = TL_bySite)
ab.mod34 <- lm(abund ~ bio13_Hdiffs * bio18_Hdiffs + Elevation_1, data = TL_bySite)
ab.mod35 <- lm(abund ~ bio13_Hdiffs * bio15_Hdiffs + Elevation_1, data = TL_bySite)

ab.mod36 <- lm(abund ~ bio15_Hdiffs * Elevation_1 + Lat, data = TL_bySite)
ab.mod37 <- lm(abund ~ bio13_Hdiffs * Elevation_1 + Lat, data = TL_bySite)
ab.mod38 <- lm(abund ~ bio13_Hdiffs * Elevation_1 + Lat, data = TL_bySite)



abmodel_list <- list(ab.mod1, ab.mod2, ab.mod3, ab.mod4, ab.mod5, ab.mod6, ab.mod7, ab.mod8, ab.mod9, ab.mod10,
                     ab.mod11, ab.mod12, ab.mod13, ab.mod14, ab.mod15, ab.mod16, ab.mod17, ab.mod18, ab.mod19, 
                     ab.mod20, ab.mod21, ab.mod22, ab.mod23, ab.mod24, ab.mod25, ab.mod26, ab.mod27, ab.mod28, 
                     ab.mod29, ab.mod30, ab.mod31, ab.mod32, ab.mod33, ab.mod34, ab.mod35, ab.mod36,
                     ab.mod37, ab.mod38)

## Create AIC model selection table
ab_AIC <- aictab(abmodel_list, sort = TRUE, second.ord = TRUE)
ab_AIC # Model 2 is AICw = 0.36

#### Results from top model: TL.mod2 ####

sjPlot::tab_model(ab.mod24,
                  show.re.var= TRUE, 
                  dv.labels= "Precip vs Telomere Length") # Conditional R2 = 0.466


##### AIC for TL vs offset and precip

pr_TL_groups$zLogTS <- log(pr_TL_groups$TS_Mean) # will use log(TL) as response variable


TL.OP.mod1 <- lmer(LogTS ~ fut.offset + (1|Loc_Num), data = pr_TL_groups, REML=T)
TL.OP.mod2 <- lmer(LogTS ~ future + (1|Loc_Num), data = pr_TL_groups, REML=T)
TL.OP.mod3 <- lmer(LogTS ~ future + fut.offset + (1|Loc_Num), data = pr_TL_groups, REML=T)

TL.OP_list <- list(TL.OP.mod1, TL.OP.mod2, TL.OP.mod3)

TL.OP.mod_AIC <- aictab(TL.OP_list, sort = TRUE, second.ord = TRUE)
TL.OP.mod_AIC # Model 44 is AICw = 0.71


