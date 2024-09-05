# R script for paper "Changes in bee (Hymenoptera: Anthophila) diversity during forest stand succession after final felling"

# libs ----
library(readxl)
library(tidyverse)
library(mgcv)
library(modelsummary)
library(flextable)
library(vegan)
library(dplyr)
library(mgcViz)
library(ggeffects)
library(ggview)

# data ----
data <- read_excel("data_bees.xlsx")
veg_data <- read_excel("data_vegetation.xlsx")


# research questions ----
# In this study, we aimed to find out how bee diversity changes in clear-cuts during succession (1–30 years),
# and which landscape (forest stand age and area) and vegetation descriptors best explain these changes. 

# Q1 -  How does forest stand age and area influence bee diversity?
# Q2 -  Are there any flowering plant genera that are related to a higher bee diversity?
# Q3 -  Which vegetation descriptors help to better explain bee diversity?


# data preprocessing ----

# calculating bee diversity
bee_data <- data[6:59]
bee_data <- as.data.frame(bee_data)
bee_div <- diversity(bee_data,"shannon")
data$bee_div <- bee_div

# dataframe preparation for GAMs - calculating mean values of bee diversity
data2 <- data %>%
  group_by(Identificator) %>%
  summarise(bee_div = ifelse(n() > 1, mean(bee_div, na.rm = TRUE),
                             bee_div[!is.na(bee_div)]),
            Area = Area[1],
            Age = Age[1])

# Calculating diversity and total coverage of flowers
flower_coverage <- veg_data[2:34]
fl_div <- diversity(flower_coverage, "shannon")
veg_data$fl_div <- fl_div
veg_data$fl_cov <- rowSums(flower_coverage)
data2 <- merge(data2, veg_data[, c("Identificator", "fl_div", "fl_cov", "Lysimachia", "Campanula")],
               by = "Identificator", all.x = TRUE)
data2[is.na(data2)] <- 0


# Question 1 ----
# How does forest stand age and area influence bee diversity?

## modelling ----
# Using GAMM - we expected that the response function of the stand may not be linear due to its 
# age-related characteristics

mod1=gamm(log1p(bee_div) ~ s(Age) + s(Area) + factor(Period),
           random=list(Identificator=~1),data=data) 

mod2=gamm(log1p(bee_div) ~ s(Age) + Area + factor(Period),
           random=list(Identificator=~1),data=data) 

mod3=gamm(log1p(bee_div) ~ s(Age,Area) + factor(Period),
          random=list(Identificator=~1),data=data) 

## model selection ----
MuMIn::AICc(mod1,mod2,mod3)
# mod3 with factor interaction has the lowest AICc, so we selected this one
summary(mod3$gam)
sjPlot::tab_model(mod3)

## prognosis and visualisation ----
data$prognosis=expm1(predict(mod3,data))
prognosis_data=expand.grid(Age=seq(1,30,by=1),
                           Area=seq(0,5,0.25),
                           Period=c(1,2))
prognosis_data$prognosis=expm1(predict(mod3,prognosis_data))

# Fig 2
labels <- c("Period 1", "Period 2")
names(labels) <- c("1", "2")
fig2 <- ggplot(prognosis_data, aes(Age, Area, fill = prognosis)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ Period, labeller = labeller(Period = as_labeller(labels))) +
  theme_bw() +
  labs(x = "Forest stand age, years", y = "Forest stand area, ha", fill = "Bee diversity, H") +
  theme(
    axis.text = element_text(size = 8),  
    axis.title = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9)
  )

ggview(plot = fig2,dpi=600,width=174,height=70,units="mm")
#ggsave(plot = fig2,filename="Fig2.png",dpi=600,width=174,height=70,units="mm")


# Question 2 ----
# Are there any flowering plant genera that are related to a higher bee diversity?

# This was carried out in PC-ORD

# The acquired stress level of NMS was 12.66. 
# The second axis was positively correlated with forest stand age (τ = 0.520, p < 0.001) and negatively correlated 
# with bee diversity (τ = - 0.624, p < 0.001). Additionally, the flower coverage of three plant genera were negatively 
# correlated with the second axis: Lysimachia (τ  = -0.389, p = 0.014), Campanula (τ = -0.226, p = 0.040) and 
# Galeopsis (τ = -0.389, p = 0.035), although the correlation was weak. Since we observed Galeopsis plant genera only in 
# two of the sample sites, we did not include it in further analysis.


# Question 3 ----
# Which vegetation descriptors help to better explain bee diversity?

# Using GAM - there is no need to consider pseudoreplication or seasonal variability, as we are using the mean values
# of bee diversity, because no significant difference between the sample periods was found (p=0.096)

## corellations ----
# Checking for correlations >0.59 to exclude them from further analysis, as is the usual procedure, since strongly 
# correlated variables reduce descriptive and out-of-sample predictive power
corel <- data2[,3:8]
cor(corel)
cor.test(data2$Lysimachia, data2$fl_cov)
# Correlation between Lysimachia and fl_cov is 0.89 (p<0.001), we will exclude this combination in the models

## modelling ----
# generating all possible factor combinations using the "dredge" function
options(na.action = na.fail)

a1 <- gam(log1p(bee_div)~s(Age,Area)+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr1 <- MuMIn::dredge(a1,subset=!(Lysimachia&&fl_cov))

b1 <- gam(log1p(bee_div)~s(Age)+Area+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr2 <- MuMIn::dredge(b1,subset=!(Lysimachia&&fl_cov))

c1 <- gam(log1p(bee_div)~s(Age)+s(Area)+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr3 <- MuMIn::dredge(c1,subset=!(Lysimachia&&fl_cov))

## model selection ----

# Based on the AICc of the best models in dr1, dr2, dr3 it is clear that dr2 and dr3
# results are identical, therefore we chose the second model (dr2) as the best, as it 
# is a simpler model

# Checking the response curves for models with delta < 2 in addition to AICc to choose 
# the best model

# 1st model from dr2
dr2_1 <- gam(log1p(bee_div)~s(Age) + fl_div, data=data2)
plot(dr2_1) # the end curves upward

# 2nd model from dr2
dr2_2 <- gam(log1p(bee_div)~s(Age) + fl_div + Area, data=data2)
plot(dr2_2) # the end curves slightly downward

# 3rd model from dr2
dr2_3 <- gam(log1p(bee_div)~s(Age) + fl_div + Lysimachia, data=data2)
plot(dr2_3) # the end curves slightly upward

# We would expect the bee diversity to rise again at some point during the forest succession
# as windthrow occurs, but not as soon as 30 years after clear-cutting, so we chose dr2_2
# as the best model based on AICc values and the response curve

## prognosis and visualisation ----

# Fig. 3
# This figure is an addition to answer Q1
gamviz <- ggpredict(dr2_2,terms="Age [n=100]")
fig3 <- ggplot(gamviz, aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(alpha = 0.25, fill = "blue") +
  geom_line() +
  theme_classic() +
  labs(x = "Forest stand age, years", y = "Bee diversity, H") +
  theme(
    axis.text = element_text(size = 8),  
    axis.title = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9)
  )

ggview(plot = fig3,dpi=600,width=84,height=60,units="mm")
#ggsave(plot = fig3, filename="Fig3.png",dpi=600,width=84,height=60,units="mm")

# Fig. 4
prognosis_data2=expand.grid(Age=seq(1,30,by=1),
                            Area=seq(0,5,0.25),
                            fl_div=c(0,0.5,1,1.5,2,2.5))
prognosis_data2$prognosis=expm1(predict(dr2_2,prognosis_data2))

dose.labs <- c("Plant diversity, H: 0", "Plant diversity, H: 0,5", "Plant diversity, H: 1",
               "Plant diversity, H: 1,5", "Plant diversity, H: 2", "Plant diversity, H: 2,5")
names(dose.labs) <- c("0", "0.5", "1", "1.5", "2", "2.5")

fig4 <- ggplot(prognosis_data2, aes(Age, Area, fill = prognosis)) +
  geom_raster() +
  facet_wrap(~ fl_div, labeller = labeller(fl_div = as_labeller(dose.labs))) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Forest stand age, years", y = "Forest stand area, ha", 
       fill = "Bee diversity, H") +
  theme(
    axis.text = element_text(size = 8),  
    axis.title = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9)
  )

ggview(plot = fig4,dpi=600,width=174,height=100,units="mm")
#ggsave(plot = fig4, filename="Fig4.png",dpi=600,width=174,height=100,units="mm")

