# R script for paper "Changes in bee (Hymenoptera: Anthophila) diversity during forest stand succession after final felling"
#install.packages("remotes")
#remotes::install_github("idmn/ggview")

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
data <- read_excel("data_clear-cut_bees.xlsx")

# research questions ----


# site selection ----
# soon to be come. I promise, this was carried out in R. 
## But in a terrifically messy file :(

# data preprocessing ----

# calculating bee diversity
bee_data <- data[6:59]
bee_data <- as.data.frame(bee_data)
bee_div <- diversity(bee_data,"shannon")
data$bee_div <- bee_div

# Question 1 ----


## modelling ----
##### GAMMs - influence of general forest characteristics on bee diversity

gamm1=gamm(log1p(bee_div) ~ s(Age) + s(Area) + factor(Period),
           random=list(Identificator=~1),data=data) 

gamm2=gamm(log1p(bee_div) ~ s(Age) + Area + factor(Period),
           random=list(Identificator=~1),data=data) 

gamm3=gamm(log1p(bee_div) ~ Age +s(Area) + factor(Period),
           random=list(Identificator=~1),data=data) 

gamm4=gamm(log1p(bee_div) ~ s(Age,Area) + factor(Period),
          random=list(Identificator=~1),data=data) 

## model selection ----
MuMIn::AICc(gamm1,gamm2,gamm3,gamm4)
# gamm4 has the lowest AICc
summary(gamm4$gam)
sjPlot::tab_model(gamm4)

## prognosis and visualisation ----
data$prognosis=expm1(predict(gamm4,data))
prognosis_data=expand.grid(Age=seq(1,30,by=1),
                           Area=seq(0,5,0.25),
                           Period=c(1,2))
prognosis_data$prognosis=expm1(predict(gamm4,prognosis_data))

labels <- c("period 1", "period 2")
names(labels) <- c("1", "2")
ggplot(prognosis_data, aes(Age, Area, fill = prognosis)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ Period, labeller = labeller(Reize = as_labeller(labels))) +
  theme_bw() +
  labs(x = "Forest stand age, years", y = "Forest stand area, ha", fill = "Bee diversity, H")

#ggview(dpi=600,width=200,height=159,units="mm")
#ggsave(filename="",dpi=600,width=200,height=159,units="mm")

# Question 2 ----

# clicked in PC-ORD

# Question 3 ----


# data frame preparation for GAMs - calculating mean values of bee diversity
data2 <- data %>%
  group_by(Identificator) %>%
  summarise(bee_div = ifelse(n() > 1, mean(bee_div, na.rm = TRUE),
                                  bee_div[!is.na(bee_div)]),
            Area = Area[1],
            Age = Age[1])

# vegetation data - flowering plants. Calculating diversity and total coverage
veg_data <- read_excel("data_vegetation.xlsx")
flower_coverage <- veg_data[2:34]
fl_div <- diversity(flower_coverage, "shannon")
veg_data$fl_div <- fl_div
veg_data$fl_cov <- rowSums(flower_coverage)

data2 <- merge(data2, veg_data[, c("Identificator", "fl_div", "fl_cov", "Lysimachia", "Campanula")],
               by = "Identificator", all.x = TRUE)
data2[is.na(data2)] <- 0



##### GAMs - influence of general forest characteristics 
#             and vegetation factors on bee diversity

# Checking for correlations >0.8
corel <- data2[,3:8]
cor(corel) 
cor.test(data2$Lysimachia, data2$fl_cov)
# Correlation between Lysimachia and fl_cov is 0.89 (p<0.001), we will exclude this combination in the models

options(na.action = na.fail)

a1 <- gam(log1p(bee_div)~Age+Area+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr1 <- MuMIn::dredge(a1,subset=!(Lysimachia&&fl_cov))
#

b1 <- gam(log1p(bee_div)~s(Age)+Area+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr2 <- MuMIn::dredge(b1,subset=!(Lysimachia&&fl_cov))

c1 <- gam(log1p(bee_div)~s(Age)+s(Area)+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr3 <- MuMIn::dredge(c1,subset=!(Lysimachia&&fl_cov))

d1 <- gam(log1p(bee_div)~Age+s(Area)+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr4 <- MuMIn::dredge(d1,subset=!(Lysimachia&&fl_cov))

gam <- gam(log1p(bee_div)~s(Age) + fl_div + Area, data=data2)
summary(gam)
sjPlot::tab_model(gam)

# prognosis and visualisation
gamviz <- ggpredict(gam,terms="Age [n=100]")
ggplot(gamviz,aes(x,predicted,ymin=conf.low,ymax=conf.high))+
  geom_ribbon(alpha=0.25, fill = "blue")+
  geom_line()+
  theme_classic() +
  labs(x="Forest stand age, years", y="Bee diversity, H")


prognosis_data2=expand.grid(Age=seq(1,30,by=1),
                            Area=seq(0,5,0.25),
                            fl_div=c(0,0.5,1,1.5,2,2.5))
prognosis_data2$prognosis=expm1(predict(gam,prognosis_data2))

dose.labs <- c("Plant diversity, H: 0", "Plant diversity, H: 0,5", "Plant diversity, H: 1",
               "Plant diversity, H: 1,5", "Plant diversity, H: 2", "Plant diversity, H: 2,5")
names(dose.labs) <- c("0", "0.5", "1", "1.5", "2", "2.5")

ggplot(prognosis_data2, aes(Age, Area, fill = prognosis)) +
  geom_raster() +
  facet_wrap(~ fl_div, labeller = labeller(fl_div = as_labeller(dose.labs))) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "Forest stand age, years", y = "Forest stand area, ha", 
       fill = "Bee diversity, H", main = "Augu daudzveidība, H")



