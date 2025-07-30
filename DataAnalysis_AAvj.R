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
library(DHARMa)
library(MuMIn)
library(gstat)
library(sp)
library(spdep)
library(openxlsx)

# data ----
data <- read_excel("data_bees.xlsx")
veg_data <- read_excel("data_vegetation.xlsx")
locations <- read_excel("locations.xlsx")

# research questions ----
# In this study, we aimed to find out how bee diversity changes in clear-cuts during succession (1–30 years),
# and which landscape (forest stand age and area) and vegetation descriptors best explain these changes. 

# Q1 -  How does forest stand age and area influence bee diversity and are there differences between sampling periods?
# Q2 -  Are there any flowering plant genera that are related to a higher bee diversity?
# Q3 -  Which vegetation descriptors help to better explain bee diversity?


# data preprocessing ----


# number of observation sessions per location
sessions=data %>% 
  summarise(period1=max(ifelse(Period==1,1,0)),
            period2=max(ifelse(Period==2,1,0)),
            periods=n(),.by = "Identificator")

# calculating bee diversity
bee_data <- data[6:59]
bee_data <- as.data.frame(bee_data)
bee_div <- diversity(bee_data,"shannon")
data$bee_div <- bee_div

# sampling period and sampling location identificators as factors
data$id=as.factor(data$Identificator)
data$fper=as.factor(data$Period)

# dataframe with mean values of bee diversity
data2 <- data %>%
  group_by(Identificator,id) %>%
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
# How does forest stand age and area influence bee diversity and are there differences between sampling periods?

hist(data$bee_div)
skewness=sum((data$bee_div-mean(data$bee_div))^3)/((length(data$bee_div)-1)*sd(data$bee_div)^3)
skewness
prop_zeros=length(data$bee_div[data$bee_div==0])/length(data$bee_div)
prop_zeros

mg0a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+(1|id),
                     data=data,family=gaussian(),REML=FALSE)
mg0b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg1a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+fper+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg1b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+fper+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg2a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+fper+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg2b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+fper+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg3a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+fper+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg3b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+fper+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg4a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg4b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg5a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg5b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)
mg6a=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+(1|id),
                      data=data,family=gaussian(),REML=FALSE)
mg6b=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+(1|id/fper),
                      data=data,family=gaussian(),REML=FALSE)

# mg6b - nonconvergence


mt0a=glmmTMB::glmmTMB(bee_div ~ 1+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt0b=glmmTMB::glmmTMB(bee_div ~ 1+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt1a=glmmTMB::glmmTMB(bee_div ~ 1+Area+Age+fper+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt1b=glmmTMB::glmmTMB(bee_div ~ 1+Area+Age+fper+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt2a=glmmTMB::glmmTMB(bee_div ~ 1+Area+s(Age)+fper+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt2b=glmmTMB::glmmTMB(bee_div ~ 1+Area+s(Age)+fper+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt3a=glmmTMB::glmmTMB(bee_div ~ 1+s(Area,Age)+fper+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt3b=glmmTMB::glmmTMB(bee_div ~ 1+s(Area,Age)+fper+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt4a=glmmTMB::glmmTMB(bee_div ~ 1+Area+Age+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt4b=glmmTMB::glmmTMB(bee_div ~ 1+Area+Age+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt5a=glmmTMB::glmmTMB(bee_div ~ 1+Area+s(Age)+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt5b=glmmTMB::glmmTMB(bee_div ~ 1+Area+s(Age)+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)
mt6a=glmmTMB::glmmTMB(bee_div ~ 1+s(Area,Age)+(1|id),
                      data=data,family=tw(link="log"),REML=FALSE)
mt6b=glmmTMB::glmmTMB(bee_div ~ 1+s(Area,Age)+(1|id/fper),
                      data=data,family=tw(link="log"),REML=FALSE)

# mt1b - nonconvergence
# mt2b - nonconvergence
# mt3b - nonconvergence
# mt6b - nonconvergence

comparison_table1=as.data.frame(MuMIn::AICc(mg0a,mg0b,mg1a,mg1b,mg2a,mg2b,mg3a,mg3b,mg4a,mg4b,mg5a,mg5b,mg6a,
            mt0a,mt0b,mt1a,mt2a,mt3a,mt4a,mt4b,mt5a,mt5b,mt6a))
comparison_table1$model=rownames(comparison_table1)
comparison_table1$delta=comparison_table1$AICc-min(comparison_table1$AICc)
comparison_table1=comparison_table1 %>% 
  dplyr::select(model,df,AICc,delta)
openxlsx::write.xlsx(comparison_table1,"Table_S2.xlsx")

# REML

mg1aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+fper+(1|id),
                      data=data,family=gaussian(),REML=TRUE)
mg2aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+fper+(1|id),
                      data=data,family=gaussian(),REML=TRUE)
mg3aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+fper+(1|id),
                      data=data,family=gaussian(),REML=TRUE)
mg4aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+Age+(1|id),
                      data=data,family=gaussian(),REML=TRUE)
mg5aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+Area+s(Age)+(1|id),
                      data=data,family=gaussian(),REML=TRUE)
mg6aR=glmmTMB::glmmTMB(log1p(bee_div) ~ 1+s(Area,Age)+(1|id),
                      data=data,family=gaussian(),REML=TRUE)


# residuals

DHARMa::testResiduals(DHARMa::simulateResiduals(mg1aR))
DHARMa::testResiduals(DHARMa::simulateResiduals(mg2aR))
DHARMa::testResiduals(DHARMa::simulateResiduals(mg3aR))
DHARMa::testResiduals(DHARMa::simulateResiduals(mg4aR))
DHARMa::testResiduals(DHARMa::simulateResiduals(mg5aR))
DHARMa::testResiduals(DHARMa::simulateResiduals(mg6aR))


res_mg1a=DHARMa::simulateResiduals(mg1aR)
plot(res_mg1a$scaledResiduals~data$Age)
plot(res_mg1a$scaledResiduals~data$Area)

res_mg2a=DHARMa::simulateResiduals(mg2aR)
plot(res_mg2a$scaledResiduals~data$Age)
plot(res_mg2a$scaledResiduals~data$Area)

res_mg3a=DHARMa::simulateResiduals(mg3aR)
plot(res_mg3a$scaledResiduals~data$Age)
plot(res_mg3a$scaledResiduals~data$Area)

res_mg4a=DHARMa::simulateResiduals(mg4aR)
plot(res_mg4a$scaledResiduals~data$Age)
plot(res_mg4a$scaledResiduals~data$Area)

res_mg5a=DHARMa::simulateResiduals(mg5aR)
plot(res_mg5a$scaledResiduals~data$Age)
plot(res_mg5a$scaledResiduals~data$Area)

res_mg6a=DHARMa::simulateResiduals(mg6aR)
plot(res_mg6a$scaledResiduals~data$Age)
plot(res_mg6a$scaledResiduals~data$Area)


# spatial autocorrelation

data_locs=left_join(data,locations,by="Identificator")
data_locs_sf=sf::st_as_sf(data_locs,coords=c("Longitude","Latitude"),crs=4326)
data_locs_proj=sf::st_transform(data_locs_sf,crs=3059)
koords=sf::st_coordinates(data_locs_proj)
data_locs_proj$x=koords[,1]
data_locs_proj$y=koords[,2]

res_mg1a=DHARMa::simulateResiduals(mg1aR)
data_locs_proj$res=res_mg1a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)


res_mg2a=DHARMa::simulateResiduals(mg2aR)
data_locs_proj$res=res_mg2a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)


res_mg3a=DHARMa::simulateResiduals(mg3aR)
data_locs_proj$res=res_mg3a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)


res_mg4a=DHARMa::simulateResiduals(mg4aR)
data_locs_proj$res=res_mg4a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)


res_mg5a=DHARMa::simulateResiduals(mg5aR)
data_locs_proj$res=res_mg5a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)

res_mg6a=DHARMa::simulateResiduals(mg6aR)
data_locs_proj$res=res_mg6a$scaledResiduals
ggplot(data_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  facet_wrap(~fper)+
  coord_fixed(ratio=1)
period1=data.frame(data_locs_proj[data_locs_proj$fper==1,])
coords1 <- cbind(period1$x, period1$y)
nb1 <- dnearneigh(coords1, 0, 20000)
lw1 <- nb2listw(nb1, style = "W")
moran.test(period1$res, lw1)
coordinates(period1)<-c("x","y")
Vario1 = gstat::variogram(scale(res) ~ 1, data=period1)
plot(Vario1)
period2=data.frame(data_locs_proj[data_locs_proj$fper==2,])
coords2 <- cbind(period2$x, period2$y)
nb2 <- dnearneigh(coords2, 0, 20000)
lw2 <- nb2listw(nb2, style = "W")
moran.test(period2$res, lw2)
coordinates(period2)<-c("x","y")
Vario2 = gstat::variogram(scale(res) ~ 1, data=period2)
plot(Vario2)

# summary
sjPlot::tab_model(mg1aR,mg2aR,mg3aR,mg4aR,mg5aR,mg6aR)


## prediction and visualisation ----
data$prognosis=expm1(predict(mg1aR,data))
prognosis_data=expand.grid(Age=seq(1,30,by=1),
                           Area=seq(0,5,0.25),
                           Period=c(1,2))
prognosis_data$prognosis=expm1(predict(mg6aR,prognosis_data, re.form = NA))

# Fig 2
labels <- c("Period 1", "Period 2")
names(labels) <- c("1", "2")
ggplot(prognosis_data, aes(Age, Area, fill = prognosis)) +
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
  )+
  ggview::canvas(dpi=600,width=174,height=70,units="mm")

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

# Using GLMs and GAMs - there is no need to consider pseudoreplication or seasonal variability, 
# as we are using the mean values of bee diversity, because of negligible effects in Q1
# we only use gaussian family with log1p(bee_div) because of the results in Q1

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

d1 <- glm(log1p(bee_div)~Age+Area+fl_div+fl_cov+Lysimachia+Campanula,data=data2)
dr4 <- MuMIn::dredge(d1,subset=!(Lysimachia&&fl_cov))


# exporting the dredge objects as excel files
openxlsx::write.xlsx(dr1,"dr1.xlsx")
openxlsx::write.xlsx(dr2,"dr2.xlsx")
openxlsx::write.xlsx(dr3,"dr3.xlsx")
openxlsx::write.xlsx(dr4,"dr4.xlsx")

## model selection ----

# Based on the AICc of the best models in dr1, dr2, dr3, and dr4 it is clear that dr2 and dr3
# results are identical and with lower AICc values even considering nonlinear effects
# We chose the second model (dr2) - nonlinear smoothed patch age as the best, as it 
# is a simpler model compared to also nonlinearly smoothed patch area

# Checking the response curves for models with delta < 2 in addition to AICc to choose 
# the best model

# 1st model from dr2
dr2_1 <- gam(log1p(bee_div)~s(Age) + fl_div, data=data2)
plot(dr2_1) # increasing uncertainty at the end with slight upwards curve

# 2nd model from dr2
dr2_2 <- gam(log1p(bee_div)~s(Age) + fl_div + Area, data=data2)
plot(dr2_2) # the end curves slightly downward

# 3rd model from dr2
dr2_3 <- gam(log1p(bee_div)~s(Age) + fl_div + Lysimachia, data=data2)
plot(dr2_3) # the end curves slightly upward

# We would expect the bee diversity to rise again at some point during the forest succession
# as windthrows or other natural disturbances occurs, but not as soon as 30 years after 
# clear-cutting, so we chose dr2_2 as the best model based on AICc values and the response curve
summary(dr2_2)

# residual evaluation

DHARMa::testResiduals(DHARMa::simulateResiduals(dr2_2))
res_dr2_2=DHARMa::simulateResiduals(dr2_2)
plot(res_dr2_2$scaledResiduals~data2$Age)
plot(res_dr2_2$scaledResiduals~data2$Area)
plot(res_dr2_2$scaledResiduals~data2$fl_div)
plot(res_dr2_2$scaledResiduals~data2$fl_cov)
plot(res_dr2_2$scaledResiduals~data2$Lysimachia)
plot(res_dr2_2$scaledResiduals~data2$Campanula)

data2_locs=left_join(data2,locations,by="Identificator")
data2_locs_sf=sf::st_as_sf(data2_locs,coords=c("Longitude","Latitude"),crs=4326)
data2_locs_proj=sf::st_transform(data2_locs_sf,crs=3059)
koords=sf::st_coordinates(data2_locs_proj)
data2_locs_proj$x=koords[,1]
data2_locs_proj$y=koords[,2]

data2_locs_proj$res=res_dr2_2$scaledResiduals
ggplot(data2_locs_proj,aes(y=y,x=x,col=res))+
  geom_point()+
  viridis::scale_color_viridis()+
  coord_fixed(ratio=1)
df_data2=data.frame(data2_locs_proj)
coords <- cbind(df_data2$x, df_data2$y)
nb <- dnearneigh(coords, 0, 20000)
lw <- nb2listw(nb, style = "W")
moran.test(df_data2$res, lw)
coordinates(df_data2)<-c("x","y")
Vario = gstat::variogram(scale(res) ~ 1, data=df_data2)
plot(Vario)

### model coefficients table ----
sjPlot::tab_model(dr2_2)
flextable::as_flextable(dr2_2)

## prognosis and visualisation ----

# Fig. 3
prognosis_data2=expand.grid(Age=seq(1,30,by=1),
                            Area=seq(0,5,0.25),
                            fl_div=c(0,0.5,1,1.5,2,2.5))
prognosis_data2$prognosis=expm1(predict(dr2_2,prognosis_data2))

dose.labs <- c("Plant diversity, H: 0", "Plant diversity, H: 0,5", "Plant diversity, H: 1",
               "Plant diversity, H: 1,5", "Plant diversity, H: 2", "Plant diversity, H: 2,5")
names(dose.labs) <- c("0", "0.5", "1", "1.5", "2", "2.5")

ggplot(prognosis_data2, aes(Age, Area, fill = prognosis)) +
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
  )+
  canvas(dpi=600,width=174,height=100,units="mm")

fig3 <- ggplot(prognosis_data2, aes(Age, Area, fill = prognosis)) +
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

#ggsave(plot = fig3, filename="Fig4.png",dpi=600,width=174,height=100,units="mm")

