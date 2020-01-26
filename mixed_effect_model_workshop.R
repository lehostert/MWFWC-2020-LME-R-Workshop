#Welcome to the intro to mixed effect modeling workshop
#Jim Lamer, Illinois River Biological Station, Illinois Natural History
   #Survey, Havana, IL  62644
     #Contact: email: lamer@illinois.edu, 309-543-6000

#Smallmouth buffalo catch per unit effort across LTRM study reaches
#bmbfcpue.csv
ACMaster<-read.csv(file.choose(),na.strings=".") # this brings your csv file in, the default of read.csv is to treat your header row as it is in the csv as the column headings/variable names 
#you can also enter directly in with ACMaster<- read.csv("ACMaster.csv")  The default for read.csv treats NA as missing values, if you have missing values represented by something else, you'll need to use na.strings="." if you missing values are represented by "."
str(ACMaster) # you should always check the struc
ACMaster$fyear<-factor(ACMaster$year)
ACMaster$ffstation<-factor(ACMaster$fstation)

#summary(ACMaster)

library(lattice)
library(lme4)
library(ggplot2)
library(dpyr)
library(tidyr)
ACMaster2 <- na.omit(ACMaster) #remove rows with missing values

#Let's start with linear model evaluation for temperature, current, gear
basic.lm <- lm(cpue ~ temp + current + factor(gear), data = ACMaster2)

summary(basic.lm)

drop1(basic.lm,test="F")#F = F test, no from normal distribution making p values relevant

step(basic.lm) #AIC comparison of model
plot(basic.lm, which = 1:4) #plot residuals

#Looks like temp is not significant in lm model
#Are we done, are we satisfied?
#What about variation between reaches, years, stratum?  Does it matter?
#Let's visualize the data and see
#first let's look at some of the factor variables
boxplot(cpue ~ ffstation*gear, data = ACMaster2)
boxplot(cpue ~ gear, data = ACMaster2)
boxplot(cpue ~ fyear, data = ACMaster2)
# certainly looks like something is going on here

#checking residuals to observe differences
E0 <- resid(basic.lm)
par(mar = c(5,5,2,2))
boxplot(E0 ~ ffstation, 
        data = ACMaster2,
        xlab = "Field station",
        ylab = "Residuals",
        cex.lab = 1.5)
abline(h = 0, lty = 2)

#Let's check out some of the continuous explanatory variables
ggplot(ACMaster2, aes(x = temp, y = cpue, colour = ffstation)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title="Temp_v_CPUE", fill = "Field Station")

ggplot(ACMaster2, aes(x = current, y = cpue, colour = ffstation)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title="Current_v_CPUE", fill = "Field Station")

ggplot(ACMaster2, aes(x = temp, y = cpue, colour = fyear)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title="Temp_v_CPUE", fill ="Year")

ggplot(ACMaster2, aes(x = current, y = cpue, colour = fyear)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title="Current_v_CPUE", fill = "Year")

plot.design(cpue~ffstation+fyear+stratum, data=ACMaster) #design plot



#run all analyses separately?  Not practical and compounds Type I error
ggplot(aes(temp, cpue), data = ACMaster2) + geom_point() +
  facet_wrap(~ fyear) +
  xlab("Temperature") + ylab("cpue")

ggplot(aes(temp, cpue), data = ACMaster2) + geom_point() +
  facet_wrap(~ ffstation) +
  xlab("Temperature") + ylab("cpue")

ggplot(aes(current, cpue), data = ACMaster2) + geom_point() +
  facet_wrap(~ fyear) +
  xlab("Current") + ylab("cpue")

ggplot(aes(current, cpue), data = ACMaster2) + geom_point() +
  facet_wrap(~ ffstation) +
  xlab("Current") + ylab("cpue")

#Could modify model to add covariates? #Let's try adding year to model
basic.lm2 <- lm(cpue ~ temp + current + factor(gear)+fyear, data = ACMaster2)
summary(basic.lm2)
#Besides adding 52 parameters to model, does that really help us 
#answer the question anyway?

#if we add these explanatory factor variables, what about sample size?
ftable(ffstation + stratum ~ fyear, data=ACMaster2)
#often need 10X sample size as parameters in model

#Next option - mixed effects model

#Let's run mixed effect model controlling for variation in cpue due to 
#differences in year (1|fyear) and differences in reach (1|ffstation)
buffmm <- lmer(cpue ~ temp + current + factor(gear)+
                 (1|fyear)+(1|ffstation), data = ACMaster2)
summary(buffmm)
#uh oh, no p values, what do we do?  
#We have a couple of options
#1.run model comparison, excluding variables and comparing models
buffmmA <- lmer(cpue ~ current + factor(gear)+
                 (1|fyear)+(1|ffstation), data = ACMaster2) #removes temp
buffmmB <- update(buffmmA, .~. - current)
buffmmC <- update(buffmmA, .~. - factor(gear))
anova(buffmm,buffmmA)
anova(buffmm,buffmmB)
anova(buffmm,buffmmC)

#or can use lmerTest (will automatically incorporate into summary after
#running package)
#library(lmerTest) #estimates df to give p values
#library(MuMIn) #This will return R2 values

#Nested variables
#Stratum may affect catch rates too? 
#make new nested variable since stratum and reach are likely correlated
#We can do it a couple ways, we can create new factor which is nested factor
ACMaster2 <- within(ACMaster2, stratatype <- factor(ffstation:stratum))
#This would end up looking like this (1|stratatype) or the same is 
#(1|ffstation/stratum)


buffmm1<-lmer(cpue~temp + current + factor(gear)+
                (1|fyear)+(1|ffstation)+(1|stratatype), data = ACMaster2)

summary(buffmm1)
(fixef(buffmm1))
(se.fixef(buffmm1))#returns fixed effect standard error
confint(buffmm1)

#which is a better model?
#we can compare random effects on model just like fixed effects 
#however, can only alter either random or fixed effects at same time
anova(buffmm,buffmm1)



#Visualize differences
#temp
ggplot(ACMaster2, aes(x = temp, y = cpue, colour = stratum)) +
  facet_wrap(~ffstation, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(ACMaster2, pred = predict(buffmm1)), aes(y = pred)) +
  theme(legend.position = "none")

#current
ggplot(ACMaster2, aes(x = current, y = cpue, colour = stratum)) +
  facet_wrap(~ffstation, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(ACMaster2, pred = predict(buffmm1)), aes(y = pred)) +
  theme(legend.position = "none")

#gives a way to present data
library(stargazer)
stargazer(buffmm1, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

plot(buffmm1, which = 1:4) # check model residuals

#Likely more appropriate to run as zero inflated model
#sum(ACMaster2$cpue == 0)/nrow(ACMaster2)
#Based on residuals, this should be run as glm, which glmer can be used

#IF trouble converging, you can standardize your covariates (mean 0, sd 1)
#using scale
ACMaster2$temp2 <- scale(ACMaster2$temp)
ACMaster2$current2 <- scale(ACMaster2$current)

buffmm2<-glmer(cpue~temp2 + current2 + factor(gear)+
                (1|fyear)+(1|ffstation)+(1|stratatype), 
               data = ACMaster2, family = poisson)
#overdispersion observed may lend itself to negative binomial distribution
#Mass is used for negative binomial with lme4
library(MASS)

buffmm2nb<-glmer.nb(cpue~temp + current + factor(gear)+
                 (1|fyear)+(1|ffstation), 
               data = ACMaster2, verbose = TRUE)
summary(buffmm2nb)

#check for overdispersion
E1<-resid(buffmm2, type = "pearson")
N<-nrow(ACMaster2)
p<-length(fixef(buffmm2)) + 1
Overdispersion<-sum(E1^2)/(N-p)
Overdispersion

#Can also model GAMM using nlme (General additive mixed modeling)

#Another example that adds dependency to model is repeated measures!
#I made up fish movement data (let's just say bowfin because they are cool)
#fishmovement.csv
move<-read.csv(file.choose(),na.strings=".")
str(move)
move$fID<-factor(move$ID)
move$Day<-factor(move$Day)
#Since the same fish is measured 25 times, we need to account for 
#that dependenceny
#What about day?
#each fish may respond differtly to each day?
#Let's look at lm first
model1<-lm(distance_km~Watertemp, data=move)
summary(model1)
plot(model1,which=1) #check distribution of residuals

#let's visualize and see what is going on.
ggplot(move, aes(x = Watertemp, y = distance_km)) +
  geom_point(size = 2) + geom_smooth(method = "lm")
  
#now let's add individual fish 
ggplot(move, aes(x = Watertemp, y = distance_km, colour = fID)) +
  geom_point(size = 2) + geom_smooth(method = "loess") +
  theme_classic() +
  theme(legend.position = "none")
#doesn't look like a linear relationship, but let's test anyway
model2<-lmer(distance_km~Watertemp + (1|fID), data=move)
summary(model2)
#let's verify with model comparison and remove temp
model2A<-lmer(distance_km~1 + (1|fID), data=move)
anova(model2,model2A)
model3<-lmer(distance_km~Watertemp + (1|fID) + (1|Day), data=move)
#does adding day as random effect help?
anova(model2,model3)
#model taking into account only individual fish ID the best

#The patterns in fish movement related to temperature is not a linear
#relationship, but should be explored further with GAMM (there 
#appears to be a pattern, but not linear)


#Larval fish example
#larvalfish.csv
larval<-read.csv(file.choose(),na.strings=".")
str(larval)
larval$Pool<-as.factor(larval$Pool)
larval$Hatch_Month<-factor(larval$Hatch_Month)
larval$SiteType<-factor(larval$SiteType)

#model for each family - we're interested in centrarchids and cyprinids
centrarchidae<-filter(larval, Family=="Centrarchidae")
cyprinidae<-filter(larval, Family=="Cyprinidae")
#How does temperature, site type, and river stage predict 
#centrarchid and cyprinid larval growth

#please work in groups to build model, and compare different random models


#Last group excercise:  Ebird data (from morning group)
#ebirdsurvey-jtl.csv
bird<-read.csv(file.choose(),na.strings=".")
str(bird)


#I pulled from several sources for the material in this workshop, 
#including UMESC fish graphical browser, UMRR, LTRM, Gabriela Hajduk, 
#Codingclub workshop, and Zuur et al. 2015 (Beginner's guid to 
#GLM and GLMM with R)

#For thos interested in ecological statistics and fitting 
#more advanced GLMMs and GAMMs, the highland statistics book collection
#is extremely useful and intuitive
