## -------------------------------------------------- ##
## EDSD 2024-2025: Demographic Forecasting
## Lecture 1
## Direct extrapolation by (generalized) linear models
## Date: 04/05/2025
## Instructor: Ugofilippo Basellini
## -------------------------------------------------- ##

## ---- EXERCISE 1 ---

## cleaning the workspace
rm(list=ls(all=TRUE))

## set up the directory where .R is saved (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading useful packages
library(tidyverse)

## loading data
load("data/FertSWE.Rdata")
range(FERT.SWE$Year)
range(FERT.SWE$Age)

## subsetting the data
my.df <- FERT.SWE %>% 
  filter(Age == 20, Year >= 1950)
range(my.df$Year)
range(my.df$Age)

## plotting the data
my.df %>% 
  ggplot(aes(x=Year,y=Births))+
  geom_point()

## EXERCISE 1 ------

## fitting a linear model
mod1 <- lm(formula=Births~Year,data=my.df)
summary(mod1)

## dimensions of the problem
t <- my.df$Year
y <- my.df$Births
plot(t,y)

## forecast years
tF <- min(my.df$Year):2050
my.df.fore <- tibble(Year=tF)

## fit and forecast of the model

## option 1
yF <- predict(object = mod1,newdata = my.df.fore)
plot(tF,yF)

plot(t,y,xlim=range(tF),ylim=range(y,yF),pch=16)
lines(tF,yF,col=2,lwd=2)
abline(h=0)

## option 2
yF2 <- coef(mod1)[1]+coef(mod1)[2]*tF

## checking option 1 == option 2
all.equal(unname(yF),yF2)

plot(yF2-yF)
plot(tF,yF,pch=16)
points(tF,yF2,col=2,pch=4,lwd=2)

## EXERCISE 2 ------

## dimensions of the problem
f <- my.df$Rates
plot(t,f)

## LM for fertility rates
mod2 <- lm(formula=Rates~Year,data=my.df)
summary(mod2)

## forecasting with LM
fF.LM <- coef(mod2)[1]+coef(mod2)[2]*tF

plot(t,f,xlim=range(tF),ylim=range(f,fF.LM),pch=16)
lines(tF,fF.LM,col=2,lwd=2)
abline(h=0)

## GLM for births with exposures
mod3 <- glm(formula=round(Births)~Year,data=my.df,
            family = poisson(),offset = log(Exposures))
summary(mod3)

## forecasting with GLM
fF.GLM <- exp(coef(mod3)[1]+coef(mod3)[2]*tF)

## this won't work as we do not have future exposures 
predict(object = mod3,newdata = my.df.fore)

## plotting
plot(t,f,xlim=range(tF),ylim=range(f,fF.LM,fF.GLM),
     pch=16)
lines(tF,fF.LM,col=2,lwd=2)
lines(tF,fF.GLM,col=4,lwd=2)
abline(h=0)

## EXTRA -----------

## GLM for births without exposures
mod4 <- glm(formula=round(Births)~Year,data=my.df,
            family = poisson())
summary(mod4)

## now this works (no need of exposures)
exp(predict(object = mod4,newdata = my.df.fore))


## END


