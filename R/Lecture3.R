## -------------------------------------------------- ##
## EDSD 2024-2025: Demographic Forecasting
## Lecture 3
## Forecasting by parametric approaches
## Date: 06/05/2025
## Instructor: Ugofilippo Basellini
## -------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading useful packages
library(tidyverse)
library(forecast)
library(viridis)

## loading data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("data/FertSWE.Rdata")

## subset of my data & create Age^2
my.df <- FERT.SWE %>% filter(Year >= 1950) %>% 
  mutate(AgeSq = Age^2)

## subset again
my.df.sub <- my.df %>% 
  filter(Year==2000)

my.df.sub %>% 
  ggplot(aes(x=Age,y=logRates))+
  geom_point()

## extract objects
x <- my.df.sub$Age
m <- length(x)
x.sq <- x^2
lf <- my.df.sub$logRates
plot(x,lf)

## GLM for births with exposures
mod1 <- glm(formula=round(Births)~Age+AgeSq,data=my.df.sub,
            family = poisson(),offset = log(Exposures))
summary(mod1)

## fitted values on the log-scale with GLM
lfF.GLM <- coef(mod1)[1]+coef(mod1)[2]*x+coef(mod1)[3]*x.sq

plot(x,lf)
lines(x,lfF.GLM,col=4,lwd=2)

## number of years in the dataset
t <- unique(my.df$Year)
n <- length(t)

## repeat this for every year, and save the estimated coefs

## empty matrix to store the values of the coefficients
COEFS <- matrix(NA,n,3)
ETA.hat <- matrix(NA,m,n)

## for loop 
i <- 2
for (i in 1:n){
  ## subset the data
  my.df.sub <- my.df %>% 
    filter(Year==t[i])
  
  ## GLM for births with exposures
  mod1 <- glm(formula=round(Births)~Age+AgeSq,data=my.df.sub,
              family = poisson(),offset = log(Exposures))
  
  ## save coefficients
  COEFS[i,] <- coef(mod1)
  
  ## fitted values on the log-scale with GLM
  ETA.hat[,i] <- coef(mod1)[1]+coef(mod1)[2]*x+coef(mod1)[3]*x.sq
  
  
  
}

par(mfrow=c(1,3))
plot(t,COEFS[,1],t="p")
plot(t,COEFS[,2],t="p")
plot(t,COEFS[,3],t="p")
par(mfrow=c(1,1))

LRATES <- matrix(my.df$logRates,m,n) 

par(mfrow=c(1,2))
matplot(x,LRATES,t="l",lty=1,col=viridis(n))
matplot(x,ETA.hat,t="l",lty=1,col=viridis(n))
par(mfrow=c(1,1))

my.year <- 2010
whi.year <- which(t==my.year)
plot(x,LRATES[,whi.year],t="p",pch=1,col=viridis(n))
lines(x,ETA.hat[,whi.year],lwd=2,col=viridis(n))






