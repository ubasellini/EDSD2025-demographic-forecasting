## -------------------------------------------------- ##
## EDSD 2024-2025: Demographic Forecasting
## Lecture 2
## Direct extrapolation by time-series methods
## Date: 06/05/2024
## Instructor: Ugofilippo Basellini
## -------------------------------------------------- ##

## EXERCISE 1 -------

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading useful packages
library(tidyverse)
library(forecast)
library(tseries)

## loading data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("data/TimeSeries.Rdata")

## inspect the data
head(sp500.data)

## subset the data
my.df <- sp500.data %>% 
  filter(Date>="1990-01-01") %>% 
  slice(1:250)

## dimensions of the problem
y.all <- my.df$SP500
t.all <- my.df$Date
y <- y.all[1:200]
t <- t.all[1:200]

## plotting
plot(t.all,y.all)
points(t,y,pch=16)

## ACF
acf(y)  ## y is non-stationary

## taking first differences
y.diff <- diff(y)
plot(t[-1],y.diff,t="l")

## ACF
acf(y.diff)  ## y is stationary

## BONUS - compute the KPSS and ADF tests
## H0: stationary
tseries::kpss.test(y)      ## reject H0
tseries::kpss.test(y.diff) ## cannot reject H0

## H0: non-stationary
tseries::adf.test(y)      ## cannot reject H0
tseries::adf.test(y.diff) ## reject H0


## fitting RW models
mod1 <- Arima(y,order = c(0,1,0))  ## no drift
mod2 <- Arima(y,order = c(0,1,0), include.drift = T)  ## with drift
summary(mod1)
summary(mod2)


## testing if getting same results by fitting 
## ARIMA on differenced series 
## (also look at https://stats.stackexchange.com/questions/32634/difference-time-series-before-arima-or-within-arima )
mod.test <- Arima(y.diff,order = c(0,0,0),
                  include.mean = T,
                  include.constant = F)
summary(mod.test)
yF3 <- forecast(mod.test,h=h)

## define lenght of forecast horizon
n.all <- length(t.all)
n <- length(t)
h <- n.all-n

## forecast period
tF <- t.all[!t.all%in%t]

## producing fort.all## producing forecast
yF1 <- forecast(mod1,h=h)
yF2 <- forecast(mod2,h=h)
plot(yF1)
plot(yF2)

## comparing with observed data
plot(t.all,y.all,ylim = range(y.all,yF2$mean))
points(t,y,pch=16)
lines(tF,yF1$mean,col=2,lwd=2)
lines(tF,yF2$mean,col=4,lwd=2)

plot(t.all[-1],diff(y.all),
     ylim = range(diff(y.all),diff(yF2$mean)))
points(t[-1],diff(y),pch=16)
lines(tF[-1],diff(yF1$mean),col=2,lwd=3)
lines(tF[-1],diff(yF2$mean),col=4,lwd=2)
lines(tF,yF3$mean,col=5,lwd=2,lty=2)

## EXERCISE 2 -------

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("data/FertSWE.Rdata")

## subsetting the data
my.df <- FERT.SWE %>% 
  filter(Age == 20, Year >= 1950)

## dimensions of the problem
t <- my.df$Year
y <- my.df$Births
e <- my.df$Exposures
f <- y/e
lf <- log(f)
plot(t,lf)

## forecast years
t.all <- min(my.df$Year):2050
tF <- t.all[!t.all%in%t]
h <- length(tF)

## GLM for births with exposures
mod1 <- glm(formula=round(Births)~Year,data=my.df,
            family = poisson(),offset = log(Exposures))

## fitting/forecasting with GLM
lfF.GLM <- coef(mod1)[1]+coef(mod1)[2]*t.all

## plotting
plot(t,lf,xlim = range(t.all),ylim=range(lfF.GLM,lf))
lines(t.all,lfF.GLM,col=2,lwd=2)

## Time-series models

## stationarity check
acf(lf)
acf(diff(lf))

## BONUS - compute the KPSS and ADF tests
## H0: stationary
tseries::kpss.test(lf)      ## reject H0
tseries::kpss.test(diff(lf)) ## cannot reject H0

## H0: non-stationary
tseries::adf.test(lf)      ## cannot reject H0
tseries::adf.test(diff(lf)) ## reject H0

## RWD 
mod2 <- Arima(lf,order = c(0,1,0), include.drift = T)  ## with drift

## forecast
lfF.RWD <- forecast(mod2,h=h)
plot(lfF.RWD)

## best arima model
mod3 <- auto.arima(lf)
summary(mod3)
lfF.AR <- forecast(mod3,h=h)
plot(lfF.AR)

## plotting
plot(t,lf,xlim = range(t.all),
     ylim=range(lfF.GLM,lf,lfF.AR$lower[,1]))
lines(t.all,lfF.GLM,col=2,lwd=2)
## RWD
lines(tF,lfF.RWD$mean,col=5,lwd=2)
lines(tF,lfF.RWD$upper[,1],col=5,lwd=2,lty=2)
lines(tF,lfF.RWD$lower[,1],col=5,lwd=2,lty=2)
## AR
lines(tF,lfF.AR$mean,col=6,lwd=2)
lines(tF,lfF.AR$upper[,1],col=6,lwd=2,lty=2)
lines(tF,lfF.AR$lower[,1],col=6,lwd=2,lty=2)


