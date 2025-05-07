## -------------------------------------------------- ##
## EDSD 2024-2025: Demographic Forecasting
## Lecture 4
## Forecasting by Lee-Carter method
## Date: 07/05/2025
## Instructor: Ugofilippo Basellini
## -------------------------------------------------- ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading useful packages
library(tidyverse)
library(forecast)
library(viridis)
library(fields)

## loading data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("data/MORTSWE.Rdata")

## subset the data
my.df <- MORT.SWE %>% 
  filter(Year>=1950,Sex=="Male",Age<=100)

## dimensions of the problem
x <- unique(my.df$Age)
t <- unique(my.df$Year)
m <- length(x)
n <- length(t)

## extract matrices from the dataset
Y <- matrix(my.df$Deaths,m,n)
E <- matrix(my.df$Exposures,m,n)
MX <- Y/E
LMX <- log(MX)
image(t,x,t(LMX),col=viridis(n))

## matrix of deaths without 0 & log-rates without inf
Y1 <- Y
Y1[Y==0] <- 1
LMX1 <- log(Y1/E)
image(t,x,t(LMX1),col=viridis(n))

## Alpha
Alpha <- apply(LMX1,1,mean)
matplot(x,LMX1,t="l",col="grey80",lty=1)
lines(x,Alpha,col="darkgreen",lwd=2)

## Beta and Kappa

## centred matrix of log-rates
LMX.cent <- LMX1-Alpha
image.plot(t,x,t(LMX1),col=viridis(n))
image.plot(t,x,t(LMX.cent),col=viridis(n))

## svd
my.svd <- svd(LMX.cent)
Beta <- my.svd$u[,1]
Kappa <- my.svd$v[,1]
par(mfrow=c(1,2))
plot(x,Beta)
plot(t,Kappa)

LMX.cent.svd <- my.svd$d[1]*(Beta%*%t(Kappa))
LMX.cent.svd3 <- my.svd$d[1]*(Beta%*%t(Kappa)) +
  my.svd$d[2]*(my.svd$u[,2]%*%t(my.svd$v[,2])) +
  my.svd$d[3]*(my.svd$u[,3]%*%t(my.svd$v[,3])) 

par(mfrow=c(1,3))
image.plot(t,x,t(LMX.cent),col=viridis(n))
image.plot(t,x,t(LMX.cent.svd),col=viridis(n))
image.plot(t,x,t(LMX.cent.svd3),col=viridis(n))


## including the constraints
## constraint 1
sum.Beta <- sum(Beta)
Beta <- Beta/sum.Beta
sum(Beta)
## constraint 2
Kappa1 <- my.svd$d[1]*Kappa*sum.Beta
sum(Kappa1)
## plotting
plot(x,Beta)
plot(t,Kappa1)
par(mfrow=c(1,1))



##---- step 3: adjust KAPPA -----
## function to compute difference between observed and fitted LC deaths
koptim <- function(par,alpha,beta,sum.dx,Exp){
  kappa <- par[1]
  lmx.lc <- alpha+beta*kappa
  z.lc <- exp(lmx.lc)*Exp
  sum.z.lc <- sum(z.lc)
  diff.lc <- abs(sum.dx-sum.z.lc)
  return(diff.lc)
}
## adjust Kappa every year
Kappa <- numeric(n)
for (i in 1:n){
  KappaSecStep <- optimize(f=koptim,interval=c(-100,100),alpha=Alpha,
                           beta=Beta,sum.dx=sum(Y[,i]),Exp=E[,i])
  Kappa[i] <- KappaSecStep$minimum
}
## plotting
plot(t,Kappa1,ylim=range(Kappa1,Kappa),t="l",lwd=2)
lines(t,Kappa,col=4,lwd=2)
legend("topright",c("from SVD","second-step adjustment"),col=c(1,4),lwd=2)


## fitted log-mortality
Ones <- matrix(1,n)
ETAlc <- Alpha%*%t(Ones) + Beta%*%t(Kappa)
## basic plot
g <- my.df %>%
  mutate(logRates=case_when(
    is.infinite(logRates)~NA,
    TRUE~logRates),
    Fitted=c(ETAlc)) %>%
  ggplot(aes(x=Age,group=Year))+
  geom_point(aes(y=logRates))+
  geom_line(aes(y=Fitted),color="darkorange",linewidth=1.2)+
  scale_color_viridis_c() +
  theme_bw(base_size = 18) +
  labs(y= "Log Mortality Rate")
## animating with gganimate
library(gganimate)
gg <- g + transition_time(Year) +
  ggtitle("Year {frame_time}")
animate(gg, fps=4)


## END



