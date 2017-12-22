##### Develop Minimum Working Example of seasonal random effects in INLA

library(ggplot2)
library(INLA)
library(dplyr)
library(scales)
library(lubridate)

# generate some data with known characteristics
startdate=as.Date("2000-1-1")
stopdate=as.Date("2012-12-31")
data$doy=as.numeric(format(data$date,format="%j"))
data$doydate=as.Date(paste("1",data$doy,"2017",sep="-"),"%w-%j-%Y")
data$week=as.numeric(as.character(format(data$date,"%W")))
# make a covariate x
data$x=rnorm(n,sd=1)
# add a trend term
data$time=1:n
# add some seasonality
data$seasonality=sin(data$doy/365*2*pi+5)
# generate the response variable as the sum of the seasonality, x*beta, a trend, and some noise.
data$y=3.4+data$seasonality+data$time*(-0.02)+(data$x*1.5)+rnorm(n,mean=0,sd=1)


# plot the raw data
ggplot(data)+
  geom_line(aes(y=y,x=date),size=1)+
  geom_line(aes(y=x,x=date),col="green")+
  geom_line(aes(y=seasonality,x=date),col="red")
  #geom_line(aes(y=beta,x=date),col="blue")

# Default parameterization
formula=y~ x + time +
  f(week, model="seasonal", season.length=53)

# Try squeezing the prior on seasonal random effects to make them smoother. Squeezing precision just squeezes out the seasonal fluctuation -- not what we want! 
formula=y~ x + time +
  f(week, model="seasonal", season.length=53, hyper=list(prec = list(prior="loggamma",param=c(10,10))))

# fit the model
result=inla(formula,data=data,family="gaussian")

summary(result)

cf=result$summary.fixed
cf$var=rownames(cf)
cf=cf%>%
  mutate(sig=ifelse((`0.025quant`<0 &
                       `0.975quant`<0)|
                      `0.025quant`>0&
                      `0.975quant`>0,1,0))

ggplot(cf,aes(y=mean,x=var,
              ymax=`0.975quant`,
              ymin=`0.025quant`,
              color=as.factor(sig)))+
  geom_pointrange()+
  coord_flip()

## Extract seasonality
res_week=select(result$summary.random$week,
                week=ID,
                week.mean=mean,week.025=`0.025quant`,week.975=`0.975quant`)%>%
  mutate(weekdate=as.Date(paste("1",sprintf("%02d",week),"2017",sep="-"),"%w-%W-%Y"))

## Plot seasonality
ggplot(res_week,aes(x=weekdate))+
  geom_ribbon(aes(ymin=week.025,ymax=week.975),fill="red",
              alpha=0.2,col="transparent")+
  geom_line(aes(y=week.mean),col="darkred")+
  geom_line(aes(y=seasonality,x=doydate),col="darkblue", data=data)+
  geom_line(aes(y=y,x=doydate),col="green", data=data)+
  ylab("Seasonal Variation")+
  xlab("Date")+
  scale_x_date(labels = date_format("%B"),date_breaks = "2 month",date_minor_breaks = "1 month")


