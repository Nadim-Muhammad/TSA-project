library(xts)
library(lmtest)
library(quantmod)
library(dplyr)
library(fUnitRoots)
library(vars)
library(tseries)
library(aTSA)
library(car)
library(seasonal)
library(forecast)
library(tidyverse)

###############################################################
#1. Getting data


setwd("C:\\Users\\Administrator\\OneDrive\\Pulpit\\studia\\MAGISTERSKIE\\Time Series Analysis\\project\\TSA-project")



options(scipen = 10)

source("functions/function_plot_ACF_PACF_resids.R")
source("functions/testdf.R")




data <- read.csv("data/TSA_2023_project_data_1.csv")

data_train <- data[1:970,]

data_test <- data[971:1000,]

data_train$X <- as.Date(data_train$X, format = "%Y-%m-%d")

data_test$X <- as.Date(data_test$X, format = "%Y-%m-%d")

# Let's also transform the `data.frame` into an `xts` object

data_train_xts <- xts(data_train[, -1], order.by = data_train$X)

data_test_xts <- xts(data_test[, -1], order.by = data_test$X)

###############################################################
#2 Choosing time series

plot(data_train_xts,
     col = c("black", "blue", "red", "purple", "green", "pink", "brown", "cyan", "magenta", "coral"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "Plotted time series",
     legend.loc = "bottomleft")


#BG test - null hypothesis: no serial correlation of any order (up to max.augmentations)
#ADF test - null hypothesis: unit root in the time series


for (i in 1:10) {
  print(i)
  print(testdf(variable = data_train_xts[,i],
  max.augmentations = 3))
}
#-4

data_train_xts_diff <- diff.xts(data_train_xts)

for (i in 1:10) {
  print(i)
  print(testdf(variable = data_train_xts_diff[,i],
               max.augmentations = 3))
}
#1 3 4 9

chosen_ones <- xts(order.by = data_train$X)

chosen_ones$first<- data_train_xts[,3]
chosen_ones$second<- data_train_xts[,9]

plot(chosen_ones)

chosen_ones$diff_first <- diff.xts(chosen_ones$first)
chosen_ones$diff_second <- diff.xts(chosen_ones$second)


#first
testdf(variable = chosen_ones$first,
       max.augmentations = 3)

testdf(variable = chosen_ones$diff_first,
       max.augmentations = 3)

#second
testdf(variable = chosen_ones$second,
       max.augmentations = 3)

testdf(variable = chosen_ones$diff_second,
       max.augmentations = 3)

#Both series integrated of order one (I(1))

###############################################################
#2.1 Checking for cointegration relation

cointegration <- lm(first ~ second, data = chosen_ones)

summary(cointegration)

testdf(variable = residuals(cointegration), max.augmentations = 3)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 152.195497   0.107127  1420.7   <2e-16 ***
#  second       -0.271123   0.002145  -126.4   <2e-16 ***
# The cointegrating vector is: 1, -152.195497, 0.271123
# which defines the cointegrating relationship as: 
# residuals = 1 * first - 152.195497 + 0.271123 * second

cointegrated_series <- as.xts(residuals(cointegration))

plot(cointegrated_series)

testdf(variable = cointegrated_series, max.augmentations = 3)


#############################################################
#3 ARIMA models

###############################################################
#3.1 ARIMA for first time series (from before we know that this series is I(1))
par(mfrow = c(2, 1)) 
acf(chosen_ones$diff_first,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),   # limits for the y axis - we give c(min, max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
pacf(chosen_ones$diff_first, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1))
#looks like we should use ARIMA(7,1,x)

#Lets check ARIMA(1,1,1)

arima111 <- Arima(chosen_ones$first, order = c(1, 1, 1))

coeftest(arima111)

plot_ACF_PACF_resids(arima111)
#ACF and PACF plots shows that some lags are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima111), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are not white noise


#Lets try ARIMA(7,1,3)
arima713 <- Arima(chosen_ones$first, 
                  order = c(7, 1, 3), 
                  include.constant = TRUE,
                  optim.control = list(maxit = 800),
                  optim.method = "L-BFGS-B")

arima713

coeftest(arima713)
#only ar7 and ma1 are statistically insignificant (and drift)

plot_ACF_PACF_resids(arima713)
#lags 15, 21 still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima713), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise


#Lets check ARIMA(7,1,15)

arima7115 <- Arima(chosen_ones$first, order = c(7, 1, 15))

arima7115

coeftest(arima7115)
#ar3, ar5 and ma7 and ma15 are statistically significant

plot_ACF_PACF_resids(arima7115)
# lags 26 is on significance line

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima7115), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise

#Lets check ARIMA(7,1,26)

arima7126 <- Arima(chosen_ones$first, 
                    order = c(7, 1, 26),
                    include.constant = TRUE,
                    optim.control = list(maxit = 800),
                    optim.method = "L-BFGS-B")

coeftest(arima7126)
#ar1,ar2, ar4, ar5, ar6 and ma3 statistically significant

plot_ACF_PACF_resids(arima7126)
# lag 26 is significant, lets leave it for now

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima7126), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise


#Now lets use previous models but lets set nonsignificant variables to 0
#(arima713fixed but without last ar term)
arima613fixed <- Arima(chosen_ones$first,
                    order = c(6, 1, 3),
                    fixed = c(NA, NA, NA, NA, NA, NA,
                              0, NA, NA))

coeftest(arima613fixed)
#nothing expect ar1 is significant

plot_ACF_PACF_resids(arima613fixed)
#we still have some other lags significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima613fixed), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are not a white noise!!!

#lets try arima7115 with nonsignificant variables set to 0 (expect ma10)
#(its arima5115)
arima5115fixed <- Arima(chosen_ones$first,
                       order = c(5, 1, 15),
                       fixed = c(0, 0, NA, 0, NA,
                                 0,0,0,0,0,0,NA,0,0,NA,0,0,0,0,NA))                                 

coeftest(arima5115fixed)
#everything is significant

plot_ACF_PACF_resids(arima5115fixed)
#we still have some other lags significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima5115fixed), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are not white noise!!!!

#lets try arima7126 with nonsignificant variables set to 0
#(its arima614) lets include ma2 and ma4
arima614fixed <- Arima(chosen_ones$first,
                       order = c(6, 1, 4),
                       fixed = c(NA, NA, 0, NA, NA, NA,
                                 0,NA,NA,NA))                                 

coeftest(arima614fixed)
#ar4 and ma4 are not significant

plot_ACF_PACF_resids(arima614fixed)
#lag 15,20,21,26 on significance line

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima614fixed), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise


#lets set to 0 nonsignificant variables from previous model
#(its arima 613 now)
arima613fixed <- Arima(chosen_ones$first,
                        order = c(6, 1, 3),
                        fixed = c(NA, NA, 0, 0, NA, NA,
                                  0, NA, NA))                                 

coeftest(arima613fixed)
#everything significant

plot_ACF_PACF_resids(arima613fixed)
#lag 4, 15, 20, 21 and 26 on line

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima613fixed), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise

#lets add AR(4) to the last model

arima613fixed2 <- Arima(chosen_ones$first,
                       order = c(6, 1, 3),
                       fixed = c(NA, NA, 0, NA, NA, NA,
                                 0, NA, NA))                                 

coeftest(arima613fixed2)
#everything significant

plot_ACF_PACF_resids(arima613fixed2)
#lag 15, 20, 21 on line

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima613fixed2), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise

###############################################################
#4 INFORMATION CRITERIAS
AIC(arima111, arima713, arima7115, arima7126,
    arima613fixed, arima613fixed2, arima614fixed, arima5115fixed)

BIC(arima111, arima713, arima7115, arima7126,
    arima613fixed, arima613fixed2, arima614fixed, arima5115fixed)

arima613fixed2
coeftest(arima613fixed2)

#Lets forecast 3 models with the lowest AIC
#arima613fixed2, arima713, arima614fixed

arima613fixed2_forecast <- forecast(arima613fixed2, h=30)
arima613fixed2_forecast

arima613fixed2_forecast_data <- data.frame(f_mean  = as.numeric(arima613fixed2_forecast$mean),
                                           f_lower = as.numeric(arima613fixed2_forecast$lower[, 2]),
                                           f_upper = as.numeric(arima613fixed2_forecast$upper[, 2]))

arima613fixed2_forecast_data_xts <- xts(arima613fixed2_forecast_data, order.by = data_test$X)


arima713_forecast <- forecast(arima713, h=30)
arima713_forecast

arima713_forecast_data <- data.frame(f_mean  = as.numeric(arima713_forecast$mean),
                                     f_lower = as.numeric(arima713_forecast$lower[, 2]),
                                     f_upper = as.numeric(arima713_forecast$upper[, 2]))

arima713_forecast_data_xts <- xts(arima713_forecast_data, order.by = data_test$X)


arima614fixed_forecast <- forecast(arima614fixed, h=30)
arima614fixed_forecast

arima614fixed_forecast_data <- data.frame(f_mean  = as.numeric(arima614fixed_forecast$mean),
                                          f_lower = as.numeric(arima614fixed_forecast$lower[, 2]),
                                          f_upper = as.numeric(arima614fixed_forecast$upper[, 2]))

arima614fixed_forecast_data_xts <- xts(arima614fixed_forecast_data, order.by = data_test$X)

###############################################################
#5 Forecast visualization
data_all_xts <- rbind(data_train_xts[,3], data_test_xts[,3])

ARIMA_forecasts_first <- merge(data_all_xts,
                               arima613fixed2_forecast_data_xts,
                               arima713_forecast_data_xts,
                               arima614fixed_forecast_data_xts)
names(ARIMA_forecasts_first) <- c("Time Series", "arima613 mean", "arima613 95% lower", "arima613 95% upper", "arima713 mean", "arima713 95% lower", "arima713 95% upper", "arima614 mean", "arima614 95% lower", "arima614 95% upper")

plot(ARIMA_forecasts_first["2020-11/",], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of chosen time series",
     col = c("black", "blue", "red", "red", "green", "pink", "pink", "cyan", "magenta", "magenta"),
     legend.loc = "bottomleft")

###############################################################
#6 ARIMA for second time series (from before we know that this series is I(1))
par(mfrow = c(2, 1)) 
acf(chosen_ones$diff_second,
    lag.max = 36, # max lag for ACF
    ylim = c(-0.1, 0.1),   # limits for the y axis - we give c(min, max)
    lwd = 5,               # line width
    col = "dark green",
    na.action = na.pass)   # do not stop if there are missing values in the data
pacf(chosen_ones$diff_second, 
     lag.max = 36, 
     lwd = 5, col = "dark green",
     na.action = na.pass)
par(mfrow = c(1, 1))
#lets start with checking ARIMA(1,1,1) and ARIMA(3,1,3)

#ARIMA(1,1,1)
arima111_second <- Arima(chosen_ones$second, order = c(1, 1, 1))

coeftest(arima111_second)
#all variables significant

plot_ACF_PACF_resids(arima111_second)
#ACF and PACF plots shows that some lags are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima111_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are not white noise

#ARIMA(3,1,3)
arima313_second <- Arima(chosen_ones$second, order = c(3, 1, 3))

coeftest(arima313_second)
#ma3 not significant

plot_ACF_PACF_resids(arima313_second)
#lags 15, 18, 20 are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima313_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise




#ARIMA(15,1,3)
arima1513_second <- Arima(chosen_ones$second, order = c(15, 1, 3))

coeftest(arima1513_second)
#ar1, ar2, ar3, ar4, ar8, ar15, ma1, ma3 significant

plot_ACF_PACF_resids(arima1513_second)
#lags 17, 20, 26 are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima1513_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise



#Now lets check previous models, but without insignificant variables

#ARIMA(3,1,3) without ma3 = ARIMA(3,1,2)
arima312_second <- Arima(chosen_ones$second, order = c(3, 1, 2))

coeftest(arima312_second)
#everything significant

plot_ACF_PACF_resids(arima312_second)
#lags 15, 18, 20, 26 are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima312_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise (?)

#ARIMA(15,1,3), but without insignificant variables
arima1513_fixed_second <- Arima(chosen_ones$second,
                                order = c(15, 1, 3),
                                fixed = c(NA, NA, NA, NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA,
                                          NA, 0, NA))  

coeftest(arima1513_fixed_second)
#ar1, ar4, ma1 not significant

plot_ACF_PACF_resids(arima1513_fixed_second)
#lag 20 is still significant, lag 6 and 17 on line

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima1513_fixed_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise


#ARIMA(15,1,3) fixed, but without insignificant variables

arima1513_fixed2_second <- Arima(chosen_ones$second,
                                order = c(15, 1, 3),
                                fixed = c(0, NA, NA, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA,
                                          0, 0, NA))  

coeftest(arima1513_fixed2_second)
#everything significant

plot_ACF_PACF_resids(arima1513_fixed2_second)
#lags 1, 6, 20, 26 are still significant

bj_pvalues = c()
for(i in c(1:100)){
  bj = Box.test(resid(arima1513_fixed2_second), type = "Ljung-Box", lag = i)
  bj_pvalues = append(bj_pvalues,bj$p.value)
}
plot(bj_pvalues, type='l')
abline(h=0.05, col='red', type='l')
#residuals are white noise

################################################
#4 INFORMATION CRITERIAS
AIC(arima111_second, arima313_second, arima1513_second, arima1513_fixed_second, arima1513_fixed2_second)

BIC(arima111_second, arima313_second, arima1513_second, arima1513_fixed_second, arima1513_fixed2_second)

#lets choose 3 with the lowest AIC
#arima1513_fixed_second, arima1513_fixed2_second, arima313_second



arima1513_fixed_second_forecast <- forecast(arima1513_fixed_second, h=30)
arima1513_fixed_second_forecast

arima1513_fixed_second_forecast_data <- data.frame(f_mean  = as.numeric(arima1513_fixed_second_forecast$mean),
                                           f_lower = as.numeric(arima1513_fixed_second_forecast$lower[, 2]),
                                           f_upper = as.numeric(arima1513_fixed_second_forecast$upper[, 2]))

arima1513_fixed_second_forecast_data_xts <- xts(arima1513_fixed_second_forecast_data, order.by = data_test$X)


arima1513_fixed2_second_forecast <- forecast(arima1513_fixed2_second, h=30)
arima1513_fixed2_second_forecast

arima1513_fixed2_second_forecast_data <- data.frame(f_mean  = as.numeric(arima1513_fixed2_second_forecast$mean),
                                     f_lower = as.numeric(arima1513_fixed2_second_forecast$lower[, 2]),
                                     f_upper = as.numeric(arima1513_fixed2_second_forecast$upper[, 2]))

arima1513_fixed2_second_forecast_data_xts <- xts(arima1513_fixed2_second_forecast_data, order.by = data_test$X)


arima313_second_forecast <- forecast(arima313_second, h=30)
arima313_second_forecast

arima313_second_forecast_data <- data.frame(f_mean  = as.numeric(arima313_second_forecast$mean),
                                          f_lower = as.numeric(arima313_second_forecast$lower[, 2]),
                                          f_upper = as.numeric(arima313_second_forecast$upper[, 2]))

arima313_second_forecast_data_xts <- xts(arima313_second_forecast_data, order.by = data_test$X)

########################################################
#VISUALIZATION OF THE SECOND TIME SERIES
data_all_second_xts <- rbind(data_train_xts[,9], data_test_xts[,9])

ARIMA_forecasts_second <- merge(data_all_second_xts,
                               arima1513_fixed_second_forecast_data_xts,
                               arima1513_fixed2_second_forecast_data_xts,
                               arima313_second_forecast_data_xts)
names(ARIMA_forecasts_first) <- c("Second Time Series", "arima1513_fixed mean", "arima1513_fixed 95% lower", "arima1513_fixed 95% upper", "arima1513_fixed2 mean", "arima1513_fixed2 95% lower", "arima1513_fixed2 95% upper", "arima313 mean", "arima313 95% lower", "arima313 95% upper")

plot(ARIMA_forecasts_second["2020-11/",], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of chosen time series",
     col = c("black", "blue", "red", "red", "green", "pink", "pink", "cyan", "magenta", "magenta"),
     legend.loc = "bottomleft")






######################################
#WE NEED A VAR LAG ORDER FIRST
johan.test.trace <- 
  ca.jo(chosen_ones[,1:2], # data 
        ecdet = "trend", # "none" for no intercept in cointegrating equation, 
        # "const" for constant term in cointegrating equation and 
        # "trend" for trend variable in cointegrating equation
        type = "trace",  # type of the test: trace or eigen
        K = 2,           # lag order of the series (levels) in the VAR
        season = 12) 
summary(johan.test.trace)
?ca.jo

