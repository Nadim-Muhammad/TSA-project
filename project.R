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
library(kableExtra)
library(formattable)
library(Metrics)
library(TSEwgt)
library(tidyverse)

###############################################################
#1. Getting data


setwd("C:\\Users\\nm412083\\Desktop\\Time Series Analysis\\project\\TSA-project-main\\TSA-project")



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

cointegrated_series <- lag.xts(residuals(cointegration))

chosen_ones$lag_resid_coint <- cointegrated_series

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
names(data_all_xts) <- "first"

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
names(data_all_second_xts) <- "second"

ARIMA_forecasts_second <- merge(data_all_second_xts,
                               arima1513_fixed_second_forecast_data_xts,
                               arima1513_fixed2_second_forecast_data_xts,
                               arima313_second_forecast_data_xts)
names(ARIMA_forecasts_second) <- c("Second Time Series", "arima1513_fixed mean", "arima1513_fixed 95% lower", "arima1513_fixed 95% upper", "arima1513_fixed2 mean", "arima1513_fixed2 95% lower", "arima1513_fixed2 95% upper", "arima313 mean", "arima313 95% lower", "arima313 95% upper")

plot(ARIMA_forecasts_second["2020-11/",], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of chosen time series",
     col = c("black", "blue", "red", "red", "green", "pink", "pink", "cyan", "magenta", "magenta"),
     legend.loc = "bottomleft")






######################################
#ECM and VECM

ecm_model <- lm(diff_first ~ diff_second + lag_resid_coint, data = chosen_ones)

#lets delete the intercept

ecm_model <- lm(diff_first ~ diff_second + lag_resid_coint -1, data = chosen_ones)

summary(ecm_model)
#diff_second is not statistically significant, but lets leave it


# The parameter -0.03568 describes a short term relationship 
# between `first` and `second`, so if `second` increases by 1 then 
# the `first` in the **short run** will decrease by 0.03568.

# The long run relationship is described by the parameter 
# 0.271123 from the cointegrating relationship: so if `ppi` 
# increases by 1 in the LONG RUN `cpi` will increase by 0.271123.

# The value of -0.55365 is the estimate of the *adjustment 
# coefficient*. As expected, its sign is negative and 
# this value means that about 55.37% of the unexpected error 
# (increase in gap) will be corrected in the next period, 
# so any unexpected deviation should be corrected finally 
# on average within about 1.8 periods.

#GRANGER CAUSALITY TEST

grangertest(first ~ second,
            data = chosen_ones,
            order = 3)
#second causes first (we rejected null hypothesis about no causality)

grangertest(second ~ first,
            data = chosen_ones,
            order = 3)
#first causes second (we rejected null hypothesis about no causality)


for(i in 1:5){
  print(i)
  print("first ~ second")
  print(grangertest(first ~ second,
                    data = chosen_ones,
                    order = i))
  print(i)
  print("second ~ first")
  print(grangertest(second ~ first,
              data = chosen_ones,
              order = i))
}
#Conclusion - we have bi-directional feedback in all cases

#VAR MODELS

#Lets select VAR model with proper lag legnht
VARselect(chosen_ones[,1:2], lag.max = 10) %>%
  .$criteria %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(nLags = 1:nrow(.)) %>%
  kbl(digits = 3) %>%
  kable_classic("striped", full_width = F)
#6 lags are the best according to all information criterias

#Lets check IC with seasonal VAR

VARselect(chosen_ones[,1:2], lag.max = 10, season = 31) %>%
  .$criteria %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(nLags = 1:nrow(.)) %>%
  kbl(digits = 3) %>%
  kable_classic("striped", full_width = F)
#All but one (FPE IC) shows that 6 lags are the best
#With ACF and PACF we showed that there is no seasonality detected

#But still let's check it (lets check season = 31, 30, 7:
VAR_model_6_lags <- VAR(chosen_ones[,1:2],
                        p = 6,
                        season = 31)

summary(VAR_model_6_lags)
#only sd10 in first is on the line of significance

VAR_model_6_lags <- VAR(chosen_ones[,1:2],
                        p = 6,
                        season = 30)

summary(VAR_model_6_lags)
#only sd12 in first is on the line of significance

VAR_model_6_lags <- VAR(chosen_ones[,1:2],
                        p = 6,
                        season = 7)

summary(VAR_model_6_lags)
#only sd3 in first is statistically significant

VAR_model_6_lags <- VAR(chosen_ones[,1:2],
                        p = 6)

summary(VAR_model_6_lags)
#some variables are not significant, but in general all lag orders have some significant part

plot(VAR_model_6_lags)
#ACF and PACF for first looks good, but for second both ACF and PACF have lag 8 on the line so lets make a model with 8 lags
VAR_model_8_lags <- VAR(chosen_ones[,1:2],
                        p = 8)

summary(VAR_model_8_lags)
#lags 7 and 8 are not significant

plot(VAR_model_8_lags)
#ACF and PACF in both cases looks good

#lets look at the residuals
serial.test(VAR_model_6_lags)

serial.test(VAR_model_8_lags)
#Both models show no autocorrelation in residuals


#AIC/BIC
AIC(VAR_model_6_lags, VAR_model_8_lags)
BIC(VAR_model_6_lags, VAR_model_8_lags)

#AIC prefers model with 8 lags, but BIC prefers model with 6 lags, lets stick to 6 lags, because as we showed earlier adding more lags gives us nonsignificant variables

#IRF
plot(irf(VAR_model_6_lags, n.ahead = 36))

#FEVD
plot(fevd(VAR_model_6_lags, n.ahead = 36))

#JOHANSON COINTEGRATION TEST

johan.test.trace <- 
  ca.jo(chosen_ones[,1:2],
        ecdet = "trend",
        type = "trace",  # type of the test: trace or eigen
        K = 6) # lags in VAR model
summary(johan.test.trace)

cbind(summary(johan.test.trace)@teststat, summary(johan.test.trace)@cval)
#exactly one cointegrated vector

#lets change type to eigen
johan.test.eigen <- 
  ca.jo(chosen_ones[,1:2],
        ecdet = "trend",
        type = "eigen",  # type of the test: trace or eigen
        K = 6) # lags in VAR model
summary(johan.test.eigen)

cbind(summary(johan.test.eigen)@teststat, summary(johan.test.eigen)@cval)
#same interpretation

#VECM

VECM_model <- cajorls(johan.test.trace, # defined specification
                                     r = 1) # number of cointegrating vectors

summary(VECM_model$rlm)
#only etc1 not significant in first equation

VECM_model$beta

# We can reparametrize the VEC model into VAR 
# (here we use the specification object):

VECM_model.asVAR <- vec2var(johan.test.trace, r = 1)

# Lets see the result:

VECM_model.asVAR

# Based on the reparametrized model, we can calculate 
# and plot Impulse Response Functions:

plot(irf(VECM_model.asVAR, n.ahead = 36))

# We can also perform variance decomposition:

plot(fevd(VECM_model.asVAR, n.ahead = 36))

# The results are pretty similar to the earlier `VAR(6)` 
# model.

# Let's also check if model residuals are autocorrelated.

# Residuals can be extracted only from the VAR 
# reparametrized model.

head(residuals(VECM_model.asVAR))
serial.test(VECM_model.asVAR)

# The null is not rejected, residuals are not autocorrelated.

# You can see the ACF and PACF functions by plotting 
# the results of the `serial.test()`

plot(serial.test(VECM_model.asVAR))
#8 lag in PACF of residuals and squared Residuals in second model seems significant

#lets plot EDF (Empirical Distribution Function) better

VECM_model.asVAR %>%
  residuals() %>%
  as_tibble() %>%
  ggplot(aes(`resids of first`)) +
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "pink") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(residuals(VECM_model.asVAR)[, 1]), 
                            sd = sd(residuals(VECM_model.asVAR)[, 1]))) +
  theme_bw() + 
  labs(
    title = "Density of PPI residuals", 
    y = "", x = "",
    caption = "source: own calculations"
  )

VECM_model.asVAR %>%
  residuals() %>%
  as_tibble() %>%
  ggplot(aes(`resids of second`)) +
  geom_histogram(aes(y =..density..),
                 colour = "black", 
                 fill = "pink") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(residuals(VECM_model.asVAR)[, 2]), 
                            sd = sd(residuals(VECM_model.asVAR)[, 2]))) +
  theme_bw() + 
  labs(
    title = "Density of CPI residuals", 
    y = "", x = "",
    caption = "source: own calculations"
  )


#lets use Jarque-Bera multivariate test to check for normality in residuals
normality.test(VECM_model.asVAR)
#null about normality is not rejected!!!

VECM_model.asVAR.fore <- 
  predict(
    VECM_model.asVAR,     # no of cointegrating vectors 
    n.ahead = 30, # forecast horizon
    ci = 0.95)


# VECM forecasts for `first`:

VECM_model.asVAR.fore$fcst$first

# VECM forecasts for `second`:

VECM_model.asVAR.fore$fcst$second


VECM_first_forecast <- xts(VECM_model.asVAR.fore$fcst$first[,-4], order.by = data_test$X)

VECM_second_forecast <- xts(VECM_model.asVAR.fore$fcst$second[,-4], order.by = data_test$X)


names(VECM_first_forecast) <- c("first_fore", "first_lower", "first_upper")

names(VECM_second_forecast) <- c("second_fore", "second_lower", "second_upper")


VECM_first_all_data <- merge(data_all_xts,
                             data_all_second_xts,
                             VECM_first_forecast,
                             VECM_second_forecast)

plot(VECM_first_all_data["2020-11/", c("first", "first_fore",
                        "first_lower", "first_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 days forecast of first",
     col = c("black", "blue", "red", "red"))

plot(VECM_first_all_data["2020-11/", c("second", "second_fore",
                        "second_lower", "second_upper")], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 days forecast of second",
     col = c("black", "blue", "red", "red"))


#WITH ARIMA MODELS

forecasts_first_all <- merge(ARIMA_forecasts_first,
                             VECM_first_forecast)

forecasts_second_all <- merge(ARIMA_forecasts_second,
                             VECM_second_forecast)



#first
plot(forecasts_first_all["2020-11/",], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of chosen time series",
     col = c("black", "blue", "red", "red", "green", "pink", "pink", "cyan", "magenta", "magenta", "orange", "yellow", "yellow"),
     legend.loc = "bottomleft")

#second
plot(forecasts_second_all["2020-11/",], 
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "30 day forecast of chosen time series",
     col = c("black", "blue", "red", "red", "green", "pink", "pink", "cyan", "magenta", "magenta", "orange", "yellow", "yellow"),
     legend.loc = "bottomleft")

#EX POST ERRORS

names(forecasts_first_all)

names(forecasts_second_all)

################# FIRST
errors_first_ts <- {}


#ARIMA(6,1,3)
errors_first_ts$arima613mae   <-  mae(forecasts_first_all[971:1000,1],forecasts_first_all$arima613.mean)
errors_first_ts$arima613mse   <-  mse(forecasts_first_all[971:1000,1],forecasts_first_all$arima613.mean)
errors_first_ts$arima613mape  <-  mape(forecasts_first_all[971:1000,1],forecasts_first_all$arima613.mean)
errors_first_ts$arima613amape <-  mean(abs((forecasts_first_all[971:1000,1] - forecasts_first_all$arima613.mean) / 
                                        (forecasts_first_all[971:1000,1] + forecasts_first_all$arima613.mean)))

#ARIMA(7,1,3)
errors_first_ts$arima713mae   <-  mae(forecasts_first_all[971:1000,1],forecasts_first_all$arima713.mean)
errors_first_ts$arima713mse   <-  mse(forecasts_first_all[971:1000,1],forecasts_first_all$arima713.mean)
errors_first_ts$arima713mape  <-  mape(forecasts_first_all[971:1000,1],forecasts_first_all$arima713.mean)
errors_first_ts$arima713amape <-  mean(abs((forecasts_first_all[971:1000,1] - forecasts_first_all$arima713.mean) / 
                                             (forecasts_first_all[971:1000,1] + forecasts_first_all$arima713.mean)))

#ARIMA(6,1,4)
errors_first_ts$arima614mae   <-  mae(forecasts_first_all[971:1000,1],forecasts_first_all$arima614.mean)
errors_first_ts$arima614mse   <-  mse(forecasts_first_all[971:1000,1],forecasts_first_all$arima614.mean)
errors_first_ts$arima614mape  <-  mape(forecasts_first_all[971:1000,1],forecasts_first_all$arima614.mean)
errors_first_ts$arima614amape <-  mean(abs((forecasts_first_all[971:1000,1] - forecasts_first_all$arima614.mean) / 
                                             (forecasts_first_all[971:1000,1] + forecasts_first_all$arima614.mean)))

#VECM
errors_first_ts$VECM_mae   <-  mae(forecasts_first_all[971:1000,1],forecasts_first_all$first_fore)
errors_first_ts$VECM_mse   <-  mse(forecasts_first_all[971:1000,1],forecasts_first_all$first_fore)
errors_first_ts$VECM_mape  <-  mape(forecasts_first_all[971:1000,1],forecasts_first_all$first_fore)
errors_first_ts$VECM_amape <-  mean(abs((forecasts_first_all[971:1000,1] - forecasts_first_all$first_fore) / 
                                             (forecasts_first_all[971:1000,1] + forecasts_first_all$first_fore)))

errors_first_ts

#MAE
errors_first_ts$arima613mae
errors_first_ts$arima713mae
errors_first_ts$arima614mae
errors_first_ts$VECM_mae
#arima(6,1,4) the best in terms of MAE
?cbind
#MSE
errors_first_ts$arima613mse
errors_first_ts$arima713mse
errors_first_ts$arima614mse
errors_first_ts$VECM_mse
#arima(6,1,4) the best in terms of MSE

#MAPE
errors_first_ts$arima613mape
errors_first_ts$arima713mape
errors_first_ts$arima614mape
errors_first_ts$VECM_mape
#arima(6,1,4) the best in terms of MAPE

#AMAPE
errors_first_ts$arima613amape
errors_first_ts$arima713amape
errors_first_ts$arima614amape
errors_first_ts$VECM_amape
#arima(6,1,4) the best in terms of AMAPE

#################################################################

################# SECOND
errors_second_ts <- {}


#ARIMA(15,1,3) fixed
errors_second_ts$arima1513_fixed_mae   <-  mae(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed.mean)
errors_second_ts$arima1513_fixed_mse   <-  mse(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed.mean)
errors_second_ts$arima1513_fixed_mape  <-  mape(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed.mean)
errors_second_ts$arima1513_fixed_amape <-  mean(abs((forecasts_second_all[971:1000,1] - forecasts_second_all$arima1513_fixed.mean) / 
                                             (forecasts_second_all[971:1000,1] + forecasts_second_all$arima1513_fixed.mean)))


#ARIMA(15,1,3) fixed2
errors_second_ts$arima1513_fixed2_mae   <-  mae(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed2.mean)
errors_second_ts$arima1513_fixed2_mse   <-  mse(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed2.mean)
errors_second_ts$arima1513_fixed2_mape  <-  mape(forecasts_second_all[971:1000,1],forecasts_second_all$arima1513_fixed2.mean)
errors_second_ts$arima1513_fixed2_amape <-  mean(abs((forecasts_second_all[971:1000,1] - forecasts_second_all$arima1513_fixed2.mean) / 
                                                     (forecasts_second_all[971:1000,1] + forecasts_second_all$arima1513_fixed2.mean)))


#ARIMA(3,1,3)
errors_second_ts$arima313_mae   <-  mae(forecasts_second_all[971:1000,1],forecasts_second_all$arima313.mean)
errors_second_ts$arima313_mse   <-  mse(forecasts_second_all[971:1000,1],forecasts_second_all$arima313.mean)
errors_second_ts$arima313_mape  <-  mape(forecasts_second_all[971:1000,1],forecasts_second_all$arima313.mean)
errors_second_ts$arima313_amape <-  mean(abs((forecasts_second_all[971:1000,1] - forecasts_second_all$arima313.mean) / 
                                                     (forecasts_second_all[971:1000,1] + forecasts_second_all$arima313.mean)))


#VECM
errors_second_ts$VECM_mae   <-  mae(forecasts_second_all[971:1000,1],forecasts_second_all$second_fore)
errors_second_ts$VECM_mse   <-  mse(forecasts_second_all[971:1000,1],forecasts_second_all$second_fore)
errors_second_ts$VECM_mape  <-  mape(forecasts_second_all[971:1000,1],forecasts_second_all$second_fore)
errors_second_ts$VECM_amape <-  mean(abs((forecasts_second_all[971:1000,1] - forecasts_second_all$second_fore) / 
                                              (forecasts_second_all[971:1000,1] + forecasts_second_all$second_fore)))


#MAE
errors_second_ts$arima1513_fixed_mae
errors_second_ts$arima1513_fixed2_mae
errors_second_ts$arima313_mae
errors_second_ts$VECM_mae
#ARIMA(15,1,3) fixed best in terms of MAE

#MSE
errors_second_ts$arima1513_fixed_mse
errors_second_ts$arima1513_fixed2_mse
errors_second_ts$arima313_mse
errors_second_ts$VECM_mse
#ARIMA(15,1,3) fixed best in terms of MSE

#MAPE
errors_second_ts$arima1513_fixed_mape
errors_second_ts$arima1513_fixed2_mape
errors_second_ts$arima313_mape
errors_second_ts$VECM_mape
#ARIMA(15,1,3) fixed best in terms of MAPE

#AMAPE
errors_second_ts$arima1513_fixed_amape
errors_second_ts$arima1513_fixed2_amape
errors_second_ts$arima313_amape
errors_second_ts$VECM_amape
#ARIMA(15,1,3) fixed best in terms of AMAPE

