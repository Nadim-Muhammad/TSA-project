library(xts)
library(forecast)
library(lmtest)
library(quantmod)
library(tidyverse)
library(dplyr)
library(fUnitRoots)
library(vars)
library(tseries)
library(aTSA)
library(car)
library(seasonal)

setwd("C:\\Users\\Administrator\\OneDrive\\Pulpit\\studia\\MAGISTERSKIE\\Time Series Analysis\\project")



options(scipen = 10)

source("functions/function_plot_ACF_PACF_resids.R")
source("functions/testdf.R")




data <- read.csv("data/TSA_2023_project_data_1.csv")

data$X <- as.Date(data$X, format = "%Y-%m-%d")

# Let's also transform the `data.frame` into an `xts` object

data_xts <- xts(data[, -1], order.by = data$X)

data_xts_in <- data_xts[1:970,]

data_xts_out <- data_xts[971:1000,]

plot(data_xts_in,
     col = c("black", "blue", "red", "purple", "green", "pink", "brown", "cyan", "magenta", "coral"),
     major.ticks = "years", 
     grid.ticks.on = "years",
     grid.ticks.lty = 3,
     main = "Something",
     legend.loc = "bottomleft")


#BG test - null hypothesis: no serial correlation of any order (up to max.augmentations)
#ADF test - null hypothesis: unit root in the time series


for (i in 1:10) {
  print(i)
  print(testdf(variable = data_xts_in[,i],
  max.augmentations = 3))
}
#-4

data_xts_in_diff <- diff.xts(data_xts_in)

for (i in 1:10) {
  print(i)
  print(testdf(variable = data_xts_in_diff[,i],
               max.augmentations = 3))
}
#1 3 4 9

chosen_ones <- xts(order.by = data$X)

chosen_ones$first<- data_xts_in[,3]
chosen_ones$second<- data_xts_in[,9]

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


cointegration <- lm(first ~ second, data = chosen_ones)

summary(cointegration)

testdf(variable = residuals(cointegration), max.augmentations = 3)
durbinWatsonTest(cointegration)




pp.test(residuals(cointegration))
kpss.test(residuals(cointegration))



#WE NEED A VAR LAG ORDER FIRST
johan.test.trace <- 
  ca.jo(chosen_ones, # data 
        ecdet = "trend", # "none" for no intercept in cointegrating equation, 
        # "const" for constant term in cointegrating equation and 
        # "trend" for trend variable in cointegrating equation
        type = "trace",  # type of the test: trace or eigen
        K = 2,           # lag order of the series (levels) in the VAR
        season = 12) 
summary(johan.test.trace)
?ca.jo

