plot_ACF_PACF_resids <- function(ARIMA_model) {

par(mfrow = c(2, 1)) 
   acf(resid(ARIMA_model), 
       lag.max = 36,
       ylim = c(-0.1, 0.1), 
       lwd = 5, col = "dark green",
       na.action = na.pass)
   pacf(resid(ARIMA_model), 
        lag.max = 36, 
        lwd = 5, col = "dark green",
        na.action = na.pass)
par(mfrow = c(1, 1))
}